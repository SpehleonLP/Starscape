#include "gpu/gpu_march.h"
#include "gpu/device.h"
#include "gpu/pack.h"
#include "march.spv.h"
#include <algorithm>
#include <array>
#include <cstring>
#include <format>
#include <span>

namespace starscape
{
namespace
{

using gpu::Buffer;
using gpu::Device;
using gpu::ScalarWriter;

constexpr uint32_t kGroupSize = 8;       // march.comp's local_size_x and local_size_y
constexpr size_t   kMaxGalaxies = 16;    // march.comp's kMaxGalaxies: the size of its per-galaxy arrays
constexpr uint32_t kSumsPerTexel = 3 + uint32_t(kBins);
constexpr int      kFirstTileRows = 8;
// march.comp's Pass block in scalar layout: mat3, ten uints, six floats. Vulkan guarantees only 128 bytes of
// push constants on every device.
constexpr uint32_t kPassBytes = 9 * 4 + 10 * 4 + 6 * 4;
static_assert(kPassBytes < 128);

// Every earlier submission's writes, made visible to everything this one does. Submit waits on a fence between
// submissions, but a fence orders execution only; the transfer and shader writes still need a memory dependency.
void WaitForEarlierWrites(VkCommandBuffer commands)
{
	VkMemoryBarrier const barrier{.sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER, .srcAccessMask = VK_ACCESS_MEMORY_WRITE_BIT,
		.dstAccessMask = VK_ACCESS_MEMORY_READ_BIT | VK_ACCESS_MEMORY_WRITE_BIT};
	vkCmdPipelineBarrier(commands, VK_PIPELINE_STAGE_ALL_COMMANDS_BIT, VK_PIPELINE_STAGE_ALL_COMMANDS_BIT, 0, 1, &barrier, 0, nullptr, 0,
		nullptr);
}

// A copy's writes, made visible to the host that maps the destination.
void ReleaseToHost(VkCommandBuffer commands)
{
	VkMemoryBarrier const barrier{.sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER, .srcAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT,
		.dstAccessMask = VK_ACCESS_HOST_READ_BIT};
	vkCmdPipelineBarrier(commands, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_PIPELINE_STAGE_HOST_BIT, 0, 1, &barrier, 0, nullptr, 0, nullptr);
}

// A buffer that lives only for one upload, freed on every path out.
class Scratch
{
public:
	explicit Scratch(Device& device) : m_device(device) {}
	~Scratch() { m_device.Destroy(buffer); }
	Scratch(Scratch const&) = delete;
	Scratch& operator=(Scratch const&) = delete;

	Buffer buffer;

private:
	Device& m_device;
};

class GpuMarcher final : public PassMarcher
{
public:
	GpuMarcher(std::unique_ptr<Device>&& device, PreparedScene const& prepared, MarchParams const& params, uint32_t faceSize,
		std::function<void(double)>&& progress)
	: m_device(std::move(device))
	, m_skyToScene(prepared.skyToScene)
	, m_sMin(prepared.sMin)
	, m_params(params)
	, m_faceSize(faceSize)
	, m_galaxyCount(uint32_t(prepared.galaxies.size()))
	, m_progress(std::move(progress))
	{
	}

	~GpuMarcher() override
	{
		m_device->DestroyPipeline(m_pipeline);
		m_device->DestroyDescriptorPool(m_pool);
		m_device->DestroyPipelineLayout(m_pipelineLayout);
		m_device->DestroyDescriptorSetLayout(m_setLayout);
		for (Buffer& buffer : m_buffers)
		{
			m_device->Destroy(buffer);
		}
		m_device->Destroy(m_staging);
	}
	GpuMarcher(GpuMarcher const&) = delete;
	GpuMarcher& operator=(GpuMarcher const&) = delete;

	// Uploads the scene once and builds the pipeline; whatever was made before a failure is freed by the destructor.
	Expected<void> Build(PreparedScene const& prepared)
	{
		std::array<ScalarWriter, 3> tables;
		for (PreparedGalaxy const& galaxy : prepared.galaxies)
		{
			tables[0].Align(gpu::kGalaxyAlignment);
			gpu::PackGalaxy(tables[0], galaxy.galaxy);
			// The observer in galaxy-local coordinates, in double: the one rounding to float happens on the way out.
			tables[1].Mat3(galaxy.toLocal);
			tables[1].Vec3(galaxy.toLocal * (prepared.origin - galaxy.position));
			tables[1].Float(galaxy.boundRadius);
		}
		for (PopulationTable const& population : prepared.populations)
		{
			gpu::PackPopulation(tables[2], population);
		}

		// Zero-size buffers are invalid Vulkan; a scene without galaxies still binds something the shader never reads.
		std::array<VkDeviceSize, 3> sizes{};
		std::array<VkDeviceSize, 3> offsets{};
		VkDeviceSize total = 0;
		for (size_t i = 0; i < tables.size(); ++i)
		{
			sizes[i] = std::max<VkDeviceSize>(tables[i].Bytes().size(), 16);
			offsets[i] = total;
			total += sizes[i];
		}
		Scratch upload(*m_device);
		STARSCAPE_TRY_ASSIGN(staging, m_device->CreateBuffer(total, VK_BUFFER_USAGE_TRANSFER_SRC_BIT, true));
		upload.buffer = staging;
		std::memset(upload.buffer.mapped, 0, total);
		for (size_t i = 0; i < tables.size(); ++i)
		{
			STARSCAPE_TRY_ASSIGN(buffer, m_device->CreateBuffer(sizes[i], VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT, false));
			m_buffers[i] = buffer;
			if (!tables[i].Bytes().empty())
			{
				std::memcpy(static_cast<uint8_t*>(upload.buffer.mapped) + offsets[i], tables[i].Bytes().data(), tables[i].Bytes().size());
			}
		}

		VkDeviceSize const sumBytes = VkDeviceSize(m_faceSize) * m_faceSize * kSumsPerTexel * sizeof(float);
		STARSCAPE_TRY_ASSIGN(sums, m_device->CreateBuffer(sumBytes,
			VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT, false));
		m_buffers[3] = sums;
		STARSCAPE_TRY_ASSIGN(sumStaging, m_device->CreateBuffer(sumBytes, VK_BUFFER_USAGE_TRANSFER_SRC_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT, true));
		m_staging = sumStaging;

		STARSCAPE_TRY(m_device->Submit([&](VkCommandBuffer commands)
		{
			for (size_t i = 0; i < tables.size(); ++i)
			{
				VkBufferCopy const region{offsets[i], 0, sizes[i]};
				vkCmdCopyBuffer(commands, upload.buffer.buffer, m_buffers[i].buffer, 1, &region);
			}
		}));

		STARSCAPE_TRY_ASSIGN(setLayout, m_device->CreateDescriptorSetLayout(uint32_t(m_buffers.size())));
		m_setLayout = setLayout;
		STARSCAPE_TRY_ASSIGN(pipelineLayout, m_device->CreatePipelineLayout(m_setLayout, kPassBytes));
		m_pipelineLayout = pipelineLayout;
		STARSCAPE_TRY_ASSIGN(pool, m_device->CreateDescriptorPool(1, uint32_t(m_buffers.size())));
		m_pool = pool;
		STARSCAPE_TRY_ASSIGN(set, m_device->AllocateDescriptorSet(m_pool, m_setLayout));
		m_set = set;
		m_device->WriteStorageBuffers(m_set, m_buffers);
		STARSCAPE_TRY_ASSIGN(pipeline, m_device->CreateComputePipeline(k_march_spirv, m_pipelineLayout));
		m_pipeline = pipeline;
		return {};
	}

	Expected<void> MarchPass(PassSpec const& spec, PassSums& sums) override
	{
		if (spec.faceSize != m_faceSize)
		{
			return InputError(std::format("the GPU march was built for face size {}, not {}", m_faceSize, spec.faceSize));
		}
		uint32_t const size = m_faceSize;
		uint32_t const faceTexels = size * size;
		float* const staged = static_cast<float*>(m_staging.mapped);
		VkBufferCopy const whole{0, 0, m_staging.size};
		for (uint32_t face = 0; face < 6; ++face)
		{
			for (uint32_t i = 0; i < faceTexels; ++i)
			{
				uint32_t const texel = face * faceTexels + i;
				float* const out = staged + size_t(i) * kSumsPerTexel;
				out[0] = float(sums.radiance[texel].x);
				out[1] = float(sums.radiance[texel].y);
				out[2] = float(sums.radiance[texel].z);
				std::memcpy(out + 3, sums.budget[texel].data(), sizeof(float) * kBins);
			}
			STARSCAPE_TRY(m_device->Submit([&](VkCommandBuffer commands)
			{
				WaitForEarlierWrites(commands);
				vkCmdCopyBuffer(commands, m_staging.buffer, m_buffers[3].buffer, 1, &whole);
			}));

			for (uint32_t rowBegin = 0; rowBegin < size;)
			{
				uint32_t const rowEnd = std::min(rowBegin + uint32_t(m_rows), size);
				STARSCAPE_TRY_ASSIGN(constants, PassConstants(spec, face, rowBegin, rowEnd));
				STARSCAPE_TRY_ASSIGN(seconds, m_device->Submit([&](VkCommandBuffer commands)
				{
					WaitForEarlierWrites(commands);
					vkCmdBindPipeline(commands, VK_PIPELINE_BIND_POINT_COMPUTE, m_pipeline);
					vkCmdBindDescriptorSets(commands, VK_PIPELINE_BIND_POINT_COMPUTE, m_pipelineLayout, 0, 1, &m_set, 0, nullptr);
					vkCmdPushConstants(commands, m_pipelineLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, kPassBytes, constants.data());
					vkCmdDispatch(commands, (size + kGroupSize - 1) / kGroupSize, (rowEnd - rowBegin + kGroupSize - 1) / kGroupSize, 1);
				}));
				m_rows = NextTileRows(m_rows, seconds, int(size));
				if (m_progress)
				{
					m_progress(double(face * size + rowEnd) / double(6 * size));
				}
				rowBegin = rowEnd;
			}

			STARSCAPE_TRY(m_device->Submit([&](VkCommandBuffer commands)
			{
				WaitForEarlierWrites(commands);
				vkCmdCopyBuffer(commands, m_buffers[3].buffer, m_staging.buffer, 1, &whole);
				ReleaseToHost(commands);
			}));
			for (uint32_t i = 0; i < faceTexels; ++i)
			{
				uint32_t const texel = face * faceTexels + i;
				float const* const in = staged + size_t(i) * kSumsPerTexel;
				sums.radiance[texel] = glm::dvec3(in[0], in[1], in[2]);
				std::memcpy(sums.budget[texel].data(), in + 3, sizeof(float) * kBins);
			}
		}
		return {};
	}

	std::string Name() const override
	{
		// The spec asks for LOW queue priority, which some drivers do not offer; the log says which one ran.
		return std::format("GPU {}, driver {}, {} queue priority", m_device->Name(), m_device->Driver(),
			m_device->LowPriority() ? "low" : "default");
	}

private:
	// march.comp's Pass block, member by member.
	Expected<std::vector<uint8_t>> PassConstants(PassSpec const& spec, uint32_t face, uint32_t rowBegin, uint32_t rowEnd) const
	{
		ScalarWriter w;
		w.Mat3(m_skyToScene);
		w.Uint(m_faceSize);
		w.Uint(face);
		w.Uint(rowBegin);
		w.Uint(rowEnd);
		w.Uint(uint32_t(spec.stratum.x));
		w.Uint(uint32_t(spec.stratum.y));
		w.Uint(uint32_t(spec.perAxis));
		w.Uint(m_galaxyCount);
		w.Uint(uint32_t(spec.seed));
		w.Uint(uint32_t(spec.seed >> 32));
		w.Float(spec.fluxThreshold);
		w.Float(1.0 / double(spec.perAxis * spec.perAxis));
		w.Float(m_params.alpha);
		w.Float(m_params.dsMin);
		w.Float(m_params.dsMax);
		w.Float(m_sMin);
		if (w.Bytes().size() != kPassBytes)
		{
			return InputError(std::format("the GPU march's push constants are {} bytes, not the shader's {}", w.Bytes().size(), kPassBytes));
		}
		return w.Bytes();
	}

	std::unique_ptr<Device>     m_device;
	glm::dmat3                  m_skyToScene;
	double                      m_sMin;
	MarchParams                 m_params;
	uint32_t                    m_faceSize;
	uint32_t                    m_galaxyCount;
	std::function<void(double)> m_progress;
	// Galaxies, frames, populations, face sums: march.comp's bindings 0 to 3.
	std::array<Buffer, 4>       m_buffers{};
	Buffer                      m_staging;
	VkDescriptorSetLayout       m_setLayout = VK_NULL_HANDLE;
	VkPipelineLayout            m_pipelineLayout = VK_NULL_HANDLE;
	VkDescriptorPool            m_pool = VK_NULL_HANDLE;
	VkDescriptorSet             m_set = VK_NULL_HANDLE;
	VkPipeline                  m_pipeline = VK_NULL_HANDLE;
	// Carried across faces and passes, so each pass starts from the size the last one settled on.
	int                         m_rows = kFirstTileRows;
};

}

Expected<std::unique_ptr<PassMarcher>> OpenGpuMarcher(PreparedScene const& prepared, MarchParams const& params, uint32_t faceSize,
	std::function<void(double)> progress)
{
	STARSCAPE_TRY_ASSIGN(device, gpu::OpenDevice({}));
	if (prepared.galaxies.size() > kMaxGalaxies)
	{
		return InputError(std::format("the GPU march takes at most {} galaxies; this scene has {}", kMaxGalaxies, prepared.galaxies.size()));
	}
	auto marcher = std::make_unique<GpuMarcher>(std::move(device), prepared, params, faceSize, std::move(progress));
	STARSCAPE_TRY(marcher->Build(prepared));
	return marcher;
}

}
