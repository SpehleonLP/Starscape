#include "gpu/probe.h"
#include "gpu/pack.h"
#include "probe.spv.h"
#include <algorithm>
#include <cstring>

namespace starscape::gpu
{
namespace
{

constexpr uint32_t kGroupSize = 64;   // probe.comp's local_size_x
constexpr uint32_t kEvaluateWords = 17;   // count, extinction, gas rgb, 6 x (density, population)

struct PushConstants
{
	uint32_t mode;
	uint32_t count;
	uint64_t seed;
};

// Owns everything one probe dispatch creates, so each early return frees what was made so far.
class Run
{
public:
	explicit Run(Device& device) : m_device(device) {}
	~Run()
	{
		m_device.DestroyPipeline(m_pipeline);
		m_device.DestroyDescriptorPool(m_pool);
		m_device.DestroyPipelineLayout(m_pipelineLayout);
		m_device.DestroyDescriptorSetLayout(m_setLayout);
		for (Buffer& buffer : m_buffers)
		{
			m_device.Destroy(buffer);
		}
	}
	Run(Run const&) = delete;
	Run& operator=(Run const&) = delete;

	// Galaxy bytes, input words and output words; none may be empty (a zero-size buffer is invalid).
	Expected<void> Execute(uint32_t mode, uint32_t count, uint64_t seed, std::span<uint8_t const> galaxy, std::span<uint32_t const> inputs,
		size_t outputWords)
	{
		std::array<size_t, 3> const sizes{std::max<size_t>(galaxy.size(), 16), std::max<size_t>(inputs.size_bytes(), 16), outputWords * sizeof(uint32_t)};
		for (size_t i = 0; i < sizes.size(); ++i)
		{
			STARSCAPE_TRY_ASSIGN(buffer, m_device.CreateBuffer(sizes[i], VK_BUFFER_USAGE_STORAGE_BUFFER_BIT, true));
			m_buffers[i] = buffer;
			std::memset(buffer.mapped, 0, sizes[i]);
		}
		std::memcpy(m_buffers[0].mapped, galaxy.data(), galaxy.size());
		std::memcpy(m_buffers[1].mapped, inputs.data(), inputs.size_bytes());

		STARSCAPE_TRY_ASSIGN(setLayout, m_device.CreateDescriptorSetLayout(3));
		m_setLayout = setLayout;
		STARSCAPE_TRY_ASSIGN(pipelineLayout, m_device.CreatePipelineLayout(m_setLayout, sizeof(PushConstants)));
		m_pipelineLayout = pipelineLayout;
		STARSCAPE_TRY_ASSIGN(pool, m_device.CreateDescriptorPool(1, 3));
		m_pool = pool;
		STARSCAPE_TRY_ASSIGN(set, m_device.AllocateDescriptorSet(m_pool, m_setLayout));
		m_device.WriteStorageBuffers(set, m_buffers);
		STARSCAPE_TRY_ASSIGN(pipeline, m_device.CreateComputePipeline(k_probe_spirv, m_pipelineLayout));
		m_pipeline = pipeline;

		PushConstants const push{mode, count, seed};
		STARSCAPE_TRY(m_device.Submit([&](VkCommandBuffer commands)
		{
			vkCmdBindPipeline(commands, VK_PIPELINE_BIND_POINT_COMPUTE, m_pipeline);
			vkCmdBindDescriptorSets(commands, VK_PIPELINE_BIND_POINT_COMPUTE, m_pipelineLayout, 0, 1, &set, 0, nullptr);
			vkCmdPushConstants(commands, m_pipelineLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof push, &push);
			vkCmdDispatch(commands, (count + kGroupSize - 1) / kGroupSize, 1, 1);
		}));
		return {};
	}

	// Valid after Execute; the memory is coherent, and Submit has waited for the GPU.
	uint32_t const* Output() const { return static_cast<uint32_t const*>(m_buffers[2].mapped); }

private:
	Device&                   m_device;
	std::array<Buffer, 3>     m_buffers{};
	VkDescriptorSetLayout     m_setLayout = VK_NULL_HANDLE;
	VkPipelineLayout          m_pipelineLayout = VK_NULL_HANDLE;
	VkDescriptorPool          m_pool = VK_NULL_HANDLE;
	VkPipeline                m_pipeline = VK_NULL_HANDLE;
};

uint32_t Low(uint64_t value)
{
	return uint32_t(value);
}

uint32_t High(uint64_t value)
{
	return uint32_t(value >> 32);
}

float FloatOf(uint32_t bits)
{
	float value;
	std::memcpy(&value, &bits, sizeof value);
	return value;
}

}

Expected<std::vector<uint32_t>> ProbeLatticeHash(Device& device, std::span<glm::i64vec3 const> cells, uint64_t seed)
{
	std::vector<uint32_t> inputs;
	inputs.reserve(cells.size() * 6);
	for (glm::i64vec3 const& cell : cells)
	{
		for (int axis = 0; axis < 3; ++axis)
		{
			uint64_t const bits = uint64_t(cell[axis]);
			inputs.push_back(Low(bits));
			inputs.push_back(High(bits));
		}
	}
	Run run(device);
	STARSCAPE_TRY(run.Execute(0, uint32_t(cells.size()), seed, {}, inputs, cells.size()));
	return std::vector<uint32_t>(run.Output(), run.Output() + cells.size());
}

Expected<std::vector<std::array<uint32_t, 4>>> ProbePhilox(Device& device, std::span<std::array<uint32_t, 4> const> counters, uint64_t key)
{
	std::vector<uint32_t> inputs;
	inputs.reserve(counters.size() * 4);
	for (std::array<uint32_t, 4> const& counter : counters)
	{
		inputs.insert(inputs.end(), counter.begin(), counter.end());
	}
	Run run(device);
	STARSCAPE_TRY(run.Execute(1, uint32_t(counters.size()), key, {}, inputs, counters.size() * 4));
	std::vector<std::array<uint32_t, 4>> blocks(counters.size());
	std::memcpy(blocks.data(), run.Output(), blocks.size() * sizeof blocks[0]);
	return blocks;
}

Expected<std::vector<GalaxySample>> ProbeEvaluate(Device& device, Galaxy const& galaxy, std::span<glm::dvec3 const> points)
{
	ScalarWriter packed;
	PackGalaxy(packed, galaxy);
	std::vector<uint32_t> inputs;
	inputs.reserve(points.size() * 3);
	for (glm::dvec3 const& point : points)
	{
		for (int axis = 0; axis < 3; ++axis)
		{
			float const narrowed = float(point[axis]);
			uint32_t bits;
			std::memcpy(&bits, &narrowed, sizeof bits);
			inputs.push_back(bits);
		}
	}
	Run run(device);
	STARSCAPE_TRY(run.Execute(2, uint32_t(points.size()), 0, packed.Bytes(), inputs, points.size() * kEvaluateWords));
	std::vector<GalaxySample> samples(points.size());
	for (size_t i = 0; i < points.size(); ++i)
	{
		uint32_t const* words = run.Output() + i * kEvaluateWords;
		GalaxySample& sample = samples[i];
		sample.count = int(words[0]);
		sample.extinction = FloatOf(words[1]);
		sample.gas = glm::dvec3(FloatOf(words[2]), FloatOf(words[3]), FloatOf(words[4]));
		for (int k = 0; k < 6; ++k)
		{
			sample.parts[size_t(k)].density = FloatOf(words[5 + 2 * k]);
			sample.parts[size_t(k)].population = int(words[6 + 2 * k]);
		}
	}
	return samples;
}

}
