#include "writer.h"
#include <kre_starscape.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <format>
#include <fstream>
#include <system_error>

namespace starscape
{

std::vector<std::vector<glm::dvec3>> BuildMips(uint32_t faceSize, std::vector<glm::dvec3> const& level0)
{
	std::vector<std::vector<glm::dvec3>> levels{level0};
	for (uint32_t size = faceSize; size > 1; size /= 2)
	{
		uint32_t const half = size / 2;
		std::vector<glm::dvec3> coarse(size_t(6) * half * half);
		std::vector<glm::dvec3> const& fine = levels.back();
		for (uint32_t face = 0; face < 6; ++face)
		{
			for (uint32_t y = 0; y < half; ++y)
			{
				for (uint32_t x = 0; x < half; ++x)
				{
					glm::dvec3 sum(0.0);
					double weight = 0.0;
					for (uint32_t dy = 0; dy < 2; ++dy)
					{
						for (uint32_t dx = 0; dx < 2; ++dx)
						{
							uint32_t const fx = 2 * x + dx;
							uint32_t const fy = 2 * y + dy;
							double const w = kre_starscape::TexelSolidAngle(size, fx, fy);
							sum += w * fine[(size_t(face) * size + fy) * size + fx];
							weight += w;
						}
					}
					coarse[(size_t(face) * half + y) * half + x] = sum / weight;
				}
			}
		}
		levels.push_back(std::move(coarse));
	}
	return levels;
}

float ChooseValueScale(double maxValue)
{
	if (!(maxValue > 0.0) || !std::isfinite(maxValue))
	{
		return 1.0f;
	}
	// Clamped while still a double: a denormal maximum makes the quotient infinite, and
	// converting an out-of-range double to int is undefined.
	double const exponent = std::clamp(std::floor(std::log2(32768.0 / maxValue)), -100.0, 100.0);
	return std::ldexp(1.0f, int(exponent));
}

std::vector<uint8_t> BuildFile(Cooked const& cooked, uint64_t sceneHash)
{
	std::vector<std::vector<glm::dvec3>> const levels = BuildMips(cooked.faceSize, cooked.diffuse);
	double maxValue = 0.0;
	// Non-finite values are skipped explicitly: the encoder zeroes NaN and clamps infinity,
	// so they must not set the scale (std::max's NaN behaviour depends on argument order).
	auto accumulate = [&maxValue](glm::dvec3 const& v)
	{
		for (int i = 0; i < 3; ++i)
		{
			if (std::isfinite(v[i]) && v[i] > maxValue)
			{
				maxValue = v[i];
			}
		}
	};
	for (glm::dvec3 const& v : cooked.diffuse)
	{
		accumulate(v);
	}
	for (CookedStar const& star : cooked.stars)
	{
		accumulate(star.flux);
	}
	float const scale = ChooseValueScale(maxValue);

	kre_starscape::Header header{};
	header.magic           = kre_starscape::kMagic;
	header.version         = kre_starscape::kVersion;
	header.sceneHash       = sceneHash;
	header.diffuseEncoding = uint32_t(kre_starscape::Encoding::RGB9E5);
	header.diffuseFaceSize = cooked.faceSize;
	header.diffuseMipCount = kre_starscape::MipCountFor(cooked.faceSize);
	header.starCount       = uint32_t(cooked.stars.size());
	header.diffuseOffset   = sizeof(kre_starscape::Header);
	header.diffuseBytes    = kre_starscape::DiffuseBytesFor(cooked.faceSize);
	header.starsOffset     = header.diffuseOffset + header.diffuseBytes;
	header.fluxThreshold   = float(cooked.fluxThreshold);
	header.valueScale      = scale;

	std::vector<uint8_t> bytes(size_t(header.starsOffset + uint64_t(header.starCount) * sizeof(kre_starscape::Star)));
	std::memcpy(bytes.data(), &header, sizeof header);
	uint8_t* cursor = bytes.data() + header.diffuseOffset;
	auto pack = [scale](glm::dvec3 const& v)
	{
		return kre_starscape::EncodeRGB9E5(float(v.x * scale), float(v.y * scale), float(v.z * scale));
	};
	for (std::vector<glm::dvec3> const& level : levels)
	{
		for (glm::dvec3 const& texel : level)
		{
			uint32_t const packed = pack(texel);
			std::memcpy(cursor, &packed, 4);
			cursor += 4;
		}
	}
	for (CookedStar const& star : cooked.stars)
	{
		kre_starscape::Star const record{{float(star.dir.x), float(star.dir.y), float(star.dir.z)}, pack(star.flux)};
		std::memcpy(cursor, &record, sizeof record);
		cursor += sizeof record;
	}
	return bytes;
}

Expected<void> WriteFileAtomic(std::filesystem::path const& path, std::span<uint8_t const> bytes)
{
	std::filesystem::path partial = path;
	partial += ".partial";
	{
		std::ofstream file(partial, std::ios::binary | std::ios::trunc);
		if (!file)
		{
			return IoError(std::format("{}: cannot create", partial.string()));
		}
		file.write(reinterpret_cast<char const*>(bytes.data()), std::streamsize(bytes.size()));
		file.close();
		if (!file)
		{
			std::error_code ignored;
			std::filesystem::remove(partial, ignored);
			return IoError(std::format("{}: write failed", partial.string()));
		}
	}
	std::error_code error;
	std::filesystem::rename(partial, path, error);
	if (error)
	{
		std::error_code ignored;
		std::filesystem::remove(partial, ignored);
		return IoError(std::format("{}: {}", path.string(), error.message()));
	}
	return {};
}

Expected<std::vector<uint8_t>> ReadFileBytes(std::filesystem::path const& path)
{
	// ifstream opens a directory successfully on Linux and tellg then returns garbage.
	std::error_code statError;
	if (!std::filesystem::is_regular_file(path, statError))
	{
		return IoError(std::format("{}: not a regular file", path.string()));
	}
	std::ifstream file(path, std::ios::binary | std::ios::ate);
	if (!file)
	{
		return IoError(std::format("{}: cannot open", path.string()));
	}
	std::streamsize const size = file.tellg();
	file.seekg(0);
	std::vector<uint8_t> bytes(size_t(std::max<std::streamsize>(size, 0)));
	file.read(reinterpret_cast<char*>(bytes.data()), size);
	if (!file)
	{
		return IoError(std::format("{}: read failed", path.string()));
	}
	return bytes;
}

}
