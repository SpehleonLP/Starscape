#pragma once
// kre_starscape: the cooked night sky written by Starscape (https://github.com/SpehleonLP/Starscape).
// Header-only and standard-library-only so engines can include it as-is. Layout: docs/spec.md, "File format".
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <span>
#include <string_view>

namespace kre_starscape
{

inline constexpr uint32_t kMagic   = 0x5254534Bu;   // "KSTR" read as little-endian bytes
inline constexpr uint32_t kVersion = 1;
inline constexpr uint32_t kMaxFaceSize = 4096;

enum class Encoding : uint32_t
{
	RGB9E5 = 1,
	BC6H   = 2,   // reserved; no writer or reader yet
};

// Read and written with one memcpy; the static_asserts pin every offset.
struct Header
{
	uint32_t magic;
	uint32_t version;
	uint64_t sceneHash;        // FNV-1a 64 of the scene file bytes, for provenance
	uint32_t diffuseEncoding;
	uint32_t diffuseFaceSize;  // power of two
	uint32_t diffuseMipCount;  // full chain down to 1×1
	uint32_t starCount;
	uint64_t diffuseOffset;
	uint64_t diffuseBytes;
	uint64_t starsOffset;
	float    fluxThreshold;    // F_res, magnitude-0 flux units (physical, not scaled)
	float    valueScale;       // every stored radiance and flux = physical × valueScale
};
static_assert(sizeof(Header) == 64);
static_assert(offsetof(Header, sceneHash) == 8);
static_assert(offsetof(Header, diffuseOffset) == 32);
static_assert(offsetof(Header, starsOffset) == 48);
static_assert(offsetof(Header, valueScale) == 60);

// Uploads as a vertex buffer without repacking.
struct Star
{
	float    dir[3];       // unit vector, sky frame
	uint32_t fluxRGB9E5;   // linear Rec.709 flux × valueScale
};
static_assert(sizeof(Star) == 16);
static_assert(std::endian::native == std::endian::little, "kre_starscape is little-endian on disk");

// ---- RGB9E5 (VK_FORMAT_E5B9G9R9_UFLOAT_PACK32): 9-bit mantissas, shared 5-bit exponent, bias 15.

inline constexpr float kRGB9E5Max = 65408.0f;   // (511 / 512) × 2^16

inline uint32_t EncodeRGB9E5(float r, float g, float b)
{
	// `c > 0` is false for NaN, so NaN and negatives both clamp to zero, as the reference does.
	auto clampChannel = [](float c) { return c > 0.0f ? std::min(c, kRGB9E5Max) : 0.0f; };
	float const rc = clampChannel(r);
	float const gc = clampChannel(g);
	float const bc = clampChannel(b);
	float const maxChannel = std::max({rc, gc, bc});
	if (maxChannel == 0.0f)
	{
		return 0u;
	}
	// frexp gives maxChannel = m × 2^e with m in [0.5, 1), so floor(log2(maxChannel)) = e − 1 exactly,
	// where the reference's floor(log2()) can round the wrong way next to a power of two.
	int exponent = 0;
	std::frexp(maxChannel, &exponent);
	int sharedExponent = std::max(-16, exponent - 1) + 16;
	double scale = std::ldexp(1.0, sharedExponent - 24);
	if (int(std::floor(maxChannel / scale + 0.5)) == 512)
	{
		sharedExponent += 1;
		scale *= 2.0;
	}
	auto quantize = [scale](float c) { return uint32_t(std::floor(c / scale + 0.5)); };
	return quantize(rc) | (quantize(gc) << 9) | (quantize(bc) << 18) | (uint32_t(sharedExponent) << 27);
}

inline std::array<float, 3> DecodeRGB9E5(uint32_t packed)
{
	float const scale = std::ldexp(1.0f, int(packed >> 27) - 24);
	return {float(packed & 511u) * scale, float((packed >> 9) & 511u) * scale, float((packed >> 18) & 511u) * scale};
}

// ---- Cube map faces, Vulkan order and orientation ("Cube Map Face Selection" in the Vulkan spec).

enum Face : uint32_t
{
	PosX,
	NegX,
	PosY,
	NegY,
	PosZ,
	NegZ,
};

struct FaceCoord
{
	uint32_t face;
	double   s;   // [0, 1], left to right
	double   t;   // [0, 1], top to bottom (row 0 is the first row in memory)
};

inline FaceCoord DirToFace(double x, double y, double z)
{
	double const ax = std::abs(x);
	double const ay = std::abs(y);
	double const az = std::abs(z);
	uint32_t face;
	double sc;
	double tc;
	double ma;
	if (ax >= ay && ax >= az)
	{
		ma = ax;
		face = x >= 0.0 ? PosX : NegX;
		sc = x >= 0.0 ? -z : z;
		tc = -y;
	}
	else if (ay >= az)
	{
		ma = ay;
		face = y >= 0.0 ? PosY : NegY;
		sc = x;
		tc = y >= 0.0 ? z : -z;
	}
	else
	{
		ma = az;
		face = z >= 0.0 ? PosZ : NegZ;
		sc = z >= 0.0 ? x : -x;
		tc = -y;
	}
	return {face, 0.5 * (sc / ma + 1.0), 0.5 * (tc / ma + 1.0)};
}

// Unnormalised direction through (s, t) on `face`; s and t may lie outside [0, 1].
inline std::array<double, 3> FaceToDir(uint32_t face, double s, double t)
{
	double const u = 2.0 * s - 1.0;
	double const v = 2.0 * t - 1.0;
	switch (face)
	{
	case PosX:
		return {1.0, -v, -u};
	case NegX:
		return {-1.0, -v, u};
	case PosY:
		return {u, 1.0, v};
	case NegY:
		return {u, -1.0, -v};
	case PosZ:
		return {u, -v, 1.0};
	default:
		return {-u, -v, -1.0};
	}
}

inline double TexelSolidAngle(uint32_t faceSize, uint32_t x, uint32_t y)
{
	auto area = [](double u, double v) { return std::atan2(u * v, std::sqrt(u * u + v * v + 1.0)); };
	double const step = 2.0 / double(faceSize);
	double const u0 = -1.0 + double(x) * step;
	double const v0 = -1.0 + double(y) * step;
	double const u1 = u0 + step;
	double const v1 = v0 + step;
	return area(u0, v0) - area(u0, v1) - area(u1, v0) + area(u1, v1);
}

inline uint32_t MipCountFor(uint32_t faceSize)
{
	return uint32_t(std::bit_width(faceSize));
}

// Bytes from the start of the diffuse payload to the first texel of `mip`.
inline uint64_t MipOffset(uint32_t faceSize, uint32_t mip)
{
	uint64_t offset = 0;
	for (uint32_t level = 0; level < mip; ++level)
	{
		uint64_t const size = faceSize >> level;
		offset += 6u * size * size * 4u;
	}
	return offset;
}

inline uint64_t DiffuseBytesFor(uint32_t faceSize)
{
	return MipOffset(faceSize, MipCountFor(faceSize));
}

// ---- Parsing

enum class Result
{
	Ok,
	TooSmall,
	BadMagic,
	BadVersion,
	BadEncoding,
	BadFaceSize,
	BadMipCount,
	BadValueScale,
	DiffuseSizeMismatch,
	DiffuseOutOfBounds,
	StarsOutOfBounds,
	RegionsOverlap,
};

inline std::string_view Describe(Result result)
{
	switch (result)
	{
	case Result::Ok:                  return "ok";
	case Result::TooSmall:            return "file is smaller than its header";
	case Result::BadMagic:            return "not a kre_starscape file (bad magic)";
	case Result::BadVersion:          return "unsupported kre_starscape version";
	case Result::BadEncoding:         return "unsupported diffuse encoding";
	case Result::BadFaceSize:         return "diffuse face size is not a power of two in [1, 4096]";
	case Result::BadMipCount:         return "mip count does not match the face size";
	case Result::BadValueScale:       return "value scale is not a positive finite number";
	case Result::DiffuseSizeMismatch: return "diffuse byte count does not match the face size";
	case Result::DiffuseOutOfBounds:  return "diffuse payload lies outside the file";
	case Result::StarsOutOfBounds:    return "star array lies outside the file";
	case Result::RegionsOverlap:      return "diffuse payload and star array overlap";
	}
	return "unknown error";
}

struct View
{
	Header                   header{};
	std::span<uint8_t const> diffuse;   // mip-major, then face, then rows top to bottom
	std::span<uint8_t const> stars;     // header.starCount × Star
};

// Never asserts: every field is checked before anything is read through it.
inline Result Parse(std::span<uint8_t const> file, View& out)
{
	if (file.size() < sizeof(Header))
	{
		return Result::TooSmall;
	}
	Header h;
	std::memcpy(&h, file.data(), sizeof h);
	if (h.magic != kMagic)
	{
		return Result::BadMagic;
	}
	if (h.version != kVersion)
	{
		return Result::BadVersion;
	}
	if (h.diffuseEncoding != uint32_t(Encoding::RGB9E5))
	{
		return Result::BadEncoding;
	}
	if (h.diffuseFaceSize == 0 || h.diffuseFaceSize > kMaxFaceSize || !std::has_single_bit(h.diffuseFaceSize))
	{
		return Result::BadFaceSize;
	}
	if (h.diffuseMipCount != MipCountFor(h.diffuseFaceSize))
	{
		return Result::BadMipCount;
	}
	if (!(h.valueScale > 0.0f) || !std::isfinite(h.valueScale))
	{
		return Result::BadValueScale;
	}
	if (h.diffuseBytes != DiffuseBytesFor(h.diffuseFaceSize))
	{
		return Result::DiffuseSizeMismatch;
	}
	uint64_t const size = file.size();
	// Written as `length <= size − offset` so a hostile offset cannot overflow the sum.
	auto inBounds = [size](uint64_t offset, uint64_t length)
	{
		return offset >= sizeof(Header) && offset <= size && length <= size - offset;
	};
	if (!inBounds(h.diffuseOffset, h.diffuseBytes))
	{
		return Result::DiffuseOutOfBounds;
	}
	uint64_t const starBytes = uint64_t(h.starCount) * sizeof(Star);
	if (!inBounds(h.starsOffset, starBytes))
	{
		return Result::StarsOutOfBounds;
	}
	bool const disjoint = h.diffuseOffset + h.diffuseBytes <= h.starsOffset || h.starsOffset + starBytes <= h.diffuseOffset;
	if (starBytes != 0 && !disjoint)
	{
		return Result::RegionsOverlap;
	}
	out.header  = h;
	out.diffuse = file.subspan(size_t(h.diffuseOffset), size_t(h.diffuseBytes));
	out.stars   = file.subspan(size_t(h.starsOffset), size_t(starBytes));
	return Result::Ok;
}

// The payload need not be 4-byte aligned in the caller's buffer, so reads go through memcpy.
// Caller contract for DiffuseTexel and StarAt: Parse validates the header and region bounds, not the indices
// passed here. The caller must pass mip < header.diffuseMipCount, face < 6, x and y < (diffuseFaceSize >> mip),
// and index < header.starCount; anything else reads outside the view. They stay unchecked because they sit on
// per-texel loops.
inline uint32_t DiffuseTexel(View const& view, uint32_t mip, uint32_t face, uint32_t x, uint32_t y)
{
	uint32_t const size = view.header.diffuseFaceSize >> mip;
	uint64_t const offset = MipOffset(view.header.diffuseFaceSize, mip) + ((uint64_t(face) * size + y) * size + x) * 4u;
	uint32_t value;
	std::memcpy(&value, view.diffuse.data() + offset, 4);
	return value;
}

inline Star StarAt(View const& view, uint32_t index)
{
	Star star;
	std::memcpy(&star, view.stars.data() + size_t(index) * sizeof(Star), sizeof star);
	return star;
}

}
