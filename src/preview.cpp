#include "preview.h"
#include "parallel.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <string>
#include <glm/geometric.hpp>

namespace starscape
{
namespace
{

constexpr double kPi = 3.14159265358979323846;
// FWHM 1.5 px: wide enough that a star never aliases, narrow enough to stay a point.
constexpr double kSplatSigma = 1.5 / 2.3548200450309493;

glm::dvec3 Fetch(kre_starscape::View const& view, uint32_t face, int x, int y)
{
	int const size = int(view.header.diffuseFaceSize);
	if (x < 0 || y < 0 || x >= size || y >= size)
	{
		// Off this face: aim through the would-be texel centre and let the face table pick the neighbour.
		std::array<double, 3> const d = kre_starscape::FaceToDir(face, (x + 0.5) / size, (y + 0.5) / size);
		kre_starscape::FaceCoord const c = kre_starscape::DirToFace(d[0], d[1], d[2]);
		face = c.face;
		x = std::clamp(int(c.s * size), 0, size - 1);
		y = std::clamp(int(c.t * size), 0, size - 1);
	}
	auto const v = kre_starscape::DecodeRGB9E5(kre_starscape::DiffuseTexel(view, 0, face, uint32_t(x), uint32_t(y)));
	return glm::dvec3(v[0], v[1], v[2]) / double(view.header.valueScale);
}

void SplatStars(kre_starscape::View const& view, double exposure, Image& image)
{
	int const width = int(image.width);
	int const height = int(image.height);
	for (uint32_t i = 0; i < view.header.starCount; ++i)
	{
		kre_starscape::Star const star = kre_starscape::StarAt(view, i);
		auto const packed = kre_starscape::DecodeRGB9E5(star.fluxRGB9E5);
		glm::dvec3 const flux = glm::dvec3(packed[0], packed[1], packed[2]) * (exposure / double(view.header.valueScale));
		glm::dvec3 const dir = glm::normalize(glm::dvec3(star.dir[0], star.dir[1], star.dir[2]));
		double const longitude = std::atan2(dir.y, dir.x);
		double const latitude = std::asin(std::clamp(dir.z, -1.0, 1.0));
		// Continuous coordinates where integer + 0.5 is a pixel centre.
		double const u = (longitude + kPi) / (2.0 * kPi) * width;
		double const v = (0.5 * kPi - latitude) / kPi * height;

		struct Tap
		{
			int    column;
			int    row;
			double weight;
		};
		std::array<Tap, 49> taps{};
		int count = 0;
		double weightSum = 0.0;
		int const column0 = int(std::floor(u));
		int const row0 = int(std::floor(v));
		for (int dy = -3; dy <= 3; ++dy)
		{
			int const row = row0 + dy;
			if (row < 0 || row >= height)
			{
				continue;
			}
			for (int dx = -3; dx <= 3; ++dx)
			{
				int const column = column0 + dx;
				double const ex = column + 0.5 - u;
				double const ey = row + 0.5 - v;
				double const weight = std::exp(-(ex * ex + ey * ey) / (2.0 * kSplatSigma * kSplatSigma));
				taps[size_t(count++)] = {((column % width) + width) % width, row, weight};
				weightSum += weight;
			}
		}
		// Normalising by the taps actually kept and dividing by each pixel's solid angle makes
		// the star's integrated flux exact at any latitude, poles included.
		for (int t = 0; t < count; ++t)
		{
			Tap const& tap = taps[size_t(t)];
			glm::dvec3 const radiance = flux * (tap.weight / weightSum / PixelSolidAngle(image.width, uint32_t(tap.row)));
			size_t const index = (size_t(tap.row) * image.width + size_t(tap.column)) * 3;
			image.rgb[index]     += float(radiance.x);
			image.rgb[index + 1] += float(radiance.y);
			image.rgb[index + 2] += float(radiance.z);
		}
	}
}

template<class Radiance>
Image Render(uint32_t width, unsigned threads, Radiance&& radiance)
{
	Image image{width, width / 2, {}};
	image.rgb.assign(size_t(image.width) * image.height * 3, 0.0f);
	ParallelFor(image.height, threads, [&](uint32_t row)
	{
		for (uint32_t column = 0; column < image.width; ++column)
		{
			glm::dvec3 const value = radiance(PixelDirection(image.width, column + 0.5, row + 0.5));
			size_t const index = (size_t(row) * image.width + column) * 3;
			image.rgb[index]     = float(value.x);
			image.rgb[index + 1] = float(value.y);
			image.rgb[index + 2] = float(value.z);
		}
	}, {});
	return image;
}

}

double PixelSolidAngle(uint32_t width, uint32_t row)
{
	uint32_t const height = width / 2;
	double const top = 0.5 * kPi - kPi * row / height;
	double const bottom = 0.5 * kPi - kPi * (row + 1) / height;
	return 2.0 * kPi / width * (std::sin(top) - std::sin(bottom));
}

glm::dvec3 PixelDirection(uint32_t width, double column, double row)
{
	uint32_t const height = width / 2;
	double const longitude = 2.0 * kPi * column / width - kPi;
	double const latitude = 0.5 * kPi - kPi * row / height;
	return {std::cos(latitude) * std::cos(longitude), std::cos(latitude) * std::sin(longitude), std::sin(latitude)};
}

glm::dvec3 SampleDiffuse(kre_starscape::View const& view, glm::dvec3 dir)
{
	int const size = int(view.header.diffuseFaceSize);
	kre_starscape::FaceCoord const c = kre_starscape::DirToFace(dir.x, dir.y, dir.z);
	double const px = c.s * size - 0.5;
	double const py = c.t * size - 0.5;
	int const x0 = int(std::floor(px));
	int const y0 = int(std::floor(py));
	double const fx = px - x0;
	double const fy = py - y0;
	return (1.0 - fx) * (1.0 - fy) * Fetch(view, c.face, x0, y0)
		+ fx * (1.0 - fy) * Fetch(view, c.face, x0 + 1, y0)
		+ (1.0 - fx) * fy * Fetch(view, c.face, x0, y0 + 1)
		+ fx * fy * Fetch(view, c.face, x0 + 1, y0 + 1);
}

Image RenderPreview(kre_starscape::View const& view, uint32_t width, double exposureEV, unsigned threads)
{
	double const exposure = std::exp2(exposureEV);
	Image image = Render(width, threads, [&](glm::dvec3 dir) { return exposure * SampleDiffuse(view, dir); });
	SplatStars(view, exposure, image);
	return image;
}

Image RenderMarchedPreview(PreparedScene const& scene, MarchParams const& params, kre_starscape::View const& view,
	uint32_t width, double exposureEV, unsigned threads)
{
	double const exposure = std::exp2(exposureEV);
	double const threshold = double(view.header.fluxThreshold);
	Image image = Render(width, threads, [&](glm::dvec3 dir)
	{
		RayResult result;
		MarchRay(scene, params, dir, threshold, 0.0, 1.0, result);
		return exposure * result.radiance;
	});
	SplatStars(view, exposure, image);
	return image;
}

std::vector<uint8_t> EncodeRadiance(Image const& image)
{
	std::string const header = std::format("#?RADIANCE\nFORMAT=32-bit_rle_rgbe\n\n-Y {} +X {}\n", image.height, image.width);
	std::vector<uint8_t> bytes(header.begin(), header.end());
	bytes.reserve(bytes.size() + size_t(image.width) * image.height * 4);
	for (size_t i = 0; i + 2 < image.rgb.size(); i += 3)
	{
		float const r = std::max(image.rgb[i], 0.0f);
		float const g = std::max(image.rgb[i + 1], 0.0f);
		float const b = std::max(image.rgb[i + 2], 0.0f);
		float const v = std::max({r, g, b});
		if (!(v >= 1e-32f))
		{
			bytes.insert(bytes.end(), {0, 0, 0, 0});
			continue;
		}
		int exponent = 0;
		float const mantissa = std::frexp(v, &exponent);
		float const scale = mantissa * 256.0f / v;
		bytes.insert(bytes.end(), {uint8_t(r * scale), uint8_t(g * scale), uint8_t(b * scale), uint8_t(exponent + 128)});
	}
	return bytes;
}

}
