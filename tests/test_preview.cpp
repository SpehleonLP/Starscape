#include "preview.h"
#include "writer.h"
#include <kre_starscape.h>
#include <doctest.h>
#include <cmath>
#include <vector>
#include <string>

using namespace starscape;

namespace
{

constexpr double kPi = 3.14159265358979323846;

std::vector<uint8_t> SkyFile(std::vector<glm::dvec3>&& diffuse, uint32_t faceSize, std::vector<CookedStar>&& stars)
{
	Cooked cooked;
	cooked.faceSize = faceSize;
	cooked.diffuse = std::move(diffuse);
	cooked.stars = std::move(stars);
	cooked.fluxThreshold = 0.01;
	return BuildFile(cooked, 0);
}

}

TEST_CASE("A splatted star keeps its exact flux at every latitude")
{
	for (double const latitude : {0.0, 60.0, 85.0})
	{
		double const b = latitude * kPi / 180.0;
		glm::dvec3 const dir(std::cos(b) * std::cos(0.3), std::cos(b) * std::sin(0.3), std::sin(b));
		std::vector<uint8_t> const bytes = SkyFile(std::vector<glm::dvec3>(6, glm::dvec3(0.0)), 1, {{dir, {1.0, 2.0, 3.0}}});
		kre_starscape::View view;
		REQUIRE(kre_starscape::Parse(bytes, view) == kre_starscape::Result::Ok);
		Image const image = RenderPreview(view, 512, 0.0, 1);
		glm::dvec3 total(0.0);
		for (uint32_t row = 0; row < image.height; ++row)
		{
			for (uint32_t column = 0; column < image.width; ++column)
			{
				size_t const i = (size_t(row) * image.width + column) * 3;
				total += glm::dvec3(image.rgb[i], image.rgb[i + 1], image.rgb[i + 2]) * PixelSolidAngle(image.width, row);
			}
		}
		CHECK(total.x == doctest::Approx(1.0).epsilon(0.005));
		CHECK(total.y == doctest::Approx(2.0).epsilon(0.005));
		CHECK(total.z == doctest::Approx(3.0).epsilon(0.005));
	}
}

TEST_CASE("Pixel solid angles cover the sphere")
{
	double total = 0.0;
	for (uint32_t row = 0; row < 128; ++row)
	{
		total += 256.0 * PixelSolidAngle(256, row);
	}
	CHECK(total == doctest::Approx(4.0 * kPi));
	glm::dvec3 const top = PixelDirection(256, 128.0, 0.5);
	CHECK(top.z > 0.999);
}

TEST_CASE("The diffuse layer is sampled without seams")
{
	// Radiance 2 + x + z everywhere (offset keeps it positive): a smooth field with a gradient across every
	// face edge, so any seam mistake shows as an outlier.
	constexpr uint32_t size = 8;
	std::vector<glm::dvec3> diffuse;
	for (uint32_t face = 0; face < 6; ++face)
	{
		for (uint32_t y = 0; y < size; ++y)
		{
			for (uint32_t x = 0; x < size; ++x)
			{
				auto const d = kre_starscape::FaceToDir(face, (x + 0.5) / size, (y + 0.5) / size);
				double const norm = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
				diffuse.push_back(glm::dvec3(2.0 + (d[0] + d[2]) / norm));
			}
		}
	}
	std::vector<uint8_t> const bytes = SkyFile(std::move(diffuse), size, {});
	kre_starscape::View view;
	REQUIRE(kre_starscape::Parse(bytes, view) == kre_starscape::Result::Ok);
	Image const image = RenderPreview(view, 256, 0.0, 4);
	for (uint32_t row = 0; row < image.height; ++row)
	{
		for (uint32_t column = 0; column < image.width; ++column)
		{
			glm::dvec3 const d = PixelDirection(image.width, column + 0.5, row + 0.5);
			REQUIRE(image.rgb[(size_t(row) * image.width + column) * 3] == doctest::Approx(2.0 + d.x + d.z).epsilon(0.03));
		}
	}
}

TEST_CASE("Radiance RGBE encoding")
{
	Image image{2, 1, {1.0f, 0.5f, 0.25f, 0.0f, 0.0f, 0.0f}};
	std::vector<uint8_t> const bytes = EncodeRadiance(image);
	std::string const header = "#?RADIANCE\nFORMAT=32-bit_rle_rgbe\n\n-Y 1 +X 2\n";
	REQUIRE(bytes.size() == header.size() + 8);
	CHECK(std::string(bytes.begin(), bytes.begin() + std::ptrdiff_t(header.size())) == header);
	CHECK(bytes[header.size() + 0] == 128);
	CHECK(bytes[header.size() + 1] == 64);
	CHECK(bytes[header.size() + 2] == 32);
	CHECK(bytes[header.size() + 3] == 129);
	CHECK(bytes[header.size() + 7] == 0);
}
