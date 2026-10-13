#include "population.h"
#include "presets.h"
#include <doctest.h>
#include <cmath>
#include <numeric>
#include <glm/geometric.hpp>

using namespace starscape;

TEST_CASE("Luminosity bins are log-spaced with geometric centres")
{
	CHECK(BinEdge(0) == doctest::Approx(1e-4));
	CHECK(BinEdge(kBins) == doctest::Approx(1e6));
	for (int bin = 0; bin < kBins; ++bin)
	{
		CHECK(BinLuminosity(bin) == doctest::Approx(std::sqrt(BinEdge(bin) * BinEdge(bin + 1))));
	}
}

TEST_CASE("Preset populations are normalised and keep their derived mean luminosity")
{
	// Mutation: a mistyped bin fraction reddens this test.
	std::vector<Population> const presets = PresetPopulations();
	REQUIRE(presets.size() == 4);
	CHECK(presets[0].name == "young-disk");
	CHECK(presets[1].name == "old-disk");
	CHECK(presets[2].name == "bulge");
	CHECK(presets[3].name == "halo");
	// Means from the Kroupa-IMF derivation recorded in presets.cpp.
	double const expectedMean[] = {29.40, 0.343, 0.419, 0.419};
	for (size_t i = 0; i < presets.size(); ++i)
	{
		double const sum = std::accumulate(presets[i].fraction.begin(), presets[i].fraction.end(), 0.0);
		CHECK(sum == doctest::Approx(1.0).epsilon(1e-5));
		PopulationTable const table = Tabulate(presets[i]);
		CHECK(Luminance(table.lightBelow[kBins]) == doctest::Approx(expectedMean[i]).epsilon(0.01).scale(0.0));
		CHECK(Luminance(table.lightBelow[0]) == 0.0);
	}
}

TEST_CASE("Blackbody chroma has unit luminance and the right hue order")
{
	for (double kelvin : {2700.0, 5800.0, 6500.0, 20000.0, 60000.0})
	{
		CHECK(Luminance(BlackbodyChroma(kelvin)) == doctest::Approx(1.0));
	}
	glm::dvec3 const white = BlackbodyChroma(6500.0);
	CHECK(white.x == doctest::Approx(1.0).epsilon(0.2));
	CHECK(white.y == doctest::Approx(1.0).epsilon(0.2));
	CHECK(white.z == doctest::Approx(1.0).epsilon(0.2));
	glm::dvec3 const red = BlackbodyChroma(3000.0);
	CHECK(red.x > red.y);
	CHECK(red.y > red.z);
	glm::dvec3 const blue = BlackbodyChroma(20000.0);
	CHECK(blue.z > blue.y);
	CHECK(blue.y > blue.x);
}

TEST_CASE("Populations are found by name")
{
	std::vector<Population> const presets = PresetPopulations();
	CHECK(FindPopulation(presets, "bulge") == 2);
	CHECK_FALSE(FindPopulation(presets, "nebula").has_value());
}

TEST_CASE("Saturation pushes star colours from white and keeps their light")
{
	// Mutation: skipping the renormalisation in Saturate reddens the luminance check on the clipped bins.
	double const s = 4.0;
	int clipped = 0;
	for (Population const& population : PresetPopulations())
	{
		PopulationTable const physical = Tabulate(population);
		PopulationTable const vivid = Tabulate(population, s);
		for (int bin = 0; bin < kBins; ++bin)
		{
			size_t const i = size_t(bin);
			CHECK(Luminance(vivid.chroma[i]) == doctest::Approx(1.0).epsilon(1e-12));
			glm::dvec3 const pushed = glm::dvec3(1.0) + s * (physical.chroma[i] - glm::dvec3(1.0));
			if (pushed.x >= 0.0 && pushed.y >= 0.0 && pushed.z >= 0.0)
			{
				CHECK(glm::length(vivid.chroma[i] - glm::dvec3(1.0)) == doctest::Approx(s * glm::length(physical.chroma[i] - glm::dvec3(1.0))).epsilon(1e-9));
			}
			else
			{
				++clipped;
				CHECK(vivid.chroma[i].x >= 0.0);
				CHECK(vivid.chroma[i].y >= 0.0);
				CHECK(vivid.chroma[i].z >= 0.0);
			}
			CHECK(Luminance(vivid.lightBelow[i + 1]) == doctest::Approx(Luminance(physical.lightBelow[i + 1])).epsilon(1e-12));
		}
	}
	CHECK(clipped > 0);   // the renormalisation is exercised, not just the in-gamut identity
}
