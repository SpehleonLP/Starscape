#include "population.h"
#include "presets.h"
#include <doctest.h>
#include <cmath>
#include <numeric>

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
	double const expectedMean[] = {29.42, 1.124, 1.295, 1.295};
	for (size_t i = 0; i < presets.size(); ++i)
	{
		double const sum = std::accumulate(presets[i].fraction.begin(), presets[i].fraction.end(), 0.0);
		CHECK(sum == doctest::Approx(1.0).epsilon(1e-5));
		PopulationTable const table = Tabulate(presets[i]);
		CHECK(Luminance(table.lightBelow[kBins]) == doctest::Approx(expectedMean[i]).epsilon(0.01));
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
