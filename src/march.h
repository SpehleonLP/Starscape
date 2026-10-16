#pragma once
#include "galaxy.h"
#include "population.h"
#include "field/step.glsl.h"
#include "scene.h"
#include <array>
#include <optional>
#include <vector>
#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>

namespace starscape
{

// Step length ds = clamp(alpha × s, dsMin, dsMax); samplesPerAxis² jittered rays per texel.
struct MarchParams
{
	int    samplesPerAxis = 4;
	double alpha          = 0.02;
	double dsMin          = 0.01;
	double dsMax          = 5.0;
};

MarchParams ParamsFor(Quality quality);

struct PreparedScene
{
	std::vector<PreparedGalaxy>  galaxies;
	std::vector<PopulationTable> populations;
	glm::dvec3                   origin{0.0};
	glm::dmat3                   skyToScene{1.0};
	double                       sMin = 0.05;
};

PreparedScene PrepareScene(Scene const& scene);

struct Interval
{
	double begin;
	double end;
};

// Sorted, merged stretches of the ray inside any galaxy's bounding sphere, clipped to [sMin, maxDistance].
void RayIntervals(PreparedScene const& scene, glm::dvec3 dirScene, double maxDistance, std::vector<Interval>& out);

glm::dvec3 TransmittanceTo(PreparedScene const& scene, MarchParams const& params, glm::dvec3 dirSky, double distance);

// Adds `weight` × this ray's diffuse radiance and its point budget over `sampleSolidAngle` to `out`.
void MarchRay(PreparedScene const& scene, MarchParams const& params, glm::dvec3 dirSky, double fluxThreshold,
	double sampleSolidAngle, double weight, RayResult& out);

// Expected star counts by log10 of apparent flux, for choosing the threshold that meets a star budget.
struct FluxHistogram
{
	static constexpr int    kBuckets = 512;
	static constexpr double kLogMin  = -20.0;
	static constexpr double kLogMax  = 10.0;
	static constexpr double kBucketWidth = (kLogMax - kLogMin) / kBuckets;

	std::array<double, kBuckets> count{};

	void Add(double flux, double stars);
	void Merge(FluxHistogram const& other);
	double CountAbove(double flux) const;
	// Smallest flux whose CountAbove is `target`; nullopt when the whole histogram holds fewer.
	std::optional<double> FluxForCount(double target) const;
};

void HistogramRay(PreparedScene const& scene, MarchParams const& params, glm::dvec3 dirSky, double sampleSolidAngle,
	FluxHistogram& out);

struct StarSample
{
	double     distance;
	glm::dvec3 transmittance;
	int        population;
};

// Re-marches `dirSky` and inverts the cumulative count of resolved stars of `bin` at quantile `u`.
// nullopt when nothing of that bin is resolved along the ray.
std::optional<StarSample> SampleStar(PreparedScene const& scene, MarchParams const& params, glm::dvec3 dirSky, int bin,
	double fluxThreshold, double u);

}
