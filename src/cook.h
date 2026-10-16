#pragma once
#include "common.h"
#include "march.h"
#include "scene.h"
#include <array>
#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace starscape
{

class RandomStream;

struct CookedStar
{
	glm::dvec3 dir;    // sky frame, unit
	glm::dvec3 flux;   // mag-0 units, linear Rec.709, physical
};

struct CookStats
{
	double   pointFlux = 0.0;      // Σ Luminance(star flux)
	double   diffuseFlux = 0.0;    // Σ Luminance(radiance) × texel solid angle
	uint64_t starCount = 0;
	uint64_t droppedStars = 0;     // drawn, but their direction held no resolved stars on re-march
	bool     budgetUnreachable = false;
};

struct Cooked
{
	uint32_t                faceSize = 0;
	std::vector<glm::dvec3> diffuse;   // mip 0, texel order, physical radiance
	std::vector<CookedStar> stars;     // brightest first
	double                  fluxThreshold = 0.0;
	CookStats               stats;
	int                     passes = 0;   // how many passes the snapshot holds
	std::string             device;       // what marched: the GPU, its driver and queue priority
	// Expected resolved stars per texel and bin, the mean of the passes so far: what the stars were drawn from.
	// Empty unless CookOptions::keepBudget asked for it; at face 2048 it is 2.4 GB.
	std::vector<std::array<float, kBins>> budget;
};

struct CookOptions
{
	unsigned                                                 threads = 1;   // calibration and star draws
	std::optional<Quality>                                   quality{};    // overrides the scene's
	std::function<void(std::string_view phase, double done)> progress{};
	// Called after every pass with a complete sky; an error aborts the cook.
	std::function<Expected<void>(Cooked const& cooked, int pass, int passes)> onPass{};
	int                                                      stopAfterPasses = 0;   // test hook: 0 means all passes
	bool                                                     keepBudget = false;   // fill Cooked::budget (tests compare it)
};

// Same bound as the scene's target_stars. Past it a cook takes hours and writes a file nobody can load into a
// sky, and the per-texel-bin star index would outgrow the 24 bits the random stream gives it.
inline constexpr double kMaxStars = 1.0e7;
inline constexpr double kMaxStarsPerTexelBin = double(1u << 24) - double(1u << 20);

struct Stratum
{
	int x;
	int y;
};
// Pass `pass` of perAxis² visits this stratum; every prefix of passes stays spread over the texel.
// perAxis must be a power of two.
Stratum PassStratum(int perAxis, int pass);

// The running sums a pass adds to: per texel, Σ weight · radiance and Σ budget. Float, because the GPU holds
// and accumulates them in float and every pass rounds them through it; double here would only cost memory.
struct PassSums
{
	std::vector<glm::vec3>                radiance;
	std::vector<std::array<float, kBins>> budget;
};

// What a pass needs to fire one ray per texel.
struct PassSpec
{
	uint32_t faceSize;
	int      perAxis;
	Stratum  stratum;
	double   fluxThreshold;
	uint64_t seed;
};

// Fires one ray per texel through `spec.stratum` and adds it to `sums`. The GPU march implements it; the interface
// keeps Vulkan out of everything that includes the cook.
class PassMarcher
{
public:
	virtual ~PassMarcher() = default;
	virtual Expected<void> MarchPass(PassSpec const& spec, PassSums& sums) = 0;
	// What marches, for the cook's last log line: a file's bytes depend on it.
	virtual std::string Name() const = 0;
};

// A GPU dispatch's target duration: short enough to keep the desktop responsive and the driver's watchdog quiet.
inline constexpr double kTileSeconds = 0.25;
// The next dispatch's rows: scaled toward `target` seconds, at most doubling, clamped to [1, faceSize].
int NextTileRows(int rows, double seconds, int faceSize, double target = kTileSeconds);

// An input error when no GPU can march the scene (saying why), or when the threshold would resolve more than
// kMaxStars points.
Expected<Cooked> Cook(Scene const& scene, CookOptions const& options);

FluxHistogram CalibrationHistogram(PreparedScene const& scene, MarchParams const& params, uint32_t faceSize, unsigned threads);
uint32_t SamplePoisson(double lambda, RandomStream& rng);
glm::dvec3 SampleTexelDirection(uint32_t faceSize, uint32_t texel, RandomStream& rng);
glm::dvec3 TexelCentre(uint32_t faceSize, uint32_t texel);

}
