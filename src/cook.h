#pragma once
#include "common.h"
#include "march.h"
#include "population.h"
#include "rng.h"
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
	std::string             device;       // what marched, e.g. "CPU, 8 threads"
	// Expected resolved stars per texel and bin, the mean of the passes so far: what the stars were drawn from.
	// Empty unless CookOptions::keepBudget asked for it; at face 2048 it is 2.4 GB.
	std::vector<std::array<float, kBins>> budget;
};

enum class DeviceChoice { Auto, Gpu, Cpu };

struct CookOptions
{
	unsigned                                                 threads = 1;
	std::optional<Quality>                                   quality{};    // overrides the scene's
	std::function<void(std::string_view phase, double done)> progress{};
	// Called after every pass with a complete sky; an error aborts the cook.
	std::function<Expected<void>(Cooked const& cooked, int pass, int passes)> onPass{};
	int                                                      stopAfterPasses = 0;   // test hook: 0 means all passes
	DeviceChoice                                             device = DeviceChoice::Auto;
	std::function<void(std::string_view)>                    log{};   // one line per notable choice, e.g. falling back to the CPU
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

// The direction of the ray through stratum `s` of `texel`, jittered by (u, v) in [0, 1).
glm::dvec3 TexelJitterDirection(uint32_t faceSize, uint32_t texel, int perAxis, Stratum s, double u, double v);

// The running sums a pass adds to: per texel, Σ weight · radiance and Σ budget.
struct PassSums
{
	std::vector<glm::dvec3>               radiance;
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

// Fires one ray per texel through `spec.stratum` and adds it to `sums`. The CPU and GPU march implement it.
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

// An input error when the threshold would resolve more than kMaxStars points.
Expected<Cooked> Cook(Scene const& scene, CookOptions const& options);

FluxHistogram CalibrationHistogram(PreparedScene const& scene, MarchParams const& params, uint32_t faceSize, unsigned threads);
uint32_t SamplePoisson(double lambda, RandomStream& rng);
glm::dvec3 SampleTexelDirection(uint32_t faceSize, uint32_t texel, RandomStream& rng);
glm::dvec3 TexelCentre(uint32_t faceSize, uint32_t texel);

}
