#pragma once
#include <array>
#include <cstdint>
#include <string>
#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>

namespace starscape
{

// Logarithmic spiral arms: arm k sits at azimuth phase + 2πk/count + ln(R / r0) / tan(pitch).
struct Arms
{
	int    count     = 0;
	double pitchDeg  = 12.0;
	double width     = 400.0;   // Gaussian σ across the arm, parsecs of arc length
	double contrast  = 1.0;     // density multiplier on the arm's crest
	double clumping  = 0.0;     // log standard deviation of the mean-one field that breaks arms into clumps and gaps
	double phaseDeg  = 0.0;
	double r0        = 1000.0;
};

// The most colours a gas palette takes: enough for any ramp, small enough to live in a fixed array.
inline constexpr int kMaxPalette = 8;

// Gas lit by a disk's young stars. Its light is `glow` times theirs, averaged over space: it gathers where
// they cluster and onto the dust's filaments, and the dust's bubble clears it too.
struct Gas
{
	double     glow = 0.0;
	glm::dvec3 colour{1.0, 0.18, 0.32};      // hydrogen red; PrepareGalaxy normalises it to luminance 1
	glm::dvec3 coreColour{0.25, 1.0, 0.85};  // the densest clumps; normalised likewise
	// A cyclic ramp replacing `colour`, chosen by a field of `paletteWavelength`, so that neighbouring
	// clouds differ in hue. Empty (count 0): `colour` alone. PrepareGalaxy normalises each to luminance 1.
	std::array<glm::dvec3, kMaxPalette> palette{};
	int        paletteCount = 0;
	double     paletteWavelength = 1500.0;   // parsecs
	double     textureMean = 1.0;            // derived by PrepareGalaxy: ⟨G⟩, divided out like Dust::filamentMean
	double     youngLight  = 0.0;            // derived by PrepareScene: the arm population's mean L_V per star
};

// Field order matters: presets and tests use designated initialisers.
struct Disk
{
	bool   enabled          = false;
	double density          = 0.0;      // stars/pc³ at the centre of the midplane
	double scaleLength      = 1.0;
	double scaleHeight      = 1.0;
	double cutoff           = 0.0;      // spherical radius beyond which the component is zero
	Arms   arms;
	int    population       = -1;
	int    armPopulation    = -1;       // -1: arms modulate density only
	double armYoungFraction = 0.0;      // share of armPopulation on an arm's crest, at the midplane
	double armScaleHeight   = 0.0;      // the arm population's own scale height; 0: the disk's
	double youngFraction    = 0.0;      // armPopulation share everywhere in the midplane, arms or not
	double associationSize  = 0.0;      // wavelength of the field gathering young stars into OB associations; 0: none
	double associationClumping = 0.0;   // its log standard deviation
	double associationMean  = 1.0;      // derived by PrepareGalaxy, divided out like Dust::filamentMean
	Gas    gas{};                       // thin disks only; needs armPopulation. The {} keeps -Wextra quiet where designated initialisers omit it
};

// Hernquist profile n_b a / (r' (r' + a)³) over a radius r' stretched along a bar.
struct Bulge
{
	bool   enabled     = false;
	double density     = 0.0;   // n_b; total stars = 2π n_b before the cutoff
	double scale       = 1.0;   // a
	double axisY       = 1.0;   // bar axis ratios relative to its long axis
	double axisZ       = 1.0;
	double barAngleDeg = 0.0;
	double cutoff      = 0.0;
	int    population  = -1;
};

struct Halo
{
	bool   enabled    = false;
	double density    = 0.0;   // n_h in n_h (1 + r / scale)^-3.5
	double scale      = 1.0;
	double cutoff     = 0.0;
	int    population = -1;
};

struct Dust
{
	bool   enabled     = false;
	double extinction  = 0.0;   // α₀, visual magnitudes per parsec at the centre of the midplane
	double scaleLength = 1.0;
	double scaleHeight = 1.0;
	double cutoff      = 0.0;
	double laneOffset  = 0.0;   // lane shift toward each arm's inner edge, in arm widths
	double contrast    = 1.0;
	double filament    = 0.0;   // log standard deviation of the field that breaks dust into filaments and holes
	// A supernova-blown cavity such as the Sun's Local Bubble: an ellipsoid in the galaxy frame, empty
	// within 70% of its semi-axes and ramping to full density at its edge. Radius 0: none.
	glm::dvec3 bubbleCenter{0.0};
	double bubbleRadius = 0.0;  // semi-axis in the plane, parsecs
	double bubbleHeight = 0.0;  // semi-axis along the spin axis; bubbles open into chimneys toward the halo
	double filamentMean = 1.0;  // derived by PrepareGalaxy: the field's spatial mean, divided out so α₀ stays the mean
};

struct Galaxy
{
	std::string preset;
	glm::dvec3  position{0.0};      // scene frame, parsecs
	glm::dvec3  rotationDeg{0.0};   // Euler angles applied Z, then Y, then X
	uint64_t    seed = 1;
	Disk        disk;
	Disk        thickDisk;
	Bulge       bulge;
	Halo        halo;
	Dust        dust;
};

struct Contribution
{
	double density    = 0.0;   // stars/pc³
	int    population = -1;
};

// At most six parts: each disk splits into base + arm population, plus bulge and halo.
struct GalaxySample
{
	std::array<Contribution, 6> parts{};
	int                         count      = 0;
	double                      extinction = 0.0;   // α_V, magnitudes per parsec
	glm::dvec3                  gas{0.0};           // glowing gas: radiance per parsec, before transmittance
};

struct PreparedGalaxy
{
	Galaxy     galaxy;
	glm::dmat3 toLocal{1.0};
	glm::dvec3 position{0.0};
	double     boundRadius = 0.0;
};

double Sech2(double x);
double ArmProximity(Arms const& arms, double radius, double azimuth, double offset);
double BoundRadius(Galaxy const& galaxy);
PreparedGalaxy PrepareGalaxy(Galaxy const& galaxy);
void EvaluateLocal(Galaxy const& galaxy, glm::dvec3 local, GalaxySample& out);

inline void Evaluate(PreparedGalaxy const& prepared, glm::dvec3 scenePosition, GalaxySample& out)
{
	EvaluateLocal(prepared.galaxy, prepared.toLocal * (scenePosition - prepared.position), out);
}

}
