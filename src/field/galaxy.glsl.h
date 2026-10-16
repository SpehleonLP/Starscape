// A galaxy's description and its per-point evaluation, shared by the host and the GPU.
#ifndef STARSCAPE_FIELD_GALAXY_GLSL_H
#define STARSCAPE_FIELD_GALAXY_GLSL_H

#include "noise.glsl.h"

#ifdef __cplusplus
#include <string>

namespace starscape
{
#endif

// Logarithmic spiral arms: arm k sits at azimuth phase + 2πk/count + ln(R / r0) / tan(pitch).
struct Arms
{
	int  count     DEFAULT(0);
	real pitchDeg  DEFAULT(12.0);
	real width     DEFAULT(400.0);   // Gaussian σ across the arm, parsecs of arc length
	real contrast  DEFAULT(1.0);     // density multiplier on the arm's crest
	real clumping  DEFAULT(0.0);     // log standard deviation of the mean-one field that breaks arms into clumps and gaps
	real phaseDeg  DEFAULT(0.0);
	real r0        DEFAULT(1000.0);
};

// The most colours a gas palette takes: enough for any ramp, small enough to live in a fixed array.
FIELD_CONST int kMaxPalette = 8;

// Gas lit by a disk's young stars. Its light is `glow` times theirs, averaged over space: it gathers where
// they cluster and onto the dust's filaments, and the dust's bubble clears it too.
struct Gas
{
	real  glow DEFAULT(0.0);
	real3 colour DEFAULT(real3(1.0, 0.18, 0.32));       // hydrogen red; PrepareGalaxy normalises it to luminance 1
	real3 coreColour DEFAULT(real3(0.25, 1.0, 0.85));   // the densest clumps; normalised likewise
	// A cyclic ramp replacing `colour`, chosen by a field of `paletteWavelength`, so that neighbouring
	// clouds differ in hue. Empty (count 0): `colour` alone. PrepareGalaxy normalises each to luminance 1.
	ARRAY(real3, palette, kMaxPalette) ZEROED;
	int   paletteCount DEFAULT(0);
	real  paletteWavelength DEFAULT(1500.0);   // parsecs
	// A second gas band at the size of giant cloud complexes. `largeGlow` is its light as a multiple of the
	// young stars' own, like `glow`.
	real  largeGlow DEFAULT(0.0);
	real  largeWavelength DEFAULT(100.0);  // parsecs
	real  largeMean DEFAULT(1.0);          // derived by PrepareGalaxy: the field's mean in noise space, divided out
	real  textureMean DEFAULT(1.0);        // derived by PrepareGalaxy: ⟨G⟩, divided out like Dust::filamentMean
	real  youngLight DEFAULT(0.0);         // derived by PrepareScene: the arm population's mean L_V per star
};

// Field order matters: presets and tests use designated initialisers.
struct Disk
{
	bool  enabled DEFAULT(false);
	real  density DEFAULT(0.0);      // stars/pc³ at the centre of the midplane
	real  scaleLength DEFAULT(1.0);
	real  scaleHeight DEFAULT(1.0);
	real  cutoff DEFAULT(0.0);       // spherical radius beyond which the component is zero
	Arms  arms;
	int   population DEFAULT(-1);
	int   armPopulation DEFAULT(-1);        // -1: arms modulate density only
	real  armYoungFraction DEFAULT(0.0);    // share of armPopulation on an arm's crest, at the midplane
	real  armScaleHeight DEFAULT(0.0);      // the arm population's own scale height; 0: the disk's
	real  youngFraction DEFAULT(0.0);       // armPopulation share everywhere in the midplane, arms or not
	real  associationSize DEFAULT(0.0);     // wavelength of the field gathering young stars into OB associations; 0: none
	real  associationClumping DEFAULT(0.0); // its log standard deviation
	real  associationMean DEFAULT(1.0);     // derived by PrepareGalaxy, divided out like Dust::filamentMean
	Gas   gas DEFAULT(Gas());               // thin disks only; needs armPopulation. The explicit Gas() keeps -Wextra quiet where designated initialisers omit it
};

// Hernquist profile n_b a / (r' (r' + a)³) over a radius r' stretched along a bar.
struct Bulge
{
	bool enabled DEFAULT(false);
	real density DEFAULT(0.0);   // n_b; total stars = 2π n_b before the cutoff
	real scale DEFAULT(1.0);     // a
	real axisY DEFAULT(1.0);     // bar axis ratios relative to its long axis
	real axisZ DEFAULT(1.0);
	real barAngleDeg DEFAULT(0.0);
	real cutoff DEFAULT(0.0);
	int  population DEFAULT(-1);
};

struct Halo
{
	bool enabled DEFAULT(false);
	real density DEFAULT(0.0);   // n_h in n_h (1 + r / scale)^-3.5
	real scale DEFAULT(1.0);
	real cutoff DEFAULT(0.0);
	int  population DEFAULT(-1);
};

struct Dust
{
	bool enabled DEFAULT(false);
	real extinction DEFAULT(0.0);   // α₀, visual magnitudes per parsec at the centre of the midplane
	real scaleLength DEFAULT(1.0);
	real scaleHeight DEFAULT(1.0);
	real cutoff DEFAULT(0.0);
	real laneOffset DEFAULT(0.0);   // lane shift toward each arm's inner edge, in arm widths
	real contrast DEFAULT(1.0);
	real filament DEFAULT(0.0);     // log standard deviation of the field that breaks dust into filaments and holes
	// A supernova-blown cavity such as the Sun's Local Bubble: an ellipsoid in the galaxy frame, empty
	// within 70% of its semi-axes and ramping to full density at its edge. Radius 0: none.
	real3 bubbleCenter DEFAULT(real3(0.0));
	real  bubbleRadius DEFAULT(0.0);  // semi-axis in the plane, parsecs
	real  bubbleHeight DEFAULT(0.0);  // semi-axis along the spin axis; bubbles open into chimneys toward the halo
	real  filamentMean DEFAULT(1.0);  // derived by PrepareGalaxy: the field's spatial mean, divided out so α₀ stays the mean
	// Gas glowing on the bubble's wall, lit in patches: arcs across the sky from inside it. Its strength
	// is a multiple of the young stars' light at the wall's radius; it needs the thin disk's arm population.
	real  bubbleRim DEFAULT(0.0);
	// Dust lit by nearby young stars scatters their light. `reflection` is that light as a multiple of the
	// young stars' own where the dust runs at 1 mag/kpc. Added light, a stylisation: the stars lose nothing.
	real  reflection DEFAULT(0.0);
	real3 reflectionColour DEFAULT(real3(0.45, 0.62, 1.0));   // PrepareGalaxy normalises it to luminance 1
};

struct Galaxy
{
#ifdef __cplusplus
	std::string preset;
#endif
	real3    position DEFAULT(real3(0.0));      // scene frame, parsecs
	real3    rotationDeg DEFAULT(real3(0.0));   // Euler angles applied Z, then Y, then X
	uint64_t seed DEFAULT(1);
	Disk     disk;
	Disk     thickDisk;
	Bulge    bulge;
	Halo     halo;
	Dust     dust;
};

struct Contribution
{
	real density DEFAULT(0.0);   // stars/pc³
	int  population DEFAULT(-1);
};

// At most six parts: each disk splits into base + arm population, plus bulge and halo.
struct GalaxySample
{
	ARRAY(Contribution, parts, 6) ZEROED;
	int   count DEFAULT(0);
	real  extinction DEFAULT(0.0);   // α_V, magnitudes per parsec
	real3 gas DEFAULT(real3(0.0));   // glowing gas: radiance per parsec, before transmittance
};

// Flux of a 1 L_V☉ source at 10 pc in magnitude-0 units: 10^(−0.4 × 4.83), the Sun's M_V. Here rather than
// with the march step because the gas's light is scaled by it too.
FIELD_CONST real kFluxAt10pc = 0.01169499;

// Salts so the dust filaments do not trace the same noise as the arm clumps.
FIELD_CONST uint64_t kDustNoiseSalt = 0xD057D057D057D057ul;

// Moments measured over 4·10⁵ random points. Ridged fbm at 6 octaves (Fbm's σ is kFbmSigma).
FIELD_CONST real kRidgedMean = 0.6307;
FIELD_CONST real kRidgedSigma = 0.1339;
FIELD_CONST int kDustOctaves = 6;
FIELD_CONST int kAssociationOctaves = 3;
FIELD_CONST uint64_t kAssociationNoiseSalt = 0xA550C1A7E5ul;

// The gas texture's noise. Fixed rather than scene keys: the owner-approved look was tuned on exactly these.
FIELD_CONST uint64_t kGasClumpSalt = 0x6A5ul;
FIELD_CONST uint64_t kGasWispSalt = 0x6A6ul;
FIELD_CONST real kGasClumpWavelength = 40.0;   // parsecs: a bright nebula's core
FIELD_CONST real kGasWispWavelength = 25.0;
FIELD_CONST int kGasOctaves = 5;
// The palette's hue field. Fixed rather than a scene key: the owner-approved nebula look was tuned on it.
FIELD_CONST uint64_t kPaletteSalt = 0x9A1E77Eul;
FIELD_CONST uint64_t kLargeGasSalt = 0xB16B16ul;
FIELD_CONST int kLargeGasOctaves = 4;
// The rim's patches and wisps. Fixed rather than scene keys: the owner-approved nebula look was tuned on them.
FIELD_CONST uint64_t kRimMaskSalt = 0x81B1ul;
FIELD_CONST uint64_t kRimWispSalt = 0x81B0ul;

FIELD_INLINE real Radians(real d)
{
	return d * 3.14159265358979323846 / 180.0;
}

// Dust is filamentary, so its field rides ridged noise; that noise is bounded above, which pulls the
// log-normal's mean below one (0.65 at σ 2). Measured here rather than tabulated, so any σ stays calibrated.
FIELD_INLINE real DustField(IN(Dust) dust, real3 p, real wavelength, uint64_t seed)
{
	return LogNormal(dust.filament, (RidgedFbmAt(p, wavelength, kDustOctaves, seed ^ kDustNoiseSalt) - kRidgedMean) / kRidgedSigma);
}

// Young stars are born in OB associations tens to a couple of hundred parsecs across, so the bright
// stars, which are mostly young, gather into groups like Orion's while the old field stays smooth.
FIELD_INLINE real AssociationField(IN(Disk) disk, real3 p, real wavelength, uint64_t seed)
{
	return LogNormal(disk.associationClumping, FbmAt(p, wavelength, kAssociationOctaves, seed ^ kAssociationNoiseSalt) / kFbmSigma);
}

FIELD_INLINE real DustWavelength(IN(Galaxy) galaxy)
{
	return 0.5 * (galaxy.disk.arms.count > 0 ? galaxy.disk.arms.width : galaxy.dust.scaleHeight);
}

// Log-normal like the clumps: dust is filaments with clear holes between them, not a modulated sheet.
// Mean one over space: the dust's filaments with the extinction's normalisation.
FIELD_INLINE real DustFilaments(IN(Galaxy) galaxy, real3 p)
{
	return DustField(galaxy.dust, p, DustWavelength(galaxy), galaxy.seed) / galaxy.dust.filamentMean;
}

// The bubble's normalised radius: 1 on its ellipsoid. Only meaningful when it has a radius.
FIELD_INLINE real BubbleRadius(IN(Dust) dust, real3 p)
{
	const real3 d = p - dust.bubbleCenter;
	return sqrt((d.x * d.x + d.y * d.y) / (dust.bubbleRadius * dust.bubbleRadius) + d.z * d.z / (dust.bubbleHeight * dust.bubbleHeight));
}

// Empty within 70% of the bubble's semi-axes, ramping smoothly to one at its edge; one without a bubble.
// Without the cavity, filaments within a few hundred parsecs of an observer in the plane blot the high
// latitudes, where the real sky's are faint cirrus.
FIELD_INLINE real BubbleFactor(IN(Dust) dust, real3 p)
{
	if (dust.bubbleRadius <= 0.0)
	{
		return 1.0;
	}
	const real t = clamp((BubbleRadius(dust, p) - 0.7) / 0.3, 0.0, 1.0);
	return t * t * (3.0 - 2.0 * t);
}

// Mean one: the soft clumps that also pick the gas colour.
FIELD_INLINE real GasClumps(IN(Galaxy) galaxy, real3 p)
{
	return LogNormal(1.2, FbmAt(p, kGasClumpWavelength, kGasOctaves, galaxy.seed ^ kGasClumpSalt) / kFbmSigma);
}

// G before its mean is divided out: clumps plus thin bright wisps, gathered onto the dust's filaments
// because gas and dust are one medium.
FIELD_INLINE real GasTexture(IN(Galaxy) galaxy, real3 p, real clumps)
{
	const real wisps = 6.0 * pow(max(0.0, RidgedFbmAt(p, kGasWispWavelength, kGasOctaves, galaxy.seed ^ kGasWispSalt)), 6.0);
	const real filaments = galaxy.dust.enabled ? DustFilaments(galaxy, p) : 1.0;
	return (0.6 * clumps + wisps) * sqrt(filaments);
}

// A one-σ change in the hue field moves about one palette entry, so a sky shows every entry without
// flickering between them. Entries have luminance 1, so any blend of two does too.
FIELD_INLINE real3 PaletteColour(IN(Galaxy) galaxy, IN(Gas) gas, real3 p)
{
	if (gas.paletteCount == 0)
	{
		return gas.colour;
	}
	const real count = real(gas.paletteCount);
	const real z = FbmAt(p, gas.paletteWavelength, 3, galaxy.seed ^ kPaletteSalt) / kFbmSigma;
	real u = FMod(0.9 * z, count);
	if (u < 0.0)
	{
		u += count;
	}
	const real index = floor(u);
	real f = u - index;
	f = f * f * (3.0 - 2.0 * f);
	const int i = int(index) % gas.paletteCount;
	return mix(gas.palette[i], gas.palette[(i + 1) % gas.paletteCount], f);
}

FIELD_INLINE real LargeGasField(real3 p, real wavelength, uint64_t seed)
{
	return LogNormal(1.0, FbmAt(p, wavelength, kLargeGasOctaves, seed ^ kLargeGasSalt) / kFbmSigma);
}

FIELD_INLINE real ArmProximity(IN(Arms) arms, real radius, real azimuth, real offset)
{
	// ln(R / r0) diverges on the axis, where arms have no meaning anyway.
	if (arms.count <= 0 || radius < 1e-6)
	{
		return 0.0;
	}
	const real period = 2.0 * 3.14159265358979323846 / arms.count;
	const real phase = Radians(arms.phaseDeg) + log(radius / arms.r0) / tan(Radians(arms.pitchDeg));
	real delta = azimuth - phase;
	delta -= period * floor(delta / period + 0.5);
	// Arc length to the nearest arm. Positive delta is the inner edge: the arm through that
	// azimuth crosses it at a larger radius.
	const real arc = delta * radius - offset;
	return exp(-arc * arc / (2.0 * arms.width * arms.width));
}

// A shell just inside the dust's ramp, so most of it sits in front of the wall's dust. Not confined to the
// young disk's thickness: the bubble opens into chimneys, so the wall wraps the whole sky.
FIELD_INLINE real3 RimGlow(IN(Galaxy) galaxy, real radius, real3 p)
{
	if (galaxy.dust.bubbleRim <= 0.0 || !galaxy.dust.enabled || galaxy.dust.bubbleRadius <= 0.0 || !galaxy.disk.enabled)
	{
		return real3(0.0);
	}
	const real offset = BubbleRadius(galaxy.dust, p) - 0.85;
	// Cut where the Gaussian is below e^-25, so the term costs nothing away from the wall.
	if (abs(offset) > 0.2)
	{
		return real3(0.0);
	}
	// Only patches of the wall glow: the eye finds arcs where a lit patch's edge runs along the line of sight.
	const real mask = clamp(FbmAt(p, 150.0, 3, galaxy.seed ^ kRimMaskSalt) / kFbmSigma - 0.5, 0.0, 1.0);
	if (mask <= 0.0)
	{
		return real3(0.0);
	}
	const real shell = exp(-(offset / 0.04) * (offset / 0.04));
	const real wisps = 0.3 + 6.0 * pow(max(0.0, RidgedFbmAt(p, 30.0, 5, galaxy.seed ^ kRimWispSalt)), 6.0);
	// The young density the disk would have at this radius, without arms, associations or its thickness.
	const real reference = galaxy.disk.density * exp(-radius / galaxy.disk.scaleLength) * (galaxy.disk.youngFraction + galaxy.disk.armYoungFraction);
	const real light = galaxy.dust.bubbleRim * galaxy.disk.gas.youngLight * reference * shell * mask * wisps;
	// Four times the palette's frequency gives the rim its own, finer colour changes.
	return (light * 100.0 * kFluxAt10pc) * PaletteColour(galaxy, galaxy.disk.gas, 4.0 * p);
}

FIELD_INLINE void Push(INOUT(GalaxySample) result, real density, int population)
{
	if (density > 0.0 && population >= 0)
	{
		result.parts[result.count].density = density;
		result.parts[result.count].population = population;
		result.count += 1;
	}
}

// Adds one disk's stars and gas; returns its young density, which the thin disk's reflection uses.
FIELD_INLINE real AddDisk(IN(Galaxy) galaxy, IN(Disk) disk, real3 p, real radius, real r, real azimuth, INOUT(GalaxySample) result)
{
	if (!disk.enabled || r > disk.cutoff)
	{
		return 0.0;
	}
	const real crest = ArmProximity(disk.arms, radius, azimuth, 0.0);
	real modulation = 1.0;
	real clumps = 1.0;
	if (disk.arms.count > 0)
	{
		// Log-normal, so clumping opens real gaps between star clouds instead of a ripple on the crest.
		clumps = LogNormal(disk.arms.clumping, FbmAt(p, disk.arms.width, 4, galaxy.seed) / kFbmSigma);
		modulation = max(0.0, 1.0 + (disk.arms.contrast - 1.0) * crest * clumps);
	}
	const real midplane = disk.density * exp(-radius / disk.scaleLength) * modulation;
	// Young stars concentrate on arms but are not confined to them: the Sun sits between arms,
	// yet half the brightest stars in Earth's sky are young B stars.
	const real young = disk.armPopulation >= 0 ? min(1.0, (disk.youngFraction + disk.armYoungFraction * crest) * clumps) : 0.0;
	Push(result, midplane * (1.0 - young) * Sech2(p.z / (2.0 * disk.scaleHeight)), disk.population);
	// Young stars have not yet been heated out of the plane: sharing the old disk's thickness turns
	// each nearby arm into a wall whose tangent draws a hard edge from horizon to zenith.
	const real youngHeight = disk.armScaleHeight > 0.0 ? disk.armScaleHeight : disk.scaleHeight;
	// Associations add young stars where they are rather than taking a share of the old ones, so the
	// field multiplies the young density, not the share capped at one.
	real associations = 1.0;
	if (disk.associationSize > 0.0)
	{
		associations = AssociationField(disk, p, disk.associationSize, galaxy.seed) / disk.associationMean;
	}
	const real youngDensity = midplane * young * associations * Sech2(p.z / (2.0 * youngHeight));
	Push(result, youngDensity, disk.armPopulation);
	if (youngDensity > 0.0 && (disk.gas.glow > 0.0 || disk.gas.largeGlow > 0.0))
	{
		const real3 base = PaletteColour(galaxy, disk.gas, p);
		const real bubble = BubbleFactor(galaxy.dust, p);
		if (disk.gas.glow > 0.0)
		{
			const real gasClumps = GasClumps(galaxy, p);
			const real gasTexture = GasTexture(galaxy, p, gasClumps) / disk.gas.textureMean;
			// Dense clumps turn from the base colour to the core colour; both have luminance 1, so the mix keeps the light.
			const real t = clamp((gasClumps - 1.0) / 2.0, 0.0, 1.0);
			const real3 colour = mix(base, disk.gas.coreColour, t * t * (3.0 - 2.0 * t));
			const real light = disk.gas.glow * disk.gas.youngLight * youngDensity * gasTexture * bubble;
			result.gas += (light * 100.0 * kFluxAt10pc) * colour;
		}
		if (disk.gas.largeGlow > 0.0)
		{
			const real field = LargeGasField(p, disk.gas.largeWavelength, galaxy.seed) / disk.gas.largeMean;
			const real light = disk.gas.largeGlow * disk.gas.youngLight * youngDensity * field * bubble;
			result.gas += (light * 100.0 * kFluxAt10pc) * base;
		}
	}
	return youngDensity;
}

FIELD_INLINE void EvaluateLocal(IN(Galaxy) galaxy, real3 p, INOUT(GalaxySample) result)
{
	result.count = 0;
	result.extinction = 0.0;
	result.gas = real3(0.0);
	const real radius = Hypot(p.x, p.y);
	const real r = length(p);
	const real azimuth = atan(p.y, p.x);

	// The thin disk's young stars light the reflection below.
	const real youngDensity = AddDisk(galaxy, galaxy.disk, p, radius, r, azimuth, result);
	AddDisk(galaxy, galaxy.thickDisk, p, radius, r, azimuth, result);
	if (r <= galaxy.disk.cutoff)
	{
		result.gas += RimGlow(galaxy, radius, p);
	}

	if (galaxy.bulge.enabled && r <= galaxy.bulge.cutoff)
	{
		const real c = cos(Radians(galaxy.bulge.barAngleDeg));
		const real s = sin(Radians(galaxy.bulge.barAngleDeg));
		const real bx = c * p.x + s * p.y;
		const real by = (-s * p.x + c * p.y) / galaxy.bulge.axisY;
		const real bz = p.z / galaxy.bulge.axisZ;
		// The Hernquist cusp diverges at r' = 0; a 1 pc floor keeps an observer at the core finite.
		const real rb = max(1.0, sqrt(bx * bx + by * by + bz * bz));
		const real a = galaxy.bulge.scale;
		Push(result, galaxy.bulge.density * a / (rb * (rb + a) * (rb + a) * (rb + a)), galaxy.bulge.population);
	}

	if (galaxy.halo.enabled && r <= galaxy.halo.cutoff)
	{
		Push(result, galaxy.halo.density * pow(1.0 + r / galaxy.halo.scale, -3.5), galaxy.halo.population);
	}

	if (galaxy.dust.enabled && r <= galaxy.dust.cutoff)
	{
		const real lane = ArmProximity(galaxy.disk.arms, radius, azimuth, galaxy.dust.laneOffset * galaxy.disk.arms.width);
		const real laneModulation = 1.0 + (galaxy.dust.contrast - 1.0) * lane;
		const real filaments = DustFilaments(galaxy, p);
		result.extinction = galaxy.dust.extinction * exp(-radius / galaxy.dust.scaleLength) * Sech2(p.z / (2.0 * galaxy.dust.scaleHeight)) * laneModulation * filaments
			* BubbleFactor(galaxy.dust, p);
		if (galaxy.dust.reflection > 0.0 && youngDensity > 0.0)
		{
			// α_V in magnitudes per kiloparsec, so `reflection` is a multiple of the young light at 1 mag/kpc.
			const real light = galaxy.dust.reflection * (result.extinction * 1000.0) * galaxy.disk.gas.youngLight * youngDensity;
			result.gas += (light * 100.0 * kFluxAt10pc) * galaxy.dust.reflectionColour;
		}
	}
}

#ifdef __cplusplus
}
#endif

#endif
