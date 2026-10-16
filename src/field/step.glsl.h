// What the march does at one step: the per-population tables and how a step's stars are accumulated.
#ifndef STARSCAPE_FIELD_STEP_GLSL_H
#define STARSCAPE_FIELD_STEP_GLSL_H

#include "galaxy.glsl.h"

#ifdef __cplusplus
namespace starscape
{
#endif

FIELD_CONST int kBins = 24;

// What the march reads per step, derived once per population.
struct PopulationTable
{
	ARRAY(real3, chroma, kBins) ZEROED;       // linear Rec.709, luminance 1
	ARRAY(real, fraction, kBins) ZEROED;
	ARRAY(real, luminosity, kBins) ZEROED;
	// lightBelow[i] = Σ_{j<i} fraction_j × luminosity_j × chroma_j: the diffuse light of every
	// bin below a resolution split, so a step costs one lookup instead of a loop over bins.
	ARRAY(real3, lightBelow, kBins + 1) ZEROED;
};

struct RayResult
{
	real3 radiance DEFAULT(real3(0.0));   // mag-0 flux per steradian, linear Rec.709
	ARRAY(real, budget, kBins) ZEROED;    // expected resolved stars per luminosity bin
	// The expected light of those resolved stars, in radiance units. Not used by the cook; it exists so the
	// split can be checked exactly: radiance + resolvedRadiance is the same with any threshold.
	real3 resolvedRadiance DEFAULT(real3(0.0));
};

// Channel extinction relative to V approximates Cardelli R_V = 3.1 at 610, 550 and 465 nm;
// 1.0857 = 2.5 / ln 10 turns magnitudes into optical depth.
FIELD_INLINE real3 TransmittanceForDepth(real magnitudes)
{
	const real tau = magnitudes / 1.0857;
	return real3(exp(-0.75 * tau), exp(-1.00 * tau), exp(-1.32 * tau));
}

// Lowest luminosity resolved at distance x: L T_V 100 C / x^2 >= F_res.
FIELD_INLINE real MinResolvedLuminosity(real fluxThreshold, real tv, real x)
{
	if (!(tv > 0.0))
	{
		return kHugeReal;
	}
	return fluxThreshold * x * x / (100.0 * kFluxAt10pc * tv);
}

// Distance out to which a star of luminosity L is resolved.
FIELD_INLINE real ResolvedDistance(real fluxThreshold, real tv, real luminosity)
{
	return sqrt(luminosity * 100.0 * kFluxAt10pc * tv / fluxThreshold);
}

// (y³ − x³)/3 factored, so a float march far from the observer does not lose the span to cancellation.
FIELD_INLINE real CubeSpan(real x, real y)
{
	return (y - x) * (x * x + x * y + y * y) / 3.0;
}

// CubeSpan(x, x + w) from the width itself: far from the observer a float b − a keeps only a few bits of a short
// step, while the step's own width is exact.
FIELD_INLINE real CubeSpanFrom(real x, real w)
{
	const real y = x + w;
	return w * (x * x + x * y + y * y) / 3.0;
}

// The first bin whose luminosity is not below `minLuminosity`; luminosities ascend.
FIELD_INLINE int FirstResolvedBin(IN(PopulationTable) population, real minLuminosity)
{
	int bin = 0;
	while (bin < kBins && population.luminosity[bin] < minLuminosity)
	{
		bin += 1;
	}
	return bin;
}

// One population's stars over a step [a, b]: density and transmittance constant, distance integrated
// exactly, so the resolved-or-diffuse split is not quantised to the step grid (a split decided at each step's
// midpoint showed as ±9% ripples in star counts). `ds` is passed rather than
// computed as b − a, which can differ from the walk's step in the last bit and would move the CPU's bytes.
FIELD_INLINE void AccumulatePart(IN(PopulationTable) population, real density, real fluxThreshold, real a, real b,
	real ds, real3 transmittance, real weight, real sampleSolidAngle, INOUT(RayResult) result)
{
	const int splitA = FirstResolvedBin(population, MinResolvedLuminosity(fluxThreshold, transmittance.y, a));
	const int splitB = FirstResolvedBin(population, MinResolvedLuminosity(fluxThreshold, transmittance.y, b));
	// n s^2 dOmega ds stars at flux 100 L C / s^2 each: the s^2 cancels, leaving radiance 100 n L C ds.
	result.radiance += (weight * density * 100.0 * kFluxAt10pc * ds) * (transmittance * population.lightBelow[splitA]);
	const real lightScale = weight * density * 100.0 * kFluxAt10pc;
	for (int bin = splitA; bin < splitB; ++bin)
	{
		// Resolved out to c, diffuse beyond it.
		const real c = clamp(ResolvedDistance(fluxThreshold, transmittance.y, population.luminosity[bin]), a, b);
		const real f = population.fraction[bin];
		const real3 binLight = transmittance * population.chroma[bin] * (f * population.luminosity[bin]);
		result.budget[bin] += density * f * sampleSolidAngle * CubeSpan(a, c);
		result.radiance += (lightScale * (b - c)) * binLight;
		result.resolvedRadiance += (lightScale * (c - a)) * binLight;
	}
	for (int bin = splitB; bin < kBins; ++bin)
	{
		const real f = population.fraction[bin];
		result.budget[bin] += density * f * sampleSolidAngle * CubeSpanFrom(a, ds);
		result.resolvedRadiance += (lightScale * ds) * (transmittance * population.chroma[bin] * (f * population.luminosity[bin]));
	}
}

#ifdef __cplusplus
}
#endif

#endif
