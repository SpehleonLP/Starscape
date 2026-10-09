#include "presets.h"
#include <utility>

namespace starscape
{

// Derivation (offline script, reproduced here so the numbers can be regenerated):
// - Kroupa (2001) IMF, ξ(m) ∝ m^-1.3 for 0.08–0.5 M☉ and m^-2.3 above, sampled 0.08–150 M☉.
// - Main sequence: L_bol = 0.23 m^2.3 (m < 0.43), m^4 (< 2), 1.4 m^3.5 (< 55), 32000 m;
//   R = m^0.8 (m < 1) or m^0.57; T = 5772 (L / R²)^0.25.
// - L_V = L_bol × (blackbody flux under the CIE ȳ curve / σT⁴), normalised so the Sun has L_V = 1.
// - A turnoff mass removes stars above it; stars within 4 % above the turnoff become giants
//   (the post-main-sequence lifetime is ~10 % of the main-sequence one, lifetime ∝ m^-2.5),
//   spread log-uniformly over L_V 10–1000 at T 4900 → 4000 K. Heavier stars are remnants and,
//   like every star below the lowest edge, are counted in bin 0.
// - young-disk: no turnoff. old-disk: 1.0 M☉. bulge and halo: 0.85 M☉.
// - T per bin is the light-weighted mean; empty bins repeat the previous bin's T.
// Mean L_V per star: young-disk 29.42, old-disk 1.124, bulge 1.295, halo 1.295.
std::vector<Population> PresetPopulations()
{
	Population young{"young-disk",
		{2.041409e-01, 1.388415e-01, 1.299330e-01, 1.212253e-01, 1.121052e-01, 5.238831e-02, 4.934281e-02, 4.202055e-02,
		 3.494578e-02, 2.863004e-02, 2.334727e-02, 1.869960e-02, 1.551932e-02, 1.111291e-02, 7.276815e-03, 4.573593e-03,
		 2.717156e-03, 1.546605e-03, 8.465797e-04, 4.846047e-04, 3.021698e-04, 0.0, 0.0, 0.0},
		{2680, 2842, 2999, 3170, 3354, 3639, 4026, 4495, 5060, 5764, 6780, 8146,
		 9957, 12263, 15408, 19718, 25664, 33876, 45214, 61254, 65806, 65806, 65806, 65806}};

	Population old{"old-disk",
		{2.963755e-01, 1.388415e-01, 1.299330e-01, 1.212253e-01, 1.121052e-01, 5.238831e-02, 4.934281e-02, 4.202055e-02,
		 3.494578e-02, 1.789629e-02, 0.0, 0.0, 1.022090e-03, 1.034405e-03, 1.022090e-03, 1.022090e-03,
		 8.250609e-04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{2680, 2842, 2999, 3170, 3354, 3639, 4026, 4495, 5060, 5567, 5567, 5567,
		 4792, 4604, 4416, 4229, 4066, 4066, 4066, 4066, 4066, 4066, 4066, 4066}};

	Population bulge{"bulge",
		{3.183034e-01, 1.388415e-01, 1.299330e-01, 1.212253e-01, 1.121052e-01, 5.238831e-02, 4.934281e-02, 4.202055e-02,
		 3.005734e-02, 0.0, 0.0, 0.0, 1.199895e-03, 1.214352e-03, 1.199895e-03, 1.199895e-03,
		 9.685900e-04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{2680, 2842, 2999, 3170, 3354, 3639, 4026, 4495, 5001, 5001, 5001, 5001,
		 4792, 4604, 4416, 4229, 4066, 4066, 4066, 4066, 4066, 4066, 4066, 4066}};

	Population halo = bulge;
	halo.name = "halo";

	return {std::move(young), std::move(old), std::move(bulge), std::move(halo)};
}

// Milky Way numbers: thin disk 0.1 stars/pc³ at the Sun's 8.2 kpc (n₀ = 0.1 e^{8200/2600});
// thick disk ~6 % of that locally; bulge n_b = 3.2 × 10⁹ (a spherical bulge would hold 2π n_b = 2 × 10¹⁰ stars; the 0.5 × 0.4 bar axes
// scale that by 0.2 to ≈ 4 × 10⁹, ≈ 3.3 × 10⁹ inside the 6 kpc cutoff) with a 25° bar;
// halo ≈ 10⁻⁴ stars/pc³ locally; dust ≈ 1 mag/kpc locally (α₀ = 0.001 e^{8200/3200}).
std::optional<Galaxy> PresetGalaxy(std::string_view name, std::span<Population const> populations)
{
	auto population = [&](std::string_view populationName) { return FindPopulation(populations, populationName).value_or(-1); };
	Galaxy galaxy;
	galaxy.preset = std::string(name);
	if (name == "milky-way")
	{
		galaxy.disk = {.enabled = true, .density = 2.34, .scaleLength = 2600.0, .scaleHeight = 300.0, .cutoff = 18000.0,
			.arms = {.count = 4, .pitchDeg = 12.0, .width = 400.0, .contrast = 1.5, .clumping = 1.2, .phaseDeg = 0.0, .r0 = 1000.0},
			.population = population("old-disk"), .armPopulation = population("young-disk"), .armYoungFraction = 0.1, .armScaleHeight = 90.0};
		galaxy.thickDisk = {.enabled = true, .density = 0.058, .scaleLength = 3600.0, .scaleHeight = 900.0, .cutoff = 22000.0,
			.arms = {}, .population = population("old-disk")};
		galaxy.bulge = {.enabled = true, .density = 3.2e9, .scale = 600.0, .axisY = 0.5, .axisZ = 0.4, .barAngleDeg = 25.0,
			.cutoff = 6000.0, .population = population("bulge")};
		galaxy.halo = {.enabled = true, .density = 0.24, .scale = 1000.0, .cutoff = 30000.0, .population = population("halo")};
		galaxy.dust = {.enabled = true, .extinction = 0.013, .scaleLength = 3200.0, .scaleHeight = 100.0, .cutoff = 18000.0,
			.laneOffset = 0.15, .contrast = 4.0, .filament = 1.2};
		return galaxy;
	}
	if (name == "small-spiral")
	{
		galaxy.disk = {.enabled = true, .density = 1.5, .scaleLength = 1500.0, .scaleHeight = 250.0, .cutoff = 10000.0,
			.arms = {.count = 2, .pitchDeg = 18.0, .width = 300.0, .contrast = 1.5, .clumping = 1.2, .phaseDeg = 0.0, .r0 = 600.0},
			.population = population("old-disk"), .armPopulation = population("young-disk"), .armYoungFraction = 0.12, .armScaleHeight = 80.0};
		galaxy.bulge = {.enabled = true, .density = 3.2e8, .scale = 300.0, .axisY = 1.0, .axisZ = 1.0, .barAngleDeg = 0.0,
			.cutoff = 3000.0, .population = population("bulge")};
		galaxy.halo = {.enabled = true, .density = 0.05, .scale = 800.0, .cutoff = 15000.0, .population = population("halo")};
		galaxy.dust = {.enabled = true, .extinction = 0.008, .scaleLength = 1800.0, .scaleHeight = 90.0, .cutoff = 10000.0,
			.laneOffset = 0.15, .contrast = 4.0, .filament = 1.2};
		return galaxy;
	}
	if (name == "elliptical")
	{
		galaxy.bulge = {.enabled = true, .density = 1.6e10, .scale = 3000.0, .axisY = 0.8, .axisZ = 0.6, .barAngleDeg = 0.0,
			.cutoff = 30000.0, .population = population("bulge")};
		galaxy.halo = {.enabled = true, .density = 0.5, .scale = 3000.0, .cutoff = 60000.0, .population = population("halo")};
		return galaxy;
	}
	if (name == "dwarf")
	{
		galaxy.bulge = {.enabled = true, .density = 1.6e7, .scale = 500.0, .axisY = 1.0, .axisZ = 0.7, .barAngleDeg = 0.0,
			.cutoff = 5000.0, .population = population("old-disk")};
		return galaxy;
	}
	return std::nullopt;
}

}
