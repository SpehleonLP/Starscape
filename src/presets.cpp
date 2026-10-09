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

}
