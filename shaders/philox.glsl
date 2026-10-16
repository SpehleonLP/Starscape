// Philox4x32-10 (Salmon et al., SC11), the CPU's src/rng.h in GLSL, bit for bit.
uvec4 Philox4x32_10(uvec4 counter, uvec2 key)
{
	const uint kMultiplier0 = 0xD2511F53u;
	const uint kMultiplier1 = 0xCD9E8D57u;
	for (int round = 0; round < 10; ++round)
	{
		uint hi0;
		uint lo0;
		uint hi1;
		uint lo1;
		umulExtended(kMultiplier0, counter.x, hi0, lo0);
		umulExtended(kMultiplier1, counter.z, hi1, lo1);
		counter = uvec4(hi1 ^ counter.y ^ key.x, lo1, hi0 ^ counter.w ^ key.y, lo0);
		key += uvec2(0x9E3779B9u, 0xBB67AE85u);
	}
	return counter;
}

// The CPU's NextDouble, ((high >> 5)·2²⁶ + (low >> 6)) / 2⁵³, in float: no shaderFloat64 needed. Both scales
// are exact powers of two; the low term only matters when high >> 5 is 0.
float WordsToUnit(uint high, uint low)
{
	return float(high >> 5) * (1.0 / 134217728.0) + float(low >> 6) * (1.0 / 9007199254740992.0);
}
