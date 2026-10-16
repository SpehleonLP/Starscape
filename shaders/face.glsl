// Cube-face geometry and the ray jitter. Include after philox.glsl.

// kre_starscape::FaceToDir: the Vulkan cube face table, sky-frame directions.
vec3 FaceToDir(uint face, float s, float t)
{
	const float u = 2.0 * s - 1.0;
	const float v = 2.0 * t - 1.0;
	switch (face)
	{
	case 0u:
		return vec3(1.0, -v, -u);
	case 1u:
		return vec3(-1.0, -v, u);
	case 2u:
		return vec3(u, 1.0, v);
	case 3u:
		return vec3(u, -1.0, -v);
	case 4u:
		return vec3(u, -v, 1.0);
	default:
		return vec3(-u, -v, -1.0);
	}
}

float TexelSolidAngle(uint faceSize, uint x, uint y)
{
	const float stepSize = 2.0 / float(faceSize);
	const float u0 = -1.0 + float(x) * stepSize;
	const float v0 = -1.0 + float(y) * stepSize;
	const float u1 = u0 + stepSize;
	const float v1 = v0 + stepSize;
	return atan(u0 * v0, sqrt(u0 * u0 + v0 * v0 + 1.0)) - atan(u0 * v1, sqrt(u0 * u0 + v1 * v1 + 1.0))
		- atan(u1 * v0, sqrt(u1 * u1 + v0 * v0 + 1.0)) + atan(u1 * v1, sqrt(u1 * u1 + v1 * v1 + 1.0));
}

// The sky-frame direction of the ray through stratum (stratumX, stratumY) of texel (x, y) on `face`, jittered by the
// host's RandomStream(seed, Jitter, texel, 0, 0, stratum): purpose Jitter = 1 in the top byte of word 2, and the
// stratum index as the block, so a pass draws its stratum's sample without replaying the ones before it. The march
// and the jitter probe both call this, so the probe tests the march's own directions.
vec3 TexelJitterDirection(uint faceSize, uint face, uint x, uint y, uint stratumX, uint stratumY, uint perAxis, uvec2 seed)
{
	const uint texel = (face * faceSize + y) * faceSize + x;
	const uint stratum = stratumY * perAxis + stratumX;
	const uvec4 block = Philox4x32_10(uvec4(texel, 0u, 1u << 24, stratum), seed);
	const float u = WordsToUnit(block.x, block.y);
	const float v = WordsToUnit(block.z, block.w);
	const float s = (float(x) + (float(stratumX) + u) / float(perAxis)) / float(faceSize);
	const float t = (float(y) + (float(stratumY) + v) / float(perAxis)) / float(faceSize);
	return normalize(FaceToDir(face, s, t));
}
