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
