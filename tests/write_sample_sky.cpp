// Writes a small synthetic sky to argv[1], so the CLI's inspect and preview tests have a file without a GPU cook.
#include "writer.h"
#include "rng.h"
#include <cstdio>

using namespace starscape;

int main(int argc, char** argv)
{
	if (argc != 2)
	{
		std::fprintf(stderr, "usage: write-sample-sky <out.kre_starscape>\n");
		return 1;
	}
	constexpr uint32_t faceSize = 16;
	Cooked cooked;
	cooked.faceSize = faceSize;
	cooked.fluxThreshold = 0.004;
	cooked.passes = 1;
	RandomStream rng(2, Purpose::Test, 0, 0, 0);
	for (uint32_t i = 0; i < 6 * faceSize * faceSize; ++i)
	{
		cooked.diffuse.push_back({rng.NextDouble(), rng.NextDouble(), rng.NextDouble()});
	}
	cooked.stars.push_back({{0.0, 0.0, 1.0}, {3.0, 2.0, 1.0}});
	cooked.stars.push_back({{0.6, 0.8, 0.0}, {0.01, 0.02, 0.03}});
	cooked.stats.starCount = cooked.stars.size();
	Expected<void> const written = WriteFileAtomic(argv[1], BuildFile(cooked, 0));
	if (!written)
	{
		std::fprintf(stderr, "write-sample-sky: %s\n", written.error().message.c_str());
		return 2;
	}
	return 0;
}
