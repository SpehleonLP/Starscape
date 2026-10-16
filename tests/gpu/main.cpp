#define DOCTEST_CONFIG_IMPLEMENT
#include "gpu/device.h"
#include <doctest.h>
#include <cstdio>

// ctest reads 77 as "skipped": a machine without a usable device must never report these tests as passed.
int main(int argc, char** argv)
{
	// The probe device is destroyed before the tests run: volkLoadDevice fills process-wide function pointers,
	// and a second device opened by a test would overwrite the first's while it is still alive.
	{
		auto device = starscape::gpu::OpenDevice({.allowSoftware = true});
		if (!device)
		{
			std::printf("SKIPPED: %s\n", device.error().message.c_str());
			return 77;
		}
	}
	doctest::Context context(argc, argv);
	return context.run();
}
