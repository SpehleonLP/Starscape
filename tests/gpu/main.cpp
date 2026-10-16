#define DOCTEST_CONFIG_IMPLEMENT
#include "gpu/device.h"
#include "gpu/require_gpu.h"
#include <doctest.h>

// ctest reads 77 as "skipped": a machine without a usable device must never report these tests as passed.
// STARSCAPE_REQUIRE_GPU turns that into a failure.
int main(int argc, char** argv)
{
	// The probe device is destroyed before the tests run: volkLoadDevice fills process-wide function pointers,
	// and a second device opened by a test would overwrite the first's while it is still alive.
	{
		auto device = starscape::gpu::OpenDevice({.allowSoftware = true});
		if (!device)
		{
			return starscape::test::NoDeviceExit(device.error().message);
		}
	}
	doctest::Context context(argc, argv);
	return context.run();
}
