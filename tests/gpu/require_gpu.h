#pragma once
#include <cstdio>
#include <cstdlib>
#include <string>

namespace starscape::test
{

// The exit for a test binary that found no device: 77, which ctest reads as "skipped", or 1 when
// STARSCAPE_REQUIRE_GPU is set (to anything but empty or 0). The owner's gate sets it, because a device check that
// wrongly refuses every GPU would otherwise skip the whole GPU suite and still report 100% passed.
inline int NoDeviceExit(std::string const& why)
{
	char const* const require = std::getenv("STARSCAPE_REQUIRE_GPU");
	if (require != nullptr && std::string(require) != "" && std::string(require) != "0")
	{
		std::printf("FAILED: STARSCAPE_REQUIRE_GPU is set and no device qualified: %s\n", why.c_str());
		return 1;
	}
	std::printf("SKIPPED: %s\n", why.c_str());
	return 77;
}

}
