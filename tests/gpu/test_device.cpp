#include "gpu/device.h"
#include <doctest.h>
#include <volk.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

using namespace starscape::gpu;

TEST_CASE("A device opens, names itself and round-trips a buffer")
{
	// Mutation: mapping the buffer without HOST_COHERENT memory makes the readback stale on some drivers; here,
	// a Submit that skips its fence wait returns before the copy and reddens the comparison.
	auto opened = OpenDevice({.allowSoftware = true});
	REQUIRE(opened.has_value());
	Device& device = **opened;
	std::printf("device: %s (driver %s, low priority %d)\n", device.Name().c_str(), device.Driver().c_str(), int(device.LowPriority()));
	CHECK_FALSE(device.Name().empty());
	auto a = device.CreateBuffer(1024, VK_BUFFER_USAGE_TRANSFER_SRC_BIT, true);
	auto b = device.CreateBuffer(1024, VK_BUFFER_USAGE_TRANSFER_DST_BIT, true);
	REQUIRE(a.has_value());
	REQUIRE(b.has_value());
	std::memset(a->mapped, 0x5A, 1024);
	std::memset(b->mapped, 0x00, 1024);
	auto elapsed = device.Submit([&](VkCommandBuffer cb)
	{
		VkBufferCopy const region{0, 0, 1024};
		vkCmdCopyBuffer(cb, a->buffer, b->buffer, 1, &region);
	});
	REQUIRE(elapsed.has_value());
	CHECK(static_cast<unsigned char const*>(b->mapped)[1023] == 0x5A);
	device.Destroy(*a);
	device.Destroy(*b);
}

TEST_CASE("STARSCAPE_GPU_DISABLE refuses every device and says so")
{
	// Mutation: never reading STARSCAPE_GPU_DISABLE (e.g. looking up another name) lets the device open, which
	// reddens the REQUIRE_FALSE.
	setenv("STARSCAPE_GPU_DISABLE", "1", 1);
	auto opened = OpenDevice({.allowSoftware = true});
	unsetenv("STARSCAPE_GPU_DISABLE");
	REQUIRE_FALSE(opened.has_value());
	CHECK(opened.error().message.find("STARSCAPE_GPU_DISABLE") != std::string::npos);
}

namespace
{

DeviceFacts Capable()
{
	return {.apiVersion = VK_API_VERSION_1_2, .computeQueue = true, .shaderInt64 = true, .scalarBlockLayout = true,
		.type = VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU};
}

}

TEST_CASE("A fully capable device has no shortfall")
{
	CHECK(Shortfall(Capable(), {}).empty());
}

TEST_CASE("A CPU-type device is refused unless software is allowed")
{
	// Mutation: deleting the CPU-type branch of Shortfall lets llvmpipe cook, reddening the first two checks.
	DeviceFacts facts = Capable();
	facts.type = VK_PHYSICAL_DEVICE_TYPE_CPU;
	CHECK_FALSE(Shortfall(facts, {}).empty());
	CHECK(Shortfall(facts, {}).find("CPU device") != std::string::npos);
	CHECK(Shortfall(facts, {.allowSoftware = true}).empty());
}

TEST_CASE("Each other requirement is named when missing")
{
	// Mutation: dropping any one branch of Shortfall reddens its case; each case misses only that one thing.
	DeviceFacts old = Capable();
	old.apiVersion = VK_API_VERSION_1_1;
	CHECK(Shortfall(old, {}).find("1.2") != std::string::npos);
	DeviceFacts noCompute = Capable();
	noCompute.computeQueue = false;
	CHECK(Shortfall(noCompute, {}).find("compute queue") != std::string::npos);
	DeviceFacts noInt64 = Capable();
	noInt64.shaderInt64 = false;
	CHECK(Shortfall(noInt64, {}).find("shaderInt64") != std::string::npos);
	DeviceFacts noScalar = Capable();
	noScalar.scalarBlockLayout = false;
	CHECK(Shortfall(noScalar, {}).find("scalarBlockLayout") != std::string::npos);
}
