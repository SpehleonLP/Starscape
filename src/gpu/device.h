#pragma once
#include "common.h"
#include <volk.h>
#include <cstdint>
#include <functional>
#include <memory>
#include <span>
#include <string>

namespace starscape::gpu
{

struct DeviceOptions
{
	bool allowSoftware = false;   // tests only: let a CPU-type device such as llvmpipe through
};

struct Buffer
{
	VkBuffer       buffer = VK_NULL_HANDLE;
	VkDeviceMemory memory = VK_NULL_HANDLE;
	VkDeviceSize   size = 0;
	void*          mapped = nullptr;   // host-visible buffers only
};

// A headless compute device: one queue, one command pool, one fence. Submissions are synchronous.
class Device
{
public:
	~Device();
	Device(Device const&) = delete;
	Device& operator=(Device const&) = delete;

	std::string const& Name() const { return m_name; }
	std::string const& Driver() const { return m_driver; }
	VkDevice Handle() const { return m_device; }
	// Nanoseconds per timestamp tick, for shaders' timestamp queries.
	float TimestampPeriod() const { return m_timestampPeriod; }

	// Host-visible buffers are HOST_COHERENT and stay mapped; the others are DEVICE_LOCAL.
	Expected<Buffer> CreateBuffer(VkDeviceSize size, VkBufferUsageFlags usage, bool hostVisible);
	void Destroy(Buffer& buffer);
	// Records `record` into a one-shot command buffer, submits it and waits; returns the wall time in seconds.
	Expected<double> Submit(std::function<void(VkCommandBuffer)> const& record);
	// Entry point "main"; the shader module is not kept.
	Expected<VkPipeline> CreateComputePipeline(std::span<uint32_t const> spirv, VkPipelineLayout layout);

private:
	Device() = default;
	friend Expected<std::unique_ptr<Device>> OpenDevice(DeviceOptions const& options);

	VkInstance       m_instance = VK_NULL_HANDLE;
	VkPhysicalDevice m_physical = VK_NULL_HANDLE;
	VkDevice         m_device = VK_NULL_HANDLE;
	VkQueue          m_queue = VK_NULL_HANDLE;
	VkCommandPool    m_pool = VK_NULL_HANDLE;
	VkFence          m_fence = VK_NULL_HANDLE;
	uint32_t         m_queueFamily = 0;
	float            m_timestampPeriod = 1.0f;
	VkPhysicalDeviceMemoryProperties m_memory{};
	std::string      m_name;
	std::string      m_driver;
};

// InputError naming the unmet requirement when no device qualifies.
Expected<std::unique_ptr<Device>> OpenDevice(DeviceOptions const& options);

}
