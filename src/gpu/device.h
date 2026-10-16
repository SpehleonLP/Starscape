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

// What a physical device offers, as far as qualifying it goes.
struct DeviceFacts
{
	uint32_t             apiVersion = 0;
	bool                 computeQueue = false;
	bool                 shaderInt64 = false;
	bool                 scalarBlockLayout = false;
	VkPhysicalDeviceType type = VK_PHYSICAL_DEVICE_TYPE_OTHER;
};

// The first requirement `facts` misses, worded to follow "lacks"; empty when the device qualifies.
std::string Shortfall(DeviceFacts const& facts, DeviceOptions const& options);

// The largest buffers a device takes. Vulkan guarantees only 2^27 bytes of storage-buffer range, and drivers
// report anything from 128 MiB up, while the face sums alone are 453 MB at face 2048.
struct BufferLimits
{
	VkDeviceSize storageRange = 0;   // maxStorageBufferRange: the most one storage-buffer binding can see
	VkDeviceSize allocation = 0;     // maxMemoryAllocationSize: the most one allocation can hold
};

// The limit a storage buffer of `size` bytes exceeds, worded to follow "over"; empty when it fits.
std::string BufferShortfall(VkDeviceSize size, BufferLimits const& limits);

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
	// True only when the queue was created with the LOW global priority request.
	bool LowPriority() const { return m_lowPriority; }
	BufferLimits const& Limits() const { return m_limits; }
	// Nanoseconds per timestamp tick, for shaders' timestamp queries.
	float TimestampPeriod() const { return m_timestampPeriod; }

	// Host-visible buffers are HOST_COHERENT and stay mapped; the others are DEVICE_LOCAL.
	Expected<Buffer> CreateBuffer(VkDeviceSize size, VkBufferUsageFlags usage, bool hostVisible);
	void Destroy(Buffer& buffer);
	// Owners call this before destroying objects a submission may use: a failed Submit can unwind with work
	// still on the queue, and ~Device's own wait comes after its users' destructors.
	void WaitIdle();
	// Records `record` into a one-shot command buffer, submits it and waits; returns the wall time in seconds.
	Expected<double> Submit(std::function<void(VkCommandBuffer)> const& record);
	// Entry point "main"; the shader module is not kept.
	Expected<VkPipeline> CreateComputePipeline(std::span<uint32_t const> spirv, VkPipelineLayout layout);
	void DestroyPipeline(VkPipeline& pipeline);

	// Compute-stage descriptors: one set of `bindings` storage buffers at bindings 0..bindings-1, and
	// `pushBytes` of push constants (0: none). Each Create has its Destroy; destroying a pool frees its sets.
	Expected<VkDescriptorSetLayout> CreateDescriptorSetLayout(uint32_t bindings);
	void DestroyDescriptorSetLayout(VkDescriptorSetLayout& layout);
	Expected<VkPipelineLayout> CreatePipelineLayout(VkDescriptorSetLayout setLayout, uint32_t pushBytes);
	void DestroyPipelineLayout(VkPipelineLayout& layout);
	Expected<VkDescriptorPool> CreateDescriptorPool(uint32_t sets, uint32_t storageBuffers);
	void DestroyDescriptorPool(VkDescriptorPool& pool);
	Expected<VkDescriptorSet> AllocateDescriptorSet(VkDescriptorPool pool, VkDescriptorSetLayout layout);
	// Binds buffers[i] whole to binding i.
	void WriteStorageBuffers(VkDescriptorSet set, std::span<Buffer const> buffers);

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
	bool             m_lowPriority = false;
	BufferLimits     m_limits;
	VkPhysicalDeviceMemoryProperties m_memory{};
	std::string      m_name;
	std::string      m_driver;
};

// InputError naming the unmet requirement when no device qualifies.
Expected<std::unique_ptr<Device>> OpenDevice(DeviceOptions const& options);

}
