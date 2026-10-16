#include "gpu/device.h"
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <format>
#include <vector>

namespace starscape::gpu
{
namespace
{

Expected<void> Check(VkResult result, char const* what)
{
	if (result != VK_SUCCESS)
	{
		return InputError(std::format("Vulkan call {} failed: {}", what, int(result)));
	}
	return {};
}

struct Candidate
{
	VkPhysicalDevice device = VK_NULL_HANDLE;
	uint32_t         queueFamily = 0;
	std::string      lacks;   // empty when it qualifies
	int              rank = 0;
	std::string      name;
	uint32_t         vendor = 0;
	uint32_t         driverVersion = 0;
	float            timestampPeriod = 1.0f;
};

int TypeRank(VkPhysicalDeviceType type)
{
	switch (type)
	{
	case VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU:
		return 4;
	case VK_PHYSICAL_DEVICE_TYPE_INTEGRATED_GPU:
		return 3;
	case VK_PHYSICAL_DEVICE_TYPE_VIRTUAL_GPU:
		return 2;
	case VK_PHYSICAL_DEVICE_TYPE_CPU:
		return 1;
	default:
		return 0;
	}
}

}

std::string Shortfall(DeviceFacts const& facts, DeviceOptions const& options)
{
	if (facts.apiVersion < VK_API_VERSION_1_2)
	{
		return "Vulkan 1.2";
	}
	if (!facts.computeQueue)
	{
		return "a compute queue";
	}
	if (!facts.shaderInt64)
	{
		return "shaderInt64";
	}
	if (!facts.scalarBlockLayout)
	{
		return "scalarBlockLayout";
	}
	if (facts.type == VK_PHYSICAL_DEVICE_TYPE_CPU && !options.allowSoftware)
	{
		return "a GPU (it is a CPU device, not used for cooks)";
	}
	return {};
}

namespace
{

Candidate Examine(VkPhysicalDevice device, DeviceOptions const& options)
{
	Candidate c;
	c.device = device;
	VkPhysicalDeviceProperties properties;
	vkGetPhysicalDeviceProperties(device, &properties);
	c.name = properties.deviceName;
	c.rank = TypeRank(properties.deviceType);
	c.vendor = properties.vendorID;
	c.driverVersion = properties.driverVersion;
	c.timestampPeriod = properties.limits.timestampPeriod;
	DeviceFacts facts;
	facts.apiVersion = properties.apiVersion;
	facts.type = properties.deviceType;
	uint32_t families = 0;
	vkGetPhysicalDeviceQueueFamilyProperties(device, &families, nullptr);
	std::vector<VkQueueFamilyProperties> family(families);
	vkGetPhysicalDeviceQueueFamilyProperties(device, &families, family.data());
	auto const compute = std::find_if(family.begin(), family.end(), [](VkQueueFamilyProperties const& f) { return (f.queueFlags & VK_QUEUE_COMPUTE_BIT) != 0; });
	facts.computeQueue = compute != family.end();
	// A device below 1.2 need not implement vkGetPhysicalDeviceFeatures2, so it is not asked.
	if (properties.apiVersion >= VK_API_VERSION_1_2)
	{
		VkPhysicalDeviceVulkan12Features features12{.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_2_FEATURES};
		VkPhysicalDeviceFeatures2 features{.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2, .pNext = &features12};
		vkGetPhysicalDeviceFeatures2(device, &features);
		facts.shaderInt64 = features.features.shaderInt64 != VK_FALSE;
		facts.scalarBlockLayout = features12.scalarBlockLayout != VK_FALSE;
	}
	c.lacks = Shortfall(facts, options);
	c.queueFamily = compute == family.end() ? 0 : uint32_t(compute - family.begin());
	return c;
}

std::string DriverString(Candidate const& c)
{
	uint32_t const v = c.driverVersion;
	if (c.vendor == 0x10DE)
	{
		// NVIDIA packs its driver version as 10.8.8.6 bits rather than the Vulkan layout.
		return std::format("{}.{}.{}", v >> 22, (v >> 14) & 0xFF, (v >> 6) & 0xFF);
	}
	return std::format("{}.{}.{}", VK_VERSION_MAJOR(v), VK_VERSION_MINOR(v), VK_VERSION_PATCH(v));
}

bool HasExtension(std::vector<VkExtensionProperties> const& list, char const* name)
{
	return std::any_of(list.begin(), list.end(), [&](VkExtensionProperties const& e) { return std::strcmp(e.extensionName, name) == 0; });
}

}

Device::~Device()
{
	if (m_device != VK_NULL_HANDLE)
	{
		vkDeviceWaitIdle(m_device);
		if (m_pool != VK_NULL_HANDLE)
		{
			vkDestroyCommandPool(m_device, m_pool, nullptr);
		}
		if (m_fence != VK_NULL_HANDLE)
		{
			vkDestroyFence(m_device, m_fence, nullptr);
		}
		vkDestroyDevice(m_device, nullptr);
	}
	if (m_instance != VK_NULL_HANDLE)
	{
		vkDestroyInstance(m_instance, nullptr);
	}
}

Expected<std::unique_ptr<Device>> OpenDevice(DeviceOptions const& options)
{
	if (std::getenv("STARSCAPE_GPU_DISABLE") != nullptr)
	{
		return InputError("the GPU is disabled by STARSCAPE_GPU_DISABLE");
	}
	if (volkInitialize() != VK_SUCCESS)
	{
		return InputError("no Vulkan loader on this machine");
	}

	std::unique_ptr<Device> device(new Device());
	VkApplicationInfo const application{.sType = VK_STRUCTURE_TYPE_APPLICATION_INFO, .pApplicationName = "starscape", .apiVersion = VK_API_VERSION_1_2};
	VkInstanceCreateInfo const instanceInfo{.sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO, .pApplicationInfo = &application};
	if (vkCreateInstance(&instanceInfo, nullptr, &device->m_instance) != VK_SUCCESS)
	{
		device->m_instance = VK_NULL_HANDLE;
		return InputError("the Vulkan loader cannot create a 1.2 instance");
	}
	volkLoadInstance(device->m_instance);

	uint32_t count = 0;
	STARSCAPE_TRY(Check(vkEnumeratePhysicalDevices(device->m_instance, &count, nullptr), "vkEnumeratePhysicalDevices"));
	std::vector<VkPhysicalDevice> physical(count);
	STARSCAPE_TRY(Check(vkEnumeratePhysicalDevices(device->m_instance, &count, physical.data()), "vkEnumeratePhysicalDevices"));
	if (physical.empty())
	{
		return InputError("no Vulkan device found");
	}

	std::vector<Candidate> candidates;
	for (VkPhysicalDevice const handle : physical)
	{
		candidates.push_back(Examine(handle, options));
	}
	Candidate const* best = nullptr;
	for (Candidate const& c : candidates)
	{
		if (c.lacks.empty() && (best == nullptr || c.rank > best->rank))
		{
			best = &c;
		}
	}
	if (best == nullptr)
	{
		std::string message = "no Vulkan device qualifies:";
		for (Candidate const& c : candidates)
		{
			message += std::format(" {} lacks {};", c.name, c.lacks);
		}
		return InputError(std::move(message));
	}

	uint32_t extensionCount = 0;
	vkEnumerateDeviceExtensionProperties(best->device, nullptr, &extensionCount, nullptr);
	std::vector<VkExtensionProperties> extensions(extensionCount);
	vkEnumerateDeviceExtensionProperties(best->device, nullptr, &extensionCount, extensions.data());
	char const* priorityExtension = nullptr;
	if (HasExtension(extensions, VK_KHR_GLOBAL_PRIORITY_EXTENSION_NAME))
	{
		priorityExtension = VK_KHR_GLOBAL_PRIORITY_EXTENSION_NAME;
	}
	else if (HasExtension(extensions, VK_EXT_GLOBAL_PRIORITY_EXTENSION_NAME))
	{
		priorityExtension = VK_EXT_GLOBAL_PRIORITY_EXTENSION_NAME;
	}

	// The cook shares the machine with the owner's desktop, so it asks the driver for the LOW global priority.
	// A driver that lists the extension may still not offer LOW for this family, and another failure code
	// would then make a usable GPU look unusable; so the request is chained only when the family lists LOW.
	// The driver may also refuse it (documented as NOT_PERMITTED); then the cook runs at the default.
	bool offersLow = false;
	if (priorityExtension != nullptr)
	{
		uint32_t familyCount = 0;
		vkGetPhysicalDeviceQueueFamilyProperties2(best->device, &familyCount, nullptr);
		std::vector<VkQueueFamilyGlobalPriorityPropertiesKHR> priorities(familyCount,
			VkQueueFamilyGlobalPriorityPropertiesKHR{.sType = VK_STRUCTURE_TYPE_QUEUE_FAMILY_GLOBAL_PRIORITY_PROPERTIES_KHR});
		std::vector<VkQueueFamilyProperties2> families(familyCount, VkQueueFamilyProperties2{.sType = VK_STRUCTURE_TYPE_QUEUE_FAMILY_PROPERTIES_2});
		for (uint32_t i = 0; i < familyCount; ++i)
		{
			families[i].pNext = &priorities[i];
		}
		vkGetPhysicalDeviceQueueFamilyProperties2(best->device, &familyCount, families.data());
		VkQueueFamilyGlobalPriorityPropertiesKHR const& chosen = priorities[best->queueFamily];
		offersLow = std::find(chosen.priorities, chosen.priorities + chosen.priorityCount, VK_QUEUE_GLOBAL_PRIORITY_LOW_KHR)
			!= chosen.priorities + chosen.priorityCount;
	}
	// The priority itself comes from the global-priority struct below; this per-queue value only orders queues
	// within the family, and this device makes one.
	float const queuePriority = 1.0f;
	VkResult created = VK_ERROR_INITIALIZATION_FAILED;
	bool lowPriority = false;
	for (bool tryLow : {offersLow, false})
	{
		VkDeviceQueueGlobalPriorityCreateInfoKHR const priority{.sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_GLOBAL_PRIORITY_CREATE_INFO_KHR,
			.globalPriority = VK_QUEUE_GLOBAL_PRIORITY_LOW_KHR};
		VkDeviceQueueCreateInfo const queueInfo{.sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO, .pNext = tryLow ? &priority : nullptr,
			.queueFamilyIndex = best->queueFamily, .queueCount = 1, .pQueuePriorities = &queuePriority};
		VkPhysicalDeviceVulkan12Features features12{.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_2_FEATURES, .scalarBlockLayout = VK_TRUE};
		VkPhysicalDeviceFeatures2 features{.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2, .pNext = &features12};
		features.features.shaderInt64 = VK_TRUE;
		VkDeviceCreateInfo const deviceInfo{.sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO, .pNext = &features,
			.queueCreateInfoCount = 1, .pQueueCreateInfos = &queueInfo,
			.enabledExtensionCount = tryLow ? 1u : 0u, .ppEnabledExtensionNames = tryLow ? &priorityExtension : nullptr};
		created = vkCreateDevice(best->device, &deviceInfo, nullptr, &device->m_device);
		lowPriority = tryLow;
		if (created == VK_SUCCESS || created != VK_ERROR_NOT_PERMITTED_KHR || !tryLow)
		{
			break;
		}
	}
	device->m_lowPriority = created == VK_SUCCESS && lowPriority;
	if (created != VK_SUCCESS)
	{
		device->m_device = VK_NULL_HANDLE;
		return InputError(std::format("Vulkan call vkCreateDevice failed: {}", int(created)));
	}
	volkLoadDevice(device->m_device);

	device->m_physical = best->device;
	device->m_queueFamily = best->queueFamily;
	device->m_timestampPeriod = best->timestampPeriod;
	device->m_name = best->name;
	device->m_driver = DriverString(*best);
	vkGetPhysicalDeviceMemoryProperties(best->device, &device->m_memory);
	vkGetDeviceQueue(device->m_device, best->queueFamily, 0, &device->m_queue);

	VkCommandPoolCreateInfo const poolInfo{.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO, .flags = VK_COMMAND_POOL_CREATE_TRANSIENT_BIT,
		.queueFamilyIndex = best->queueFamily};
	STARSCAPE_TRY(Check(vkCreateCommandPool(device->m_device, &poolInfo, nullptr, &device->m_pool), "vkCreateCommandPool"));
	VkFenceCreateInfo const fenceInfo{.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO};
	STARSCAPE_TRY(Check(vkCreateFence(device->m_device, &fenceInfo, nullptr, &device->m_fence), "vkCreateFence"));
	return device;
}

Expected<Buffer> Device::CreateBuffer(VkDeviceSize size, VkBufferUsageFlags usage, bool hostVisible)
{
	Buffer buffer;
	buffer.size = size;
	VkBufferCreateInfo const info{.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO, .size = size, .usage = usage, .sharingMode = VK_SHARING_MODE_EXCLUSIVE};
	STARSCAPE_TRY(Check(vkCreateBuffer(m_device, &info, nullptr, &buffer.buffer), "vkCreateBuffer"));
	VkMemoryRequirements requirements;
	vkGetBufferMemoryRequirements(m_device, buffer.buffer, &requirements);

	// Coherent memory spares every mapped buffer its flush and invalidate calls.
	VkMemoryPropertyFlags const wanted = hostVisible ? (VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT)
		: VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT;
	uint32_t type = m_memory.memoryTypeCount;
	for (uint32_t i = 0; i < m_memory.memoryTypeCount; ++i)
	{
		if ((requirements.memoryTypeBits & (1u << i)) != 0 && (m_memory.memoryTypes[i].propertyFlags & wanted) == wanted)
		{
			type = i;
			break;
		}
	}
	if (type == m_memory.memoryTypeCount)
	{
		vkDestroyBuffer(m_device, buffer.buffer, nullptr);
		return InputError(std::format("{} has no {} memory type for a {} byte buffer", m_name, hostVisible ? "host-visible coherent" : "device-local", size));
	}
	VkMemoryAllocateInfo const allocation{.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO, .allocationSize = requirements.size, .memoryTypeIndex = type};
	auto fail = [&](VkResult result, char const* what) -> Expected<Buffer>
	{
		Destroy(buffer);
		return std::unexpected(Check(result, what).error());
	};
	if (VkResult const r = vkAllocateMemory(m_device, &allocation, nullptr, &buffer.memory); r != VK_SUCCESS)
	{
		buffer.memory = VK_NULL_HANDLE;
		return fail(r, "vkAllocateMemory");
	}
	if (VkResult const r = vkBindBufferMemory(m_device, buffer.buffer, buffer.memory, 0); r != VK_SUCCESS)
	{
		return fail(r, "vkBindBufferMemory");
	}
	if (hostVisible)
	{
		if (VkResult const r = vkMapMemory(m_device, buffer.memory, 0, VK_WHOLE_SIZE, 0, &buffer.mapped); r != VK_SUCCESS)
		{
			buffer.mapped = nullptr;
			return fail(r, "vkMapMemory");
		}
	}
	return buffer;
}

void Device::Destroy(Buffer& buffer)
{
	if (buffer.mapped != nullptr)
	{
		vkUnmapMemory(m_device, buffer.memory);
	}
	if (buffer.buffer != VK_NULL_HANDLE)
	{
		vkDestroyBuffer(m_device, buffer.buffer, nullptr);
	}
	if (buffer.memory != VK_NULL_HANDLE)
	{
		vkFreeMemory(m_device, buffer.memory, nullptr);
	}
	buffer = Buffer{};
}

Expected<double> Device::Submit(std::function<void(VkCommandBuffer)> const& record)
{
	auto const start = std::chrono::steady_clock::now();
	VkCommandBufferAllocateInfo const allocation{.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO, .commandPool = m_pool,
		.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY, .commandBufferCount = 1};
	VkCommandBuffer commands = VK_NULL_HANDLE;
	STARSCAPE_TRY(Check(vkAllocateCommandBuffers(m_device, &allocation, &commands), "vkAllocateCommandBuffers"));
	auto run = [&]() -> Expected<void>
	{
		VkCommandBufferBeginInfo const begin{.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO, .flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT};
		STARSCAPE_TRY(Check(vkBeginCommandBuffer(commands, &begin), "vkBeginCommandBuffer"));
		record(commands);
		STARSCAPE_TRY(Check(vkEndCommandBuffer(commands), "vkEndCommandBuffer"));
		VkSubmitInfo const submit{.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO, .commandBufferCount = 1, .pCommandBuffers = &commands};
		STARSCAPE_TRY(Check(vkQueueSubmit(m_queue, 1, &submit, m_fence), "vkQueueSubmit"));
		Expected<void> const waited = Check(vkWaitForFences(m_device, 1, &m_fence, VK_TRUE, UINT64_MAX), "vkWaitForFences");
		// The fence is reset even when the wait failed, so a later submission does not inherit a signal.
		vkResetFences(m_device, 1, &m_fence);
		return waited;
	};
	Expected<void> const result = run();
	vkFreeCommandBuffers(m_device, m_pool, 1, &commands);
	if (!result)
	{
		return std::unexpected(result.error());
	}
	return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

Expected<VkPipeline> Device::CreateComputePipeline(std::span<uint32_t const> spirv, VkPipelineLayout layout)
{
	VkShaderModuleCreateInfo const moduleInfo{.sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO, .codeSize = spirv.size_bytes(), .pCode = spirv.data()};
	VkShaderModule module = VK_NULL_HANDLE;
	STARSCAPE_TRY(Check(vkCreateShaderModule(m_device, &moduleInfo, nullptr, &module), "vkCreateShaderModule"));
	VkComputePipelineCreateInfo const info{.sType = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO,
		.stage = {.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO, .stage = VK_SHADER_STAGE_COMPUTE_BIT, .module = module, .pName = "main"},
		.layout = layout};
	VkPipeline pipeline = VK_NULL_HANDLE;
	VkResult const result = vkCreateComputePipelines(m_device, VK_NULL_HANDLE, 1, &info, nullptr, &pipeline);
	vkDestroyShaderModule(m_device, module, nullptr);
	STARSCAPE_TRY(Check(result, "vkCreateComputePipelines"));
	return pipeline;
}

}
