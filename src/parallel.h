#pragma once
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace starscape
{

// Runs body(i) for every i in [0, count) on `threads` threads, the caller included. The schedule
// is not deterministic, so bodies must write only to per-index slots; callers reduce in index order.
inline void ParallelFor(uint32_t count, unsigned threads, std::function<void(uint32_t)> const& body,
	std::function<void(uint32_t done)> const& progress)
{
	std::atomic<uint32_t> next{0};
	std::atomic<uint32_t> done{0};
	std::mutex progressMutex;
	auto worker = [&]
	{
		for (uint32_t index = next++; index < count; index = next++)
		{
			body(index);
			uint32_t const finished = ++done;
			if (progress)
			{
				std::lock_guard lock(progressMutex);
				progress(finished);
			}
		}
	};
	std::vector<std::jthread> pool;
	for (unsigned t = 1; t < std::max(1u, threads); ++t)
	{
		pool.emplace_back(worker);
	}
	worker();
}

}
