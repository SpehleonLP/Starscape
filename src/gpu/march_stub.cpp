// Compiled only when STARSCAPE_GPU is off; src/gpu/march.cpp is the real one.
#include "gpu/gpu_march.h"

namespace starscape
{

Expected<std::unique_ptr<PassMarcher>> OpenGpuMarcher(PreparedScene const&, MarchParams const&, uint32_t, std::function<void(double)>)
{
	return InputError("this build has no GPU march (built without glslangValidator)");
}

}
