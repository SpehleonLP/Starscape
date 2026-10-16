#include "gpu/gpu_march.h"

namespace starscape
{

Expected<std::unique_ptr<PassMarcher>> OpenGpuMarcher(PreparedScene const&, MarchParams const&, uint32_t, std::function<void(double)>)
{
#ifdef STARSCAPE_GPU
	return InputError("the GPU march is not built yet");
#else
	return InputError("this build has no GPU march (built without glslangValidator)");
#endif
}

}
