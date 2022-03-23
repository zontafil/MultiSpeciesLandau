#include "coulomb_utils.h"

#if defined CUDA && defined INTELLISENSE
#include "cuda_runtime.h"
#include "crt/device_functions.h"
#endif

namespace Kernel {
    void computedSdv(
        VectorXd* ret,
        Particle2d* p,
        Config* config
    );
}