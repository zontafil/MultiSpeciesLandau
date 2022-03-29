#include "coulomb_utils.h"

#if defined CUDA && defined INTELLISENSE
#include "cuda_runtime.h"
#include "crt/device_functions.h"
#endif

using namespace Coulomb;

namespace Kernel {
    void computedSdv(
        VectorXd* ret,
        Particle2d* p,
        Config* config
    );
    void f_eqmotion_dv(
        VectorXd* dv,
        Particle2d* p0,
        Particle2d* p1,
        VectorXd* dSdV,
        Config* config
    );
}