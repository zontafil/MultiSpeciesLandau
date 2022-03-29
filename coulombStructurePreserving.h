#include "coulomb_kernel.h"
#include "coulomb_utils.h"

namespace Coulomb {

    void Run(Config* config);
    double f(Vector2d v, Config* config);
    void initMarkers(Particle2d* p, Config* config);
    double psi(Vector2d v, double eps);
    void mesh_distribution(
        double* ret,
        Particle2d* p_mesh,
        Particle2d* p,
        Config* config
    );
    int pushForward_dv(
        Particle2d* p0,
        Particle2d* p1,
        VectorXd* dSdV,
        Config* config
    );
    int pushForward_iteration(
        Particle2d* p0,
        Particle2d* p1,
        VectorXd* dSdV,
        Config* config
    );
    void f_eqmotion(
        VectorXd* f,
        Particle2d* p0,
        Particle2d* p1,
        VectorXd* dSdV,
        Config* config
    );
    double K(Particle2d* p, Config* config);
    Vector2d Momentum(Particle2d* p, Config* config);
}