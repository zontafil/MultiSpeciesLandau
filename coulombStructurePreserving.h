#include "coulomb_kernel.h"
#include "coulomb_utils.h"

namespace Coulomb {

    void Run(Config* config);
    double f(Vector2d v, int specie, Config* config);
    Particle2d* initMarkers(int specie, Config* config, DistributionType type);
    double psi(Vector2d v, double eps);
    void mesh_distribution(
        double** ret,
        Particle2d** p_mesh,
        Particle2d** p,
        Config* config
    );
    int pushForward_dv(
        Particle2d** p0,
        Particle2d** p1,
        VectorXd* dSdV,
        Config* config
    );
    int pushForwardNewtonIteration(
        Particle2d** p0,
        Particle2d** p1,
        VectorXd* dSdV,
        Config* config
    );
    void f_eqmotion(
        VectorXd* f,
        Particle2d** p0,
        Particle2d** p1,
        VectorXd* dSdV,
        Config* config
    );
    double K(Particle2d** p, Config* config);
    double Kspecie(Particle2d** p, int s, Config* config);
    Vector2d Momentum(Particle2d** p, Config* config);
    Vector2d MomentumSpecie(Particle2d** p, int s, Config* config);
    void printState(
        double** f_mesh,
        Particle2d** p_mesh,
        Particle2d** p1,
        Config* config,
        int t,
        double E0,
        Vector2d P0
    );
}