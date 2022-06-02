#include "coulomb_kernel.h"
#include "coulomb_utils.h"

namespace Coulomb {

    void Run(Config* config);
    Config* normalizeConfig(Config* config);
    Config* copyConfig(Config* config);
    double f(Vector2d v, int specie, Config* config);
    Particle2d* initMarkers(int specie, Config* config, DistributionType type);
    Particle2d* initOutputPrintMesh(Config* config);
    double psi(Vector2d v, double eps);
    void mesh_distribution(
        double** ret,
        Particle2d* p_mesh,
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
    double nSpecie(
        Particle2d** p,
        int s,
        Config* config
    );
    double TemperatureSpecieSingleAxis(
        Particle2d** p,
        int s,
        int axis,
        Config* config
    );
    double TemperatureSpecie(
        Particle2d** p,
        int s,
        Config* config
    );
    Vector2d Momentum(Particle2d** p, Config* config);
    Vector2d MomentumSpecie(Particle2d** p, int s, Config* config);
    void printState(
        double** f_mesh,
        Particle2d* p_mesh,
        Particle2d** p1,
        Config* config,
        Config* config0,
        int t,
        double E0,
        Vector2d P0
    );
    double mccc_coefs_clog(int s1, int s2, Config* config);
    double coefs_nu(int s1, int s2, Config* config);
    void format_duration(int milliseconds, char* ret);
    Vector2d averageVelocitySpecie(
        Particle2d** p,
        int s,
        Config* config
    );
    double thermalizationTime(Config* config);
    double getTime();
}