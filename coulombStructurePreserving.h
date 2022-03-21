#include "coulomb_kernel.h"
#include "coulomb_utils.h"


double f(Vector2d v, Config* config);
void initMarkers(Particle2d* p, Config* config);
double psi(Vector2d v, double eps);
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
void f_eqmotion_dv(
    VectorXd* dv,
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
void buildQGamma(
    MatrixXd* ret,
    Particle2d* p0,
    Particle2d* p1,
    VectorXd* dSdV,
    Config* config
);
void Q(Matrix2d* ret, Vector2d v);
double K(Particle2d* p, Config* config);
Vector2d Momentum(Particle2d* p, Config* config);