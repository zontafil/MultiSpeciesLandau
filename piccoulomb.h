typedef double real;  /**< Double precision float   */
#include <math.h>
#include <stdio.h>
#include <iomanip>
#include <limits>
#include <stdio.h>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

#define CONST_PI     3.1415926535897932384626
#define CONST_2PI    6.2831853071795862319959
#define PI2E0_5 2.50662827463

enum {
    VERBOSE_DEBUG = 2,
    VERBOSE_NORMAL = 1,
    VERBOSE_MINIMAL = 0,
    VERBOSE_IO = 1
};

/**
 * @brief Print to standard output
 */
#define print_out(v,...) { if(VERBOSE_LEVEL >= (v)) printf(__VA_ARGS__); }

typedef struct {
    double xmin;
    double xmax;
    int nx;
    double ymin;
    double ymax;
    int ny;
    int nmarkers;
    double dt;
    double dx;
    double nu;
    double m;
    double h;
    double eps;
    int n_timesteps;
    int n_newton;
    Vector2d u1;
    Vector2d u2;
    double* kHermite;
    double* wHermite;
    int nHermite;
} Config;

typedef struct {
    double weight;
    Vector2d z;
} Particle2d;

double f(Vector2d v, Config* config);
void initMarkers(Particle2d* p, Config* config);
double psi(Vector2d v, double eps);
void pushForward_iteration(
    Particle2d* p0,
    Particle2d* p1,
    Config* config
);
void f_eqmotion(
    VectorXd* f,
    Particle2d* p0,
    Particle2d* p1,
    Config* config
);
void buildQGamma(
    MatrixXd* ret,
    Particle2d* p0,
    Particle2d* p1,
    Config* config
);
void computedSdv(
    VectorXd* ret,
    Particle2d* p,
    Config* config
);
void Q(Matrix2d* ret, Vector2d v);
double K(Particle2d* p, Config* config);