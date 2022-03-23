#include "coulombStructurePreserving.h"

using std::setprecision;

Config buildConfig() {
    Config config;

    double L = 10;
    int MARKERS_PER_DIM = 60;

    // peaks of the double Maxwellian
    config.u1 = Vector2d(-2, 1);
    config.u2 = Vector2d(0, -1);

    config.dx = 1E-4; // dx for finite difference derivative
    config.dt = 1/16.;
    config.n_timesteps = 1; // number of timesteps
    config.newtonTolerance = 1E-14; // minimum target for error of eq. of motions
    config.useNewton = 0;
    config.maxEOMIterations = 20;
    config.nu = 1;
    config.m = 1;
    config.h = 2.*L/MARKERS_PER_DIM;
    config.eps = 0.64*pow(config.h, 1.98);
    config.xmin = -L;
    config.xmax = L;
    config.ymin = -L;
    config.ymax = L;
    config.nx = MARKERS_PER_DIM;
    config.ny = MARKERS_PER_DIM;
    config.recordAtStep = max(config.n_timesteps - 1, 1);

    // double kHermite[5] = {-2.020183, -0.958572, 0.000000, 0.958572, 2.020183};
    // double wHermite[5] = {0.019953, 0.393619, 0.945309, 0.393619, 0.019953};
    double kHermite[6] = {-2.3506049736745, -1.3358490740137, -0.43607741192762, 0.43607741192762, 1.3358490740137, 2.3506049736745};
    double wHermite[6] = {0.0045300099055088, 0.15706732032286, 0.72462959522439, 0.72462959522439, 0.15706732032286, 0.0045300099055088};
    config.nHermite = 6;
    config.kHermite = new double[6];
    config.wHermite = new double[6];
    copy(kHermite, kHermite+config.nHermite, config.kHermite);
    copy(wHermite, wHermite+config.nHermite, config.wHermite);

    config.nmarkers = config.nx * config.ny;
    config.cudaThreadsPerBlock = 1024;

    return config;
}

int main() {
    Config config = buildConfig();
    Particle2d* p0 = new Particle2d[config.nmarkers];
    Particle2d* p1 = new Particle2d[config.nmarkers];
    initMarkers(p1, &config);

    // initial energy
    double E0 = K(p1, &config), E;
    Vector2d P0 = Momentum(p1, &config), P;
    VectorXd dSdV(2*config.nmarkers);

    print_out(VERBOSE_NORMAL, "Markers: %d dT: %f\n", config.nmarkers, config.dt);
    print_out(VERBOSE_NORMAL, "Initial Energy: %e\n", E0);
    print_out(VERBOSE_NORMAL, "Initial Momentum: %e %e\n", P0[0], P0[1]);
    for (int i=0; i<config.nmarkers; i++) {
        print_out(VERBOSE_DEBUG, "Init state: Marker ID: %d Vx: %.15e Vy: %.15e W: %.15e\n", i, p1[i].z[0], p1[i].z[1], p1[i].weight);
    }

    for (int t=0; t<config.n_timesteps; t++) {

        print_out(VERBOSE_NORMAL, "\nTimestep: %d\n", t);

        copy(p1, p1+config.nmarkers, p0);

        // precompute entropy gradient
        Kernel::computedSdv(&dSdV, p0, &config);
        if (VERBOSE_LEVEL >= VERBOSE_SILLY) {
            cout << "==== dSdV" << endl;
            cout << dSdV << endl;
            cout << "==== dSdV end" << endl;
        }

        // fixed point newton iterations
        for (int j=0; j<config.maxEOMIterations; j++) {
            if (config.useNewton) {
                if (pushForward_iteration(p0, p1, &dSdV, &config)) {
                    break;
                }
            } else {
                if (pushForward_dv(p0, p1, &dSdV, &config)) {
                    break;
                }
            }
        }

        // print system state and debug
        if (t%config.recordAtStep == 0) {

            E = K(p1, &config);
            print_out(VERBOSE_NORMAL, "Energy: %.15e Error: %.15e\n", E, (E-E0)/E0);

            P = Momentum(p1, &config);
            print_out(VERBOSE_NORMAL, "Momentum: %.15e %.15e\n", P[0], P[1]);
            print_out(VERBOSE_NORMAL, "Momentum Error: %.15e %.15e\n", (P[0]-P0[0])/P0[0], (P[1]-P0[1])/P0[1]);

            if (VERBOSE_LEVEL >= VERBOSE_DEBUG) {
                for (int i=0; i<config.nmarkers; i++) {
                    print_out(VERBOSE_DEBUG, "Marker ID: %d Vx: %.15e Vy: %.15e W: %.15e\n", i, p1[i].z[0], p1[i].z[1], p1[i].weight);
                }
            }
        }
    }

    free(p0);
    free(p1);
    return 0;
}

/**
 * @brief Distribution function of initial states
 * 
 * @param v 
 * @param u1 
 * @param u2 
 * @return double 
 */
double f(Vector2d v, Config* config) {
    return pow(config->h,2)/(4.*CONST_PI) * (exp(-(v-config->u1).squaredNorm()/2.) + exp(-(v-config->u2).squaredNorm()/2.));
}

/**
 * @brief Radial basis function
 * 
 * @param v 
 * @param eps 
 * @return double 
 */
double psi(Vector2d v, double eps) {
    return exp(-(pow(v[0],2) + pow(v[1],2)) /2./eps) / CONST_2PI / eps;
}

/**
 * @brief Init markers in a 2D mesh
 * 
 * @param p 
 * @param config 
 */
void initMarkers(Particle2d* p, Config* config) {
    int idx;
    for (int i = 0; i<config->ny; i++) {
        for (int j = 0; j<config->nx; j++) {
            idx = i*config->nx + j;
            p[idx].z[0] = double(j+0.5) / (config->nx) * (config->xmax-config->xmin) + config->xmin;
            p[idx].z[1] = double(i+0.5) / (config->nx) * (config->ymax-config->ymin) + config->ymin;
            p[idx].weight = f(p[idx].z, config);
        }
    }
}

int pushForward_dv(
    Particle2d* p0,
    Particle2d* p1,
    VectorXd* dSdV,
    Config* config
) {
    VectorXd f(2*config->nmarkers);
    Kernel::f_eqmotion_dv(&f, p0, p1, dSdV, config);
    if (VERBOSE_LEVEL >= VERBOSE_SILLY) {
        printf("=== dv\n");
        cout << f << endl;
        printf("=== dv end\n");
    }

    double znew, err = 0;
    for (int i=0; i<config->nmarkers; i++) {
        for (int j=0; j<2; j++) {
            znew = p0[i].z[j] + config->dt * f(i*2+j);
            err += pow(p1[i].z[j] - znew, 2);
            p1[i].z[j] = znew;
        }
    }

    print_out(VERBOSE_DEBUG, "Eq. of motions precision: %e\n", sqrt(err));

    if (sqrt(err) < config->newtonTolerance) {
        return 1;
    }

    return 0;
}

/**
 * @brief Markers push forward using a fixed point iteration with Newton method.
 * 
 * @param p0 
 * @param p1 
 * @param config 
 * @return int 1 if the equations of motion norm is below minimum tolerance
 */
int pushForward_iteration(
    Particle2d* p0,
    Particle2d* p1,
    VectorXd* dSdV,
    Config* config
) {
    VectorXd f = VectorXd::Zero(2*config->nmarkers);
    VectorXd f1 = VectorXd::Zero(2*config->nmarkers);
    VectorXd dv = VectorXd::Zero(2*config->nmarkers);
    MatrixXd Jf = MatrixXd::Zero(2*config->nmarkers, 2*config->nmarkers);

    f_eqmotion(&f, p0, p1, dSdV, config);
    double fnorm = f.norm();
    print_out(VERBOSE_DEBUG, "Eq. of motions precision: %e max: %e\n", fnorm, f.maxCoeff());
    if (fnorm < config->newtonTolerance) {
        return 1;
    }

    for (int i=0;i<config->nmarkers; i++){
        for (int j=0; j<2; j++) {
            p1[i].z[j] += config->dx;

            f_eqmotion(&f1, p0, p1, dSdV, config);

            Jf.col(i*2+j) = (f1 - f)/config->dx;

            p1[i].z[j] -= config->dx;
        }
    }

    dv = -Jf.inverse()*f;

    // copy back from dv to the markers array
    for (int i=0;i<config->nmarkers; i++){
        for (int j=0; j<2; j++) {
            p1[i].z[j] += dv[i*2+j];
        }
    }

    return 0;
}

/**
 * @brief Equations of motion for marker trajectories: f(z0, z1) = 0
 * where z0, z1 are two consecutive time steps
 * 
 * @param f 
 * @param p0 
 * @param p1 
 * @param config 
 */
void f_eqmotion(
    VectorXd* f,
    Particle2d* p0,
    Particle2d* p1,
    VectorXd* dSdV,
    Config* config
) {

    Kernel::f_eqmotion_dv(f, p0, p1, dSdV, config);

    for (int i=0; i<config->nmarkers; i++) {
        for (int j=0; j<2; j++) {
            f->coeffRef(i*2+j) = (p1[i].z[j] - p0[i].z[j]) / config->dt - f->coeffRef(i*2+j);
        }
    }

    if (VERBOSE_LEVEL >= VERBOSE_SILLY) {
        printf("=== fcoeff\n");
        for (int i=0; i<config->nmarkers; i++) {
            printf("%e %e\n", f->coeffRef(2*i), f->coeffRef(2*i+1));
        }
        printf("=== fcoeff end\n");
    }
}

/**
 * @brief Compute the kinetic energy of the system
 * 
 * @param p 
 * @param config 
 * @return double 
 */
double K(Particle2d* p, Config* config) {
    double ret = 0;
    for (int i = 0; i<config->nmarkers; i++) {
        ret += p[i].weight * config->m * 0.5 * p[i].z.squaredNorm();
    }
    return ret;
}

/**
 * @brief Compute momentum of the system
 * 
 * @param p 
 * @param config 
 * @return Vector2d 
 */
Vector2d Momentum(Particle2d* p, Config* config) {
    Vector2d ret(0,0);
    for (int i=0; i<config->nmarkers; i++) {
        ret += p[i].weight * config->m * p[i].z;
    }
    return ret;
}