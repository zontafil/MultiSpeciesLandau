#include "piccoulomb.h"

using std::setprecision;

Config buildConfig() {
    Config config;

    double L = 10;
    int MARKERS_PER_DIM = 60;

    // peaks of the double Maxwellian
    Vector2d u1, u2;
    u1 << -2, 1;
    u2 << 0, -1;
    config.u1 = u1;
    config.u2 = u2;

    config.dx = 1E-4; // dx for finite difference derivative
    config.dt = 1./16.; // time step
    config.n_timesteps = 2; // number of timesteps
    config.n_newton = 4; // number of Newton iterations for particle push ofrward
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

    double kHermite[5] = {-2.856970, -1.355626, 0.000000, 1.355626, 2.856970};
    double wHermite[5] = {0.028218, 0.556662, 1.336868, 0.556662, 0.028218};
    config.nHermite = 5;
    config.kHermite = new double[5];
    config.wHermite = new double[5];
    copy(kHermite, kHermite+config.nHermite, config.kHermite);
    copy(wHermite, wHermite+config.nHermite, config.wHermite);

    config.nmarkers = config.nx * config.ny;

    return config;
}

int main() {
    Config config = buildConfig();
    Particle2d* p0 = new Particle2d[config.nmarkers];
    Particle2d* p1 = new Particle2d[config.nmarkers];
    initMarkers(p1, &config);

    // initial energy
    double E0 = K(p1, &config), E;
    print_out(VERBOSE_NORMAL, "Initial Energy: %e\n", E0);
    
    for (int i=0; i<config.nmarkers; i++) {
        print_out(VERBOSE_DEBUG, "Init State, ID: %d Vx: %.15e Vy: %.15e W: %.15e\n", i, p1[i].z[0], p1[i].z[1], p1[i].weight);
    }

    for (int t=0; t<config.n_timesteps; t++) {
        print_out(VERBOSE_NORMAL, "Timestep: %d\n", t);

        // print coordinates, energy error
        E = K(p1, &config);
        print_out(VERBOSE_NORMAL, "Energy: %e Error: %e\n", E, (E-E0)/E0);

        copy(p1, p1+config.nmarkers, p0);

        for (int j=0; j<config.n_newton; j++) {
            pushForward_iteration(p0, p1, &config);
        }
    }

    for (int i=0; i<config.nmarkers; i++) {
        print_out(VERBOSE_DEBUG, "End State, ID: %d Vx: %.15e Vy: %.15e W: %.15e\n", i, p1[i].z[0], p1[i].z[1], p1[i].weight);
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
    for (int i = 0; i<config->nx; i++) {
        for (int j = 0; j<config->ny; j++) {
            idx = i*config->ny + j;
            p[idx].z[0] = double(i+0.5) / (config->nx) * (config->xmax-config->xmin) + config->xmin;
            p[idx].z[1] = double(j+0.5) / (config->nx) * (config->ymax-config->ymin) + config->ymin;
            p[idx].weight = f(p[idx].z, config);
        }
    }
}

/**
 * @brief Markers push forward using a fixed point iteration with Newton method.
 * 
 * @param p0 
 * @param p1 
 * @param config 
 */
void pushForward_iteration(
    Particle2d* p0,
    Particle2d* p1,
    Config* config
) {
    VectorXd f = VectorXd::Zero(2*config->nmarkers);
    VectorXd f1 = VectorXd::Zero(2*config->nmarkers);
    VectorXd dv = VectorXd::Zero(2*config->nmarkers);
    MatrixXd Jf = MatrixXd::Zero(2*config->nmarkers, 2*config->nmarkers);

    f_eqmotion(&f, p0, p1, config);
    print_out(VERBOSE_DEBUG, "Eq. of motions precision: %e:\n", f.norm());

    for (int i=0;i<config->nmarkers; i++){
        for (int j=0; j<2; j++) {
            p1[i].z[j] += config->dx;

            f_eqmotion(&f1, p0, p1, config);

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
    Config* config
) {
    MatrixXd QGamma(2*config->nmarkers, config->nmarkers);
    buildQGamma(&QGamma, p0, p1, config);

    int idx;
    for (int i=0; i<config->nmarkers; i++) {
        for (int j=0; j<2; j++) {
            idx = 2*i+j;
            f->coeffRef(idx) = (p1[i].z[j] - p0[i].z[j]) / config->dt;

            for (int k=0; k<config->nmarkers; k++) {
                f->coeffRef(idx) -= config->nu / config->m * p1[k].weight * QGamma(idx, k);
            }
                            
        }
    }
}

/**
 * @brief Build the product of Q(v_p^{n+1/2} - v_p'^{n+1/2}) and /Gamma(S_eps^n, p, p')
 * 
 * @param ret 
 * @param p0 
 * @param p1 
 * @param config 
 */
void buildQGamma(
    MatrixXd* ret,
    Particle2d* p0,
    Particle2d* p1,
    Config* config
) {
    Vector2d gammaTmp;
    Matrix2d qTmp;
    VectorXd dSdv(2*config->nmarkers);
    computedSdv(&dSdv, p0, config);
    for (int i=0; i<config->nmarkers; i++) {
        for (int j=0; j<config->nmarkers; j++) {
            if (i==j) {
                // diagonal, set to 0 (Q*Gamma is antisysmmetric)
                (*ret).block(2*i, j, 2, 1).setZero();
            } else if (j<i) {
                (*ret).block(2*i, j, 2, 1) = (*ret).block(2*j, i, 2, 1);
            } else {
                if ((i==18) && (j==24)) {
                    int asd = 1;
                }
                Q(&qTmp, (p1[i].z + p0[i].z - p1[j].z - p0[j].z) / 2);
                gammaTmp = dSdv.segment(2*i, 2) - dSdv.segment(2*j, 2);
                (*ret).block(2*i, j, 2, 1) = qTmp * gammaTmp;
            }
        }
    }
}

/**
 * @brief Compute 1 / (m*w_p) * dS / dvp
 * 
 * @param ret 
 * @param p 
 * @param config 
 */
void computedSdv(
    VectorXd* ret,
    Particle2d* p,
    Config* config
) {
    ret->setZero();
    double tmp, k;
    for (int i_p1 = 0; i_p1<config->nmarkers; i_p1++)
    for (int i_x = 0; i_x < 2; i_x++)
    for (int i=0; i<config->nHermite; i++)
    for (int j=0; j<config->nHermite; j++) {
        tmp = 0;
        for (int i_p2 = 0; i_p2<config->nmarkers; i_p2++) {
            tmp += p[i_p2].weight * (
                        pow(config->kHermite[i] + (p[i_p1].z[0] - p[i_p2].z[0])/ sqrt(2*config->eps), 2)
                       + pow(config->kHermite[j] + (p[i_p1].z[1] - p[i_p2].z[1])/ sqrt(2*config->eps), 2)
                    );
        }
        if (i_x == 0) {
            k = config->kHermite[i];
        } else {
            k = config->kHermite[j];
        }
        ret->coeffRef(2*i_p1+i_x) += k * config->wHermite[i] * config->wHermite[j] * (1. - tmp / (2*CONST_PI*config->eps));
    }
    (*ret) /= -(config->m * CONST_PI * config->eps);
}


/**
 * @brief Perpendicular projection operator Q (eq. 10)
 * 
 * @param ret 
 * @param v 
 */
void Q(Matrix2d* ret, Vector2d v) {
    double norm2 = v.squaredNorm();
    if (norm2 < 1E-14) {
        ret->setZero();
    } else {
        ret->coeffRef(0,0) = - v(0)*v(0) / norm2 + 1;
        ret->coeffRef(1,1) = - v(1)*v(1) / norm2 + 1;
        ret->coeffRef(1,0) = - v(1)*v(0) / norm2;
        ret->coeffRef(0,1) = ret->coeffRef(1,0);
        *ret /= sqrt(norm2);
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