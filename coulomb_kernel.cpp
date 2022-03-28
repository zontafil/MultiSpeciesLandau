#ifndef CUDA

#include "coulomb_kernel.h"

namespace Kernel {

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
    VectorXd* dSdV,
    Config* config
) {
    Vector2d gammaTmp;
    Matrix2d qTmp;

    for (int i=0; i<config->nmarkers; i++) {
        for (int j=0; j<config->nmarkers; j++) {
            if (i==j) {
                // diagonal, set to 0 (Q*Gamma is antisysmmetric)
                (*ret).block(2*i, j, 2, 1).setZero();
            } else if (j<i) {
                (*ret).block(2*i, j, 2, 1) = -(*ret).block(2*j, i, 2, 1);
            } else {
                Q(&qTmp, (p1[i].z + p0[i].z - p1[j].z - p0[j].z) / 2);
                gammaTmp = dSdV->segment(2*i, 2) - dSdV->segment(2*j, 2);
                (*ret).block(2*i, j, 2, 1) = qTmp * gammaTmp; // WHY - SIGN ?
            }
        }
    }
}

/**
 * @brief Compute dv of the equations of motion
 * 
 * @param dv 
 * @param p0 
 * @param p1 
 * @param dSdV 
 * @param config 
 */
void f_eqmotion_dv(
    VectorXd* dv,
    Particle2d* p0,
    Particle2d* p1,
    VectorXd* dSdV,
    Config* config
) {
    MatrixXd QGamma(2*config->nmarkers, config->nmarkers);
    buildQGamma(&QGamma, p0, p1, dSdV, config);

    int idx;
    for (int i=0; i<config->nmarkers; i++) {
        for (int j=0; j<2; j++) {
            idx = 2*i+j;

            dv->coeffRef(idx) = 0;
            for (int k=0; k<config->nmarkers; k++) {
                dv->coeffRef(idx) -= config->nu / config->m * p1[k].weight * QGamma(idx, k);
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
    double logsum, k;
    for (int i_p1 = 0; i_p1<config->nmarkers; i_p1++)
    for (int i_x = 0; i_x < 2; i_x++)
    for (int i=0; i<config->nHermite; i++)
    for (int j=0; j<config->nHermite; j++) {
        logsum = 0;
        for (int i_p2 = 0; i_p2<config->nmarkers; i_p2++) {
            logsum += p[i_p2].weight / (CONST_2PI*config->eps)* ( exp(
                        -pow(config->kHermite[i] + (p[i_p1].z[0] - p[i_p2].z[0])/ sqrt(2*config->eps), 2)
                        -pow(config->kHermite[j] + (p[i_p1].z[1] - p[i_p2].z[1])/ sqrt(2*config->eps), 2))
                    );
        }
        logsum = log(logsum);
        if (i_x == 0) {
            k = config->kHermite[i];
        } else {
            k = config->kHermite[j];
        }
        ret->coeffRef(2*i_p1+i_x) += k * config->wHermite[i] * config->wHermite[j] * (1. + logsum);
    }
    (*ret) *= sqrt(2.*config->eps) / (config->m * CONST_PI * config->eps);
}

}

#endif