#ifndef CUDA

#include "coulomb_kernel.h"

namespace Kernel {


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

    if (VERBOSE_LEVEL >= VERBOSE_SILLY) {
        cout << "==== dSdV" << endl;
        cout << (*ret) << endl;
        cout << "==== dSdV end" << endl;
    }
}

}

#endif