#ifdef CUDA

#include "coulomb_kernel.h"

static void HandleError( cudaError_t err, const char *file, int line )
{
	// CUDA error handeling from the "CUDA by example" book
	if (err != cudaSuccess)
    {
		printf( "%s in %s at line %d\n", cudaGetErrorString( err ), file, line );
		exit( EXIT_FAILURE );
	}
}

#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

namespace Kernel {

    /**
     * @brief Perpendicular projection operator Q (eq. 10)
     * 
     * @param ret 
     * @param v 
     */
    __device__ void Q(Matrix2d* ret, Vector2d v) {
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
    __global__ void cuda_f_eq_motion(
        double* ret,
        Particle2d* p0,
        Particle2d* p1,
        double* dSdV,
        Config* config
    ) {
        Vector2d gammaTmp;
        Matrix2d qTmp;

        int i = blockIdx.x*config->cudaThreadsPerBlock + threadIdx.x;
        if (i < config->nmarkers) {
            ret[2*i] = 0;
            ret[2*i+1] = 0;
            for (int j=0; j<config->nmarkers; j++) {
                if (i!=j) {
                    Q(&qTmp, (p1[i].z + p0[i].z - p1[j].z - p0[j].z) / 2);
                    gammaTmp(0) = dSdV[2*i] - dSdV[2*j];
                    gammaTmp(1) = dSdV[2*i+1] - dSdV[2*j+1];
                    gammaTmp = qTmp*gammaTmp;
                    ret[2*i] += config->nu / config->m * p1[j].weight * gammaTmp(0);
                    ret[2*i+1] += config->nu / config->m * p1[j].weight * gammaTmp(1);
                }
            }
        }
    }

    __global__ void cuda_dSdv(
        double* ret,
        Particle2d* p,
        Config* config 
    ) {
        int i_p1 = blockIdx.x*config->cudaThreadsPerBlock + threadIdx.x;
        if (i_p1 < config->nmarkers) {
            double logsum, k;
            for (int i_x = 0; i_x < 2; i_x++) {
                ret[2*i_p1+i_x] = 0;
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
                    ret[2*i_p1+i_x] += k * config->wHermite[i] * config->wHermite[j] * (1. + logsum);
                }
                ret[2*i_p1+i_x] *= sqrt(2.*config->eps) / (config->m * CONST_PI * config->eps);
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
        // CUDA blocks configuration
        int nblocks = ceil(float(config->nmarkers) / config->cudaThreadsPerBlock);

        // // Allocate device arrays
        double *d_ret, *d_dSdv;
        Config* d_config;
        Particle2d *d_p0, *d_p1;
        HANDLE_ERROR(cudaMalloc((void **)&d_ret, sizeof(double)*2*config->nmarkers));
        HANDLE_ERROR(cudaMalloc((void **)&d_dSdv, sizeof(double)*2*config->nmarkers));
        HANDLE_ERROR(cudaMalloc((void **)&d_config, sizeof(Config)));
        HANDLE_ERROR(cudaMalloc((void **)&d_p0, sizeof(Particle2d)*config->nmarkers));
        HANDLE_ERROR(cudaMalloc((void **)&d_p1, sizeof(Particle2d)*config->nmarkers));

        // Copy to device
        HANDLE_ERROR(cudaMemcpy(d_config, config, sizeof(Config), cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(d_dSdv, dSdV->data(), sizeof(double)*dSdV->size(), cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(d_p0, p0, sizeof(Particle2d)*config->nmarkers, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(d_p1, p1, sizeof(Particle2d)*config->nmarkers, cudaMemcpyHostToDevice));

        // Compute the entropy gradient
        // cudaThreadSynchronize();
        cuda_f_eq_motion<<<nblocks, config->cudaThreadsPerBlock>>>(d_ret, d_p0, d_p1, d_dSdv, d_config);
        HANDLE_ERROR( cudaPeekAtLastError() );
        // cudaThreadSynchronize();
        
        // Copy to host
        HANDLE_ERROR(cudaMemcpy(dv->data(), d_ret, sizeof(double)*dv->size(), cudaMemcpyDeviceToHost));
    }


    void computedSdv(
        VectorXd* ret,
        Particle2d* p,
        Config* config
    ) {

        // CUDA blocks configuration
        int nblocks = ceil(float(config->nmarkers) / config->cudaThreadsPerBlock);

        // // Allocate device arrays
        double *d_ret;
        Config* d_config, *l_config = new Config;
        Particle2d* d_p;
        double* d_kHermite, *d_wHermite;
        Vector2d *d_u1, *d_u2;
        HANDLE_ERROR(cudaMalloc((void **)&d_ret, sizeof(double)*2*config->nmarkers));
        HANDLE_ERROR(cudaMalloc((void **)&d_config, sizeof(Config)));
        HANDLE_ERROR(cudaMalloc((void **)&d_p, sizeof(Particle2d)*config->nmarkers));
        HANDLE_ERROR(cudaMalloc((void **)&d_kHermite, sizeof(double)*config->nHermite));
        HANDLE_ERROR(cudaMalloc((void **)&d_wHermite, sizeof(double)*config->nHermite));
        HANDLE_ERROR(cudaMalloc((void **)&d_u1, sizeof(Vector2d)));
        HANDLE_ERROR(cudaMalloc((void **)&d_u2, sizeof(Vector2d)));

        // // Copy to device
        memcpy(l_config, config, sizeof(Config));
        HANDLE_ERROR(cudaMemcpy(d_kHermite, config->kHermite, sizeof(double)*config->nHermite, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(d_wHermite, config->wHermite, sizeof(double)*config->nHermite, cudaMemcpyHostToDevice));
        l_config->kHermite = d_kHermite;
        l_config->wHermite = d_wHermite;
        HANDLE_ERROR(cudaMemcpy(d_config, l_config, sizeof(Config), cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(d_p, p, sizeof(Particle2d)*config->nmarkers, cudaMemcpyHostToDevice));

        // // Compute the entropy gradient
        // cudaThreadSynchronize();
        // printf("computing entropy\n");
        cuda_dSdv<<<nblocks, config->cudaThreadsPerBlock>>>(d_ret, d_p, d_config);
        HANDLE_ERROR( cudaPeekAtLastError() );
        // cudaThreadSynchronize();
        
        // Copy to host
        HANDLE_ERROR(cudaMemcpy(ret->data(), d_ret, sizeof(double)*ret->size(), cudaMemcpyDeviceToHost));
    }
}

#endif