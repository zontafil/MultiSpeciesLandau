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
        cudaThreadSynchronize();
        printf("computing entropy\n");
        cuda_dSdv<<<nblocks, config->cudaThreadsPerBlock>>>(d_ret, d_p, d_config);
        HANDLE_ERROR( cudaPeekAtLastError() );
        cudaThreadSynchronize();
        
        // Copy to host
        HANDLE_ERROR(cudaMemcpy(ret->data(), d_ret, sizeof(double)*ret->size(), cudaMemcpyDeviceToHost));
    }
}

#endif