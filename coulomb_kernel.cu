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
        double** ret,
        Particle2d** p0,
        Particle2d** p1,
        double** dSdV,
        Config* config,
        Specie* species
    ) {
        Vector2d gammaTmp, qGamma, z0_is, z1_is, dSdV_is;
        Matrix2d qTmp;
        Particle2d *p0_s2, *p1_s2;
        double nu, m;
        double *dSdV_s2;
        double *ret_s;

        // this thread is computing the EOM contribution of marker i in specie s
        // we perform a sum in all the markers j of all the species s2
        int global_idx = blockIdx.x*config->cudaThreadsPerBlock + threadIdx.x;
        int s = global_idx / config->nmarkers;
        int i = global_idx % config->nmarkers;
        if (s < config->nspecies) {

            // set some pointer to s block, to avoid doing additional computations
            ret[s][2*i] = 0;
            ret[s][2*i+1] = 0;
            ret_s = ret[s];
            z0_is = p0[s][i].z;
            z1_is = p1[s][i].z;
            dSdV_is(0) = dSdV[s][2*i];
            dSdV_is(1) = dSdV[s][2*i+1];
            m = species[s].m;
            for (int s2=0; s2<config->nspecies; s2++) {
                // set some pointer to s2 block, to avoid doing additional computations
                p0_s2 = p0[s2];
                p1_s2 = p1[s2];
                nu = species[s].nu[s2];
                dSdV_s2 = dSdV[s2];
                for (int j=0; j<config->nmarkers; j++) {
                    if (i!=j || s!=s2) {
                        Q(&qTmp, (z1_is + z0_is - p1_s2[j].z - p0_s2[j].z) / 2);
                        gammaTmp(0) = dSdV_is(0) - dSdV_s2[2*j];
                        gammaTmp(1) = dSdV_is(1) - dSdV_s2[2*j+1];
                        gammaTmp = qTmp*gammaTmp;
                        ret_s[2*i] -= nu / m * p1_s2[j].weight * gammaTmp(0);
                        ret_s[2*i+1] -= nu / m * p1_s2[j].weight * gammaTmp(1);
                    }
                }
            }
        }
    }

    /**
     * @brief CUDA version of entropy gradient
     */
    __global__ void cuda_dSdv(
        double** ret,
        Particle2d** p,
        Config* config,
        Specie* species
    ) {
        int global_idx = blockIdx.x*config->cudaThreadsPerBlock + threadIdx.x;
        int s = global_idx / config->nmarkers;
        int s2 = s;
        int i_p1 = global_idx % config->nmarkers;
        if (s < config->nspecies) {
            Particle2d* ps1 = p[s];
            double* rets1 = ret[s];
            double m = species[s].m;
            double logsum, dx, dy, kpx, kpy;
            double SQRT2EPSM1 = 1./sqrt(2.*config->eps);
            double PI2EPSM1 = 1./(CONST_2PI * config->eps);
            rets1[2*i_p1] = 0;
            rets1[2*i_p1+1] = 0;
            for (int i=0; i<config->nHermite; i++)
            for (int j=0; j<config->nHermite; j++) {
                logsum = 0;

                // TODO: Normalize z to SQRT2EPSM1 --> ~10% performance boost
                kpx = config->kHermite[i] + ps1[i_p1].z[0] * SQRT2EPSM1;
                kpy = config->kHermite[j] + ps1[i_p1].z[1] * SQRT2EPSM1;

                for (int i_p2 = 0; i_p2<config->nmarkers; i_p2++) {
                    dx = kpx - p[s2][i_p2].z[0] * SQRT2EPSM1;
                    dy = kpy - p[s2][i_p2].z[1] * SQRT2EPSM1;
                    logsum+=p[s2][i_p2].weight* exp(-dx*dx - dy*dy);
                }
                logsum = config->wHermite[i]*config->wHermite[j] * (1. + log(logsum * PI2EPSM1));
                rets1[2*i_p1] += logsum * config->kHermite[i];
                rets1[2*i_p1+1] += logsum * config->kHermite[j];
            }
            rets1[2*i_p1] *= sqrt(2.*config->eps) / (m * CONST_PI * config->eps);
            rets1[2*i_p1+1] *= sqrt(2.*config->eps) / (m * CONST_PI * config->eps);
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
        Particle2d** p0,
        Particle2d** p1,
        VectorXd* dSdV,
        Config* config
    ) {
        // CUDA blocks configuration
        int nblocks = ceil(float(config->nmarkers * config->nspecies) / config->cudaThreadsPerBlock);

        // init particles
        Particle2d** h_p0 = new Particle2d*[config->nspecies];
        Particle2d** h_p1 = new Particle2d*[config->nspecies];
        Particle2d** d_p0, **d_p1;
        HANDLE_ERROR(cudaMalloc((void **)&d_p0, sizeof(Particle2d*)*config->nspecies));    
        HANDLE_ERROR(cudaMalloc((void **)&d_p1, sizeof(Particle2d*)*config->nspecies));    
        for (int s=0; s<config->nspecies; s++) {
            HANDLE_ERROR(cudaMalloc((void **)&(h_p0[s]), sizeof(Particle2d)*config->nmarkers));
            HANDLE_ERROR(cudaMalloc((void **)&(h_p1[s]), sizeof(Particle2d)*config->nmarkers));
            HANDLE_ERROR(cudaMemcpy(h_p0[s], p0[s], sizeof(Particle2d)*config->nmarkers, cudaMemcpyHostToDevice));
            HANDLE_ERROR(cudaMemcpy(h_p1[s], p1[s], sizeof(Particle2d)*config->nmarkers, cudaMemcpyHostToDevice));
        }
        HANDLE_ERROR(cudaMemcpy (d_p0, h_p0, config->nspecies*sizeof(Particle2d*), cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy (d_p1, h_p1, config->nspecies*sizeof(Particle2d*), cudaMemcpyHostToDevice));

        // init config
        Config* d_config;
        HANDLE_ERROR(cudaMalloc((void **)&d_config, sizeof(Config)));
        HANDLE_ERROR(cudaMemcpy(d_config, config, sizeof(Config), cudaMemcpyHostToDevice));

        // init spieces config
        Specie h_species[config->nspecies], *d_species;
        HANDLE_ERROR(cudaMalloc((void **)&d_species, sizeof(Specie)*config->nspecies));    
        for (int s=0; s<config->nspecies; s++) {
            memcpy(&h_species[s], &config->species[s], sizeof(Specie));
            HANDLE_ERROR(cudaMalloc((void **)&(h_species[s].nu), sizeof(double)*config->nspecies));    
            HANDLE_ERROR(cudaMemcpy(h_species[s].nu, config->species[s].nu, sizeof(double)*config->nspecies, cudaMemcpyHostToDevice));
        }
        HANDLE_ERROR(cudaMemcpy (d_species, h_species, config->nspecies*sizeof(Specie), cudaMemcpyHostToDevice));

        // init entropy
        double** h_dSdv = new double*[config->nspecies], **d_dSdv;
        HANDLE_ERROR(cudaMalloc((void **)&d_dSdv, sizeof(double*)*config->nspecies));    
        for (int s=0; s<config->nspecies; s++) {
            HANDLE_ERROR(cudaMalloc((void **)&h_dSdv[s], sizeof(double)*2*config->nmarkers));
            HANDLE_ERROR(cudaMemcpy(h_dSdv[s], dSdV[s].data(), sizeof(double)*dSdV[s].size(), cudaMemcpyHostToDevice));
        }
        HANDLE_ERROR(cudaMemcpy (d_dSdv, h_dSdv, config->nspecies*sizeof(double*), cudaMemcpyHostToDevice));

        // init ret
        double** h_ret = new double*[config->nspecies], **d_ret;
        HANDLE_ERROR(cudaMalloc((void **)&d_ret, sizeof(double*)*config->nspecies));    
        for (int s=0; s<config->nspecies; s++) {
            HANDLE_ERROR(cudaMalloc((void **)&h_ret[s], sizeof(double)*2*config->nmarkers));
        }
        HANDLE_ERROR(cudaMemcpy (d_ret, h_ret, config->nspecies*sizeof(double*), cudaMemcpyHostToDevice));

        // Compute the entropy gradient
        cuda_f_eq_motion<<<nblocks, config->cudaThreadsPerBlock>>>(d_ret, d_p0, d_p1, d_dSdv, d_config, d_species);
        HANDLE_ERROR( cudaPeekAtLastError() );
        
        // Copy to host
        HANDLE_ERROR(cudaMemcpy(h_ret, d_ret, sizeof(double)*config->nspecies, cudaMemcpyDeviceToHost));
        for (int s=0; s<config->nspecies; s++) {
            HANDLE_ERROR(cudaMemcpy(dv[s].data(), h_ret[s], sizeof(double)*dv[s].size(), cudaMemcpyDeviceToHost));
        }
    }


    void computedSdv(
        VectorXd* ret,
        Particle2d** p,
        Config* config
    ) {

        // CUDA blocks configuration
        int nblocks = ceil(float(config->nmarkers) / config->cudaThreadsPerBlock);

        // // Allocate device arrays
        Config* d_config, *l_config = new Config;
        Particle2d** d_p;
        Particle2d** h_p = new Particle2d*[config->nspecies];

        double* d_kHermite, *d_wHermite;
        HANDLE_ERROR(cudaMalloc((void **)&d_p, sizeof(Particle2d*)*config->nspecies));    
        HANDLE_ERROR(cudaMalloc((void **)&d_config, sizeof(Config)));
        HANDLE_ERROR(cudaMalloc((void **)&d_kHermite, sizeof(double)*config->nHermite));
        HANDLE_ERROR(cudaMalloc((void **)&d_wHermite, sizeof(double)*config->nHermite));
        for (int s=0; s<config->nspecies; s++) {
            HANDLE_ERROR(cudaMalloc((void **)&(h_p[s]), sizeof(Particle2d)*config->nmarkers));
            HANDLE_ERROR(cudaMemcpy(h_p[s], p[s], sizeof(Particle2d)*config->nmarkers, cudaMemcpyHostToDevice));
        }
        HANDLE_ERROR(cudaMemcpy (d_p, h_p, config->nspecies*sizeof(Particle2d*), cudaMemcpyHostToDevice));

        // init ret
        double** h_ret = new double*[config->nspecies], **d_ret;
        HANDLE_ERROR(cudaMalloc((void **)&d_ret, sizeof(double*)*config->nspecies));    
        for (int s=0; s<config->nspecies; s++) {
            HANDLE_ERROR(cudaMalloc((void **)&h_ret[s], sizeof(double)*2*config->nmarkers));
        }
        HANDLE_ERROR(cudaMemcpy (d_ret, h_ret, config->nspecies*sizeof(double*), cudaMemcpyHostToDevice));

        // init spieces config
        Specie h_species[config->nspecies], *d_species;
        HANDLE_ERROR(cudaMalloc((void **)&d_species, sizeof(Specie)*config->nspecies));    
        for (int s=0; s<config->nspecies; s++) {
            memcpy(&h_species[s], &config->species[s], sizeof(Specie));
            HANDLE_ERROR(cudaMalloc((void **)&(h_species[s].nu), sizeof(double)*config->nspecies));    
            HANDLE_ERROR(cudaMemcpy(h_species[s].nu, config->species[s].nu, sizeof(double)*config->nspecies, cudaMemcpyHostToDevice));
        }
        HANDLE_ERROR(cudaMemcpy (d_species, h_species, config->nspecies*sizeof(Specie), cudaMemcpyHostToDevice));

        // // Copy to device
        memcpy(l_config, config, sizeof(Config));
        HANDLE_ERROR(cudaMemcpy(d_kHermite, config->kHermite, sizeof(double)*config->nHermite, cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(d_wHermite, config->wHermite, sizeof(double)*config->nHermite, cudaMemcpyHostToDevice));
        l_config->kHermite = d_kHermite;
        l_config->wHermite = d_wHermite;
        HANDLE_ERROR(cudaMemcpy(d_config, l_config, sizeof(Config), cudaMemcpyHostToDevice));

        // for (int s=0; s<config->nspecies; s++) {
        // // Compute the entropy gradient
        // HANDLE_ERROR(cudaMemcpy(d_p, p[s], sizeof(Particle2d)*config->nmarkers, cudaMemcpyHostToDevice));
        // HANDLE_ERROR(cudaMemcpy(d_p2, p[1], sizeof(Particle2d)*config->nmarkers, cudaMemcpyHostToDevice));
        cuda_dSdv<<<nblocks, config->cudaThreadsPerBlock>>>(d_ret, d_p, d_config, d_species);
        HANDLE_ERROR( cudaPeekAtLastError() );

        // Copy to host
        for (int s=0; s<config->nspecies; s++) {
            HANDLE_ERROR(cudaMemcpy(ret[s].data(), h_ret[s], sizeof(double)*ret[s].size(), cudaMemcpyDeviceToHost));
        }
        // HANDLE_ERROR(cudaMemcpy(ret[s].data(), d_ret, sizeof(double)*ret[s].size(), cudaMemcpyDeviceToHost));
        // }
    }
}

#endif