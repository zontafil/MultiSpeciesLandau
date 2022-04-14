#include "coulombStructurePreserving.h"

using std::setprecision;

namespace Coulomb {

/**
 * @brief Run the simulation and write output to file (step_C*.txt)
 * 
 * @param config 
 */
void Run(Config* config) {
    Particle2d** p_mesh = new Particle2d*[config->nspecies];
    Particle2d** p0 = new Particle2d*[config->nspecies];
    Particle2d** p1 = new Particle2d*[config->nspecies];
    for (int i=0; i<config->nspecies; i++) {
        p_mesh[i] = initMarkers(i, config, MESH);
        p0[i] = initMarkers(i, config, config->distributionType);
        p1[i] = initMarkers(i, config, config->distributionType);
    }

    // initial energy
    double E0 = K(p1, config), E;
    double** f_mesh = new double*[config->nspecies];
    Vector2d P0 = Momentum(p1, config), P;
    VectorXd* dSdV = new VectorXd[config->nspecies];
    for (int s=0; s<config->nspecies; s++) {
        f_mesh[s] = new double[config->nmarkers];
        dSdV[s] = VectorXd(2*config->nmarkers);
    }

    print_out(VERBOSE_NORMAL, "Markers: %d dT: %f\n", config->nmarkers, config->dt);
    print_out(VERBOSE_NORMAL, "Initial Energy: %e\n", E0);
    print_out(VERBOSE_NORMAL, "Initial Momentum: %e %e\n", P0[0], P0[1]);

    for (int t=0; t<config->n_timesteps; t++) {

        print_out(VERBOSE_NORMAL, "Timestep: %d\n", t);

        for (int s=0; s<config->nspecies; s++) {
            copy(p1[s], p1[s]+config->nmarkers, p0[s]);
        }

        if (t%config->recordAtStep == 0) {
            printState(f_mesh, p_mesh, p1, config, t, E0, P0);
        }

        // precompute entropy gradient
        Kernel::computedSdv(dSdV, p0, config);
        if (VERBOSE_LEVEL >= VERBOSE_SILLY) {
            cout << "==== dSdV" << endl;
            for (int s=0; s<config->nspecies; s++) {
                for (int i=0; i<config->nmarkers; i++) {
                    printf("%e\t%e\n", dSdV[s](2*i), dSdV[s](2*i+1));
                }
            }
            cout << "==== dSdV end" << endl;
        }

        // fixed point newton iterations
        for (int j=0; j<config->maxEOMIterations; j++) {
            print_out(VERBOSE_DEBUG, "Iteration %d ", j);
            if (config->useNewton) {
                if (pushForwardNewtonIteration(p0, p1, dSdV, config)) {
                    break;
                }
            } else {
                if (pushForward_dv(p0, p1, dSdV, config)) {
                    break;
                }
            }
        }
    }

    printState(f_mesh, p_mesh, p1, config, config->n_timesteps, E0, P0);

    for (int s=0; s<config->nspecies; s++) {
        free(p0[s]);
        free(p1[s]);
        free(p_mesh[s]);
        free(f_mesh[s]);
    }
}

/**
 * @brief Print the current state to screen and to file
 * 
 * @param f_mesh 
 * @param p_mesh 
 * @param p1 
 * @param config 
 * @param t 
 * @param E0 
 * @param P0 
 */
void printState(
    double** f_mesh,
    Particle2d** p_mesh,
    Particle2d** p1,
    Config* config,
    int t,
    double E0,
    Vector2d P0
) {
    // build distribution at mesh nodes and print to file
    mesh_distribution(f_mesh, p_mesh, p1, config);
    char filename[30];
    snprintf (filename, sizeof filename, "out/data/step_C_%05d.txt", t);
    FILE* fout = fopen(filename, "w+");

    // print system state and debug info to screen
    double E = K(p1, config);
    print_out(VERBOSE_NORMAL, "Energy: %.15e Error: %.15e\n", E, (E-E0)/E0);

    Vector2d P = Momentum(p1, config);
    print_out(VERBOSE_NORMAL, "Momentum: %.15e %.15e\n", P[0], P[1]);
    print_out(VERBOSE_NORMAL, "Momentum Error: %.15e %.15e\n", (P[0]-P0[0])/P0[0], (P[1]-P0[1])/P0[1]);

    double distmin = 1000000;
    for (int s1=0; s1<config->nspecies; s1++)
    for (int i1=0; i1<config->nmarkers; i1++)
    for (int s2=0; s2<config->nspecies; s2++)
    for (int i2=0; i2<config->nmarkers; i2++) {
        if (i1!=i2 || s1!=s2) {
            distmin = min(distmin, (p1[s1][i1].z - p1[s2][i2].z).norm());
        }
    }

    fprintf(fout, "%d %d\n", config->nspecies, config->nmarkers);
    fprintf(fout, "%e %e %e %e %e %e %e\n", E, (E-E0)/E0, P[0], P[1], (P[0]-P0[0])/P0[0], (P[1]-P0[1])/P0[1], distmin);
    for (int s=0; s<config->nspecies; s++) {
        E = Kspecie(p1, s, config);
        P = MomentumSpecie(p1, s, config);
        fprintf(fout, "%e %e %e %s\n", E, P[0], P[1], config->species[s].name);
    }
    for (int s=0; s<config->nspecies; s++) {
        for (int i=0; i<config->nmarkers; i++) {
            fprintf(fout, "%d %d %e %e %e\n", s, i, p_mesh[s][i].z[0], p_mesh[s][i].z[1], f_mesh[s][i]);
        }
    }


    if (VERBOSE_LEVEL >= VERBOSE_SILLY) {
        for (int s=0; s<config->nspecies; s++) {
            for (int i=0; i<config->nmarkers; i++) {
                print_out(VERBOSE_SILLY, "Specie %d ID: %d Vx: %.15e Vy: %.15e W: %.15e\n", s, i, p1[s][i].z[0], p1[s][i].z[1], p1[s][i].weight);
            }
        }
    }
}

/**
 * @brief Distribution function of initial states
 * 
 * @param v 
 * @param specie index of specie
 * @param config
 * @return double 
 */
double f(Vector2d v, int specie, Config* config) {
    double ret = 0;
    for (int i=0; i<config->species[specie].npeaks; i++) {
        ret += exp(-(v-config->species[specie].peaks[i]).squaredNorm()/2.);
    }

    ret *= pow(config->h,2)/(4.*CONST_PI);
    return ret;
}

/**
 * @brief Build distribution of states at the mesh nodes
 * 
 * @param ret 
 * @param p_mesh 
 * @param p 
 * @param config 
 */
void mesh_distribution(
    double** ret,
    Particle2d** p_mesh,
    Particle2d** p,
    Config* config
) {
    for (int s=0; s<config->nspecies; s++) {
        for (int i=0; i<config->nmarkers; i++) {
            ret[s][i] = 0;
            double CONST_2EPS_M1 = 1./(2.*config->eps);
            double CONST_2PIEPS_M1 = 1./(CONST_2PI*config->eps);
            for (int j=0; j<config->nmarkers; j++) {
                ret[s][i] += exp(-(p_mesh[s][i].z - p[s][j].z).squaredNorm()*CONST_2EPS_M1)*p[s][j].weight;
            }
            ret[s][i] *= CONST_2PIEPS_M1;
        }
    }
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
Particle2d* initMarkers(int specie, Config* config, DistributionType type) {
    int nmarkers = config->nmarkers;
    Particle2d* ret;
    ret = new Particle2d[config->nmarkers];
    int idx;
    for (int i = 0; i<config->ny; i++) {
        for (int j = 0; j<config->nx; j++) {
            idx = i*config->nx + j;
            if (type == UNIFORM) {
                ret[idx].z[0] = double(rand()) / RAND_MAX * (config->xmax - config->xmin) + config->xmin;
                ret[idx].z[1] = double(rand()) / RAND_MAX * (config->ymax - config->ymin) + config->ymin;
            } else if (type == MESH) {
                ret[idx].z[0] = double(j+0.5) / (config->nx) * (config->xmax-config->xmin) + config->xmin;
                ret[idx].z[1] = double(i+0.5) / (config->nx) * (config->ymax-config->ymin) + config->ymin;
            } else if (type == MESH_SHIFT) {
                double shift = 0;
                if (specie == 1) {
                    shift = 0.5;
                }
                ret[idx].z[0] = double(j+shift) / (config->nx) * (config->xmax-config->xmin) + config->xmin;
                ret[idx].z[1] = double(i+shift) / (config->nx) * (config->ymax-config->ymin) + config->ymin;
            } else if (type == MESH_PEAK_CENTERED) {
                double dx = 0.1;
                Vector2d peak = config->species[specie].peaks[0];
                ret[idx].z[0] = double(j+0.5) / (config->nx) * (config->xmax-config->xmin) + peak(0) - dx/2;
                ret[idx].z[1] = double(i+0.5) / (config->ny) * (config->ymax-config->ymin) + peak(1) - dx/2;
            }
            ret[idx].weight = f(ret[idx].z, specie, config);
        }
    }
    return ret;
}

int pushForward_dv(
    Particle2d** p0,
    Particle2d** p1,
    VectorXd* dSdV,
    Config* config
) {
    VectorXd* f = new VectorXd[config->nspecies];
    for (int s=0; s<config->nspecies; s++) {
        f[s] = VectorXd(2*config->nmarkers);
    }
    Kernel::f_eqmotion_dv(f, p0, p1, dSdV, config);
    if (VERBOSE_LEVEL >= VERBOSE_SILLY) {
        printf("=== dv\n");
        for (int s=0; s<config->nspecies; s++) {
            for (int i=0; i<config->nmarkers; i++) {
                printf("%d %e\t%e\n", s, f[s](2*i), f[s](2*i+1));
            }
        }
        printf("=== dv end\n");
    }

    double znew, err = 0, dz;
    for (int s=0; s<config->nspecies; s++) {
        for (int i=0; i<config->nmarkers; i++) {
            for (int j=0; j<2; j++) {
                znew = p0[s][i].z[j] + config->dt * f[s](i*2+j);
                dz = p1[s][i].z[j] - znew;
                err += abs(dz);
                p1[s][i].z[j] = znew;
            }
        }
    }
    err /= (config->nmarkers * config->nspecies);

    print_out(VERBOSE_DEBUG, "Eq. of motions precision: %e\n", err);

    if (err < config->newtonTolerance) {
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
int pushForwardNewtonIteration(
    Particle2d** p0,
    Particle2d** p1,
    VectorXd* dSdV,
    Config* config
) {
    VectorXd* f = new VectorXd[config->nspecies];
    VectorXd* f1 = new VectorXd[config->nspecies];
    VectorXd dv(2*config->nmarkers);
    MatrixXd Jf(2*config->nmarkers, 2*config->nmarkers);
    for (int s=0; s<config->nspecies; s++) {
        f[s] = VectorXd(2*config->nmarkers);
        f1[s] = VectorXd(2*config->nmarkers);
    }

    f_eqmotion(f, p0, p1, dSdV, config);
    double fnorm = 0, fmax = 0;
    for (int s=0; s<config->nspecies; s++) {
        fnorm += f[s].squaredNorm();
        fmax = max(f->maxCoeff(), fmax);
    }
    fnorm = sqrt(fnorm) / config->nspecies / config->nmarkers;
    print_out(VERBOSE_DEBUG, "Eq. of motions precision: %e max: %e\n", fnorm, fmax);
    if (fnorm < config->newtonTolerance) {
        return 1;
    }

    for (int s=0; s<config->nspecies; s++) {
        for (int i=0;i<config->nmarkers; i++){
            for (int j=0; j<2; j++) {
                p1[s][i].z[j] += config->dx;

                f_eqmotion(f1, p0, p1, dSdV, config);

                Jf.col(i*2+j) = (f1[s] - f[s])/config->dx;

                p1[s][i].z[j] -= config->dx;
            }
        }
        dv = -Jf.inverse()*f[s];

        // copy back from dv to the markers array
        for (int i=0;i<config->nmarkers; i++){
            for (int j=0; j<2; j++) {
                p1[s][i].z[j] += dv[i*2+j];
            }
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
    Particle2d** p0,
    Particle2d** p1,
    VectorXd* dSdV,
    Config* config
) {

    Kernel::f_eqmotion_dv(f, p0, p1, dSdV, config);

    for (int s=0; s<config->nspecies; s++) {
        for (int i=0; i<config->nmarkers; i++) {
            for (int j=0; j<2; j++) {
                f[s](i*2+j) = (p1[s][i].z[j] - p0[s][i].z[j]) / config->dt - f[s](i*2+j);
            }
        }
    }

    if (VERBOSE_LEVEL >= VERBOSE_SILLY) {
        printf("=== fcoeff\n");
        for (int s=0; s<config->nspecies; s++) {
            for (int i=0; i<config->nmarkers; i++) {
                printf("%d %e %e\n", s, f[s](2*i), f[s](2*i+1));
            }
        }
        printf("=== fcoeff end\n");
    }
}

/**
 * @brief Compute the kinetic energy of a specie
 * 
 * @param p 
 * @param s index of the specie
 * @param config 
 * @return double 
 */
double Kspecie(Particle2d** p, int s, Config* config) {
    double ret = 0;
    for (int i = 0; i<config->nmarkers; i++) {
        ret += p[s][i].weight * config->species[s].m * 0.5 * p[s][i].z.squaredNorm();
    }
    return ret;
}

/**
 * @brief Compute the kinetic energy of the system
 * 
 * @param p 
 * @param config 
 * @return double 
 */
double K(Particle2d** p, Config* config) {
    double ret = 0;
    for (int s=0; s<config->nspecies; s++) {
        ret += Kspecie(p, s, config);
    }
    return ret;
}

/**
 * @brief Compute momentum of a specie
 * 
 * @param p 
 * @param s 
 * @param config 
 * @return Vector2d 
 */
Vector2d MomentumSpecie(Particle2d** p, int s, Config* config) {
    Vector2d ret(0,0);
    for (int i=0; i<config->nmarkers; i++) {
        ret += p[s][i].weight * config->species[s].m * p[s][i].z;
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
Vector2d Momentum(Particle2d** p, Config* config) {
    Vector2d ret(0,0);
    for (int s=0; s<config->nspecies; s++) {
        ret += MomentumSpecie(p, s, config);
    }
    return ret;
}

}