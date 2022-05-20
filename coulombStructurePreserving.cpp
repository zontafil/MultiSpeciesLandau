#include "coulombStructurePreserving.h"

using std::setprecision;

namespace Coulomb {

/**
 * @brief Run the simulation and write output to file (step_C*.txt)
 * 
 * @param config 
 */
void Run(Config* config0) {
    Config* config;
    if (config0->normalize) {
        config = normalizeConfig(config0);
    } else {
        config = config0;
    }
    Particle2d** p_mesh = new Particle2d*[config->nspecies];
    Particle2d** p0 = new Particle2d*[config->nspecies];
    Particle2d** p1 = new Particle2d*[config->nspecies];
    for (int i=0; i<config->nspecies; i++) {
        Specie specie = config->species[i];
        p_mesh[i] = initMarkers(i, config, MESH);
        p0[i] = initMarkers(i, config, config->distributionType);
        p1[i] = initMarkers(i, config, config->distributionType);
        printf("Specie %d, n [1]: %e\n", i, nSpecie(p1, i, config));
        printf("Specie %d, T0x [eV]: %e\n", i, specie.Tx);
        printf("Specie %d, T0y [eV]: %e\n", i, specie.Ty);
        printf("Specie %d, T0comp [eV]: %e\n", i, TemperatureSpecie(p0, i, config));
        printf("Specie %d, m [kg]: %e\n", i, specie.m);
        printf("Specie %d, vmax [ms^-1] %e vmin [ms^-1] %e\n", i, specie.xmax, specie.xmin);
        for (int s=0; s<config->nspecies; s++) {
            printf("nu%d%d %e\n", i,s, specie.nu[s]);
        }
    }

    

    // initial energy
    double E0 = K(p1, config), E;
    double** f_mesh = new double*[config->nspecies];
    Vector2d P0 = Momentum(p1, config), P;
    VectorXd* dSdV = new VectorXd[config->nspecies];
    for (int s=0; s<config->nspecies; s++) {
        f_mesh[s] = new double[config->nmarkers];
        dSdV[s] = VectorXd(2*config->nmarkers);
        print_out(VERBOSE_NORMAL, "Specie %d Epsilon: %e\n", s, config->species[s].eps);
    }

    print_out(VERBOSE_NORMAL, "Markers: %d dT: %e\n", config->nmarkers, config->dt);
    print_out(VERBOSE_NORMAL, "Initial Energy [eV]: %e\n", E0);
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
 * @brief Normalize config object parameters
 * 
 * @param config 
 * @return Config* 
 */
Config* normalizeConfig(Config* config) {
    Config* ret = copyConfig(config);

    double v0 = 0;
    double n0 = 0;
    for (int s=0; s<ret->nspecies; s++) {
        v0 = fmax(v0, ret->species[s].xmax);
        v0 = fmax(v0, ret->species[s].ymax);
        n0 = fmax(n0, ret->species[s].n);
    }
    v0 /= 10.; // make the box = [-10,10]
    double nu0 = CONST_E * CONST_E * CONST_E * CONST_E / (8. * CONST_PI * CONST_E0 * CONST_E0);
    nu0 *= mccc_coefs_clog(0, 0, config);
    printf("NU0 %e\n", nu0);
    double t0 = v0 * v0 * v0 * CONST_ME * CONST_ME / (n0 * nu0);
    ret->dt /= t0;

    // double T0x = ret->species[0].Tx; // use first specie T as base for normalization
    // double T0y = ret->species[0].Ty; // use first specie T as base for normalization
    double T0 = CONST_ME * v0 * v0 / CONST_E;
    for (int s=0; s<ret->nspecies; s++) {
        ret->species[s].eps /= (v0*v0);
        printf("Specie %d eps %e eps_normalized %e\n", s, ret->species[s].eps*v0*v0, ret->species[s].eps);
        ret->species[s].Tx /= T0;
        ret->species[s].Ty /= T0;
        ret->species[s].m /= CONST_ME;
        ret->species[s].q /= CONST_E;
        ret->species[s].n /= n0;
        ret->species[s].xmin /= v0;
        ret->species[s].xmax /= v0;
        ret->species[s].ymin /= v0;
        ret->species[s].ymax /= v0;
        for (int i=0; i<ret->nspecies; i++) {
            ret->species[s].nu[i] /= nu0;
            ret->species[s].nu[i] = 1.; // FIXME ====
        }
        for (int i=0; i<ret->species[s].npeaks; i++) {
            ret->species[s].peaks[i] /= v0;
        }
    }

    ret->n0 = n0;
    ret->v0 = v0;
    ret->t0 = t0;
    ret->T0 = T0;
    
    return ret;
}

/**
 * @brief Deep copy the config object
 * 
 * @param config 
 * @return Config* 
 */
Config* copyConfig(Config* config) {
    Config* ret = (Config*)malloc(sizeof(Config));
    memcpy (ret, config, sizeof (Config));
    memcpy (ret->species, config->species, config->nspecies * sizeof (Specie));
    memcpy (ret->kHermite, config->kHermite, config->nHermite * sizeof (double));
    memcpy (ret->wHermite, config->wHermite, config->nHermite * sizeof (double));
    for (int i=0; i<ret->nspecies; i++) {
        memcpy (ret->species[i].nu, config->species[i].nu, config->nspecies * sizeof (double));
        memcpy (ret->species[i].name, config->species[i].name, 20 * sizeof (char));
    }
    return ret;
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

    double dt = config->dt;
    if (config->normalize) {
        dt *= config->t0;
    }
    fprintf(fout, "%d %d %e\n", config->nspecies, config->nmarkers, dt);
    fprintf(fout, "%e %e %e %e %e %e %e\n", E, (E-E0)/E0, P[0], P[1], (P[0]-P0[0])/P0[0], (P[1]-P0[1])/P0[1], distmin);
    for (int s=0; s<config->nspecies; s++) {
        E = Kspecie(p1, s, config);
        P = MomentumSpecie(p1, s, config);
        double T = TemperatureSpecie(p1, s, config);
        fprintf(fout, "%e %e %e %e %s\n", E, P[0], P[1], T, config->species[s].name);
    }
    for (int s=0; s<config->nspecies; s++) {
        for (int i=0; i<config->nmarkers; i++) {
            Vector2d z = p_mesh[s][i].z;
            if (config->normalize) {
                z *= config->v0;
            }
            fprintf(fout, "%d %d %e %e %e\n", s, i, z(0), z(1), f_mesh[s][i]);
        }
    }


    if (VERBOSE_LEVEL >= VERBOSE_SILLY) {
        for (int s=0; s<config->nspecies; s++) {
            for (int i=0; i<config->nmarkers; i++) {
                print_out(VERBOSE_SILLY, "Specie %d ID: %d Vx: %.15e Vy: %.15e W: %.15e\n", s, i, p1[s][i].z[0], p1[s][i].z[1], p1[s][i].weight);
            }
        }
    }

    fclose(fout);
}

/**
 * @brief Distribution function of initial states
 * 
 * @param v 
 * @param specie index of specie
 * @param config
 * @return double 
 */
double f(Vector2d v, int s, Config* config) {
    Specie specie = config->species[s];
    double ret = 0;
    double m = specie.m;
    double n = specie.n;
    double Tx = specie.Tx * CONST_E; // compute T in Kelvin
    double Ty = specie.Ty * CONST_E; // compute T in Kelvin
    if (config->normalize) {
        Tx = specie.Tx;
        Ty = specie.Ty;
    }
    for (int i=0; i<specie.npeaks; i++) {
        ret += exp(-(pow(v(0)-specie.peaks[i](0),2)/Tx + pow(v(1)-specie.peaks[i](1),2)/Ty) *m*0.5) / sqrt(Tx*Ty);
    }

    ret *= n*(specie.ymax-specie.ymin)*(specie.xmax-specie.xmin)/(config->nx*config->ny)*m/(CONST_2PI);
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
            double CONST_2EPS_M1 = 1./(2.*config->species[s].eps);
            double CONST_2PIEPS_M1 = 1./(CONST_2PI*config->species[s].eps);
            for (int j=0; j<config->nmarkers; j++) {
                ret[s][i] += exp(-(p_mesh[s][i].z - p[s][j].z).squaredNorm()*CONST_2EPS_M1)*p[s][j].weight;
            }
            ret[s][i] *= CONST_2PIEPS_M1;
            if (config->normalize) {
                ret[s][i] *= config->n0;
            }
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
Particle2d* initMarkers(int s, Config* config, DistributionType type) {
    int nmarkers = config->nmarkers;
    Particle2d* ret;
    ret = new Particle2d[config->nmarkers];
    int idx;
    Specie specie = config->species[s];
    for (int i = 0; i<config->ny; i++) {
        for (int j = 0; j<config->nx; j++) {
            idx = i*config->nx + j;
            if (type == UNIFORM) {
                ret[idx].z[0] = double(rand()) / RAND_MAX * (specie.xmax - specie.xmin) + specie.xmin;
                ret[idx].z[1] = double(rand()) / RAND_MAX * (specie.ymax - specie.ymin) + specie.ymin;
            } else if (type == MESH) {
                ret[idx].z[0] = double(j+0.5) / (config->nx) * (specie.xmax-specie.xmin) + specie.xmin;
                ret[idx].z[1] = double(i+0.5) / (config->nx) * (specie.ymax-specie.ymin) + specie.ymin;
            } else if (type == MESH_SHIFT) {
                double shift = 0;
                if (s == 1) {
                    shift = 0.5;
                }
                ret[idx].z[0] = double(j+shift) / (config->nx) * (specie.xmax-specie.xmin) + specie.xmin;
                ret[idx].z[1] = double(i+shift) / (config->nx) * (specie.ymax-specie.ymin) + specie.ymin;
            } else if (type == MESH_PEAK_CENTERED) {
                double dx = 0.1;
                Vector2d peak = specie.peaks[0];
                ret[idx].z[0] = double(j+0.5) / (config->nx) * (specie.xmax-specie.xmin) + peak(0) - dx/2;
                ret[idx].z[1] = double(i+0.5) / (config->ny) * (specie.ymax-specie.ymin) + peak(1) - dx/2;
            }
            ret[idx].weight = f(ret[idx].z, s, config);
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
 * @brief Compute density of a specie, a.k.a. integral of f over velocity
 * 
 * @param p 
 * @param s 
 * @param config 
 * @return double 
 */
double nSpecie(
    Particle2d** p,
    int s,
    Config* config
) {
    double rho = 0;
    for (int i=0; i<config->nmarkers; i++) {
        rho += p[s][i].weight;
    }
    return rho;
}

/**
 * @brief Compute specie temperature: 0.5*m*<(v-<v>)^2>
 * 
 * @param p 
 * @param s 
 * @param config 
 * @return double 
 */
double TemperatureSpecie(
    Particle2d** p,
    int s,
    Config* config
) {
    double T = 0;
    Vector2d V(0,0);

    double rho = nSpecie(p, s, config);

    // compute <v>
    for (int i=0; i<config->nmarkers; i++) {
        V += p[s][i].weight * p[s][i].z;
    }
    V /= rho;

    // compute Energy
    for (int i=0; i<config->nmarkers; i++) {
        T += p[s][i].weight * (p[s][i].z - V).squaredNorm();
    }
    T *= 0.5 * config->species[s].m / rho;
    if (config->normalize) {
        T *= config->T0; // compute T in real units [eV]
    }
    return T;
}

/**
 * @brief Compute the kinetic energy of a specie [eV]
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
    if (config->normalize) {
        ret *= config->n0 * CONST_ME * config->v0 * config->v0;
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
    if (config->normalize) {
        ret *= config->n0 * CONST_ME * config->v0;
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

/**
 * @brief Evaluate Coulomb logarithm.
 *
 * Coulomb logarithm is evaluated separately with respect to each plasma
 * species. It is calculated as a logarithm of the ratio of maximum and
 * minimum impact parameters. Maximum impact parameter is the Debye length
 * and minimum impact parameter is classical particle radius
 *
 * @param clogab array where evaluated values for Coulomb logarithm are stored.
 * @param ma test particle mass [kg]
 * @param qa test particle charge [C]
 * @param va test particle velocity [m/s]
 * @param nspec number of plasma species
 * @param mb plasma species masses [kg]
 * @param qb plasma species charges [C]
 * @param nb plasma species densities [m^-3]
 * @param Tb plasma species temperatures [J]
 */
double mccc_coefs_clog(int s1, int s2, Config* config) {

    /* Evaluate Debye length */
    double sum = 0;
    for(int i = 0; i < config->nspecies; i++){
        double qb = config->species[i].q;
        sum += config->species[i].n * qb * qb / sqrt(config->species[i].Tx*config->species[i].Ty) * CONST_E0;
    }
    double debyeLength = sqrt(CONST_E0/sum);

    /* Evaluate classical and quantum mechanical impact parameter. The one *
     * that is larger is used to evaluate Coulomb logarithm.               */
    double va = config->species[s1].peaks[0].squaredNorm();
    double qa = config->species[s1].q;
    double qb = config->species[s2].q;
    double Tb = sqrt(config->species[s2].Tx*config->species[s2].Ty) * CONST_E0;
    double ma = config->species[s1].m;
    double mb = config->species[s2].m;
    double vbar = va * va + 2 * Tb / mb;
    double mr   = ma * mb / ( ma + mb );
    double bcl  = fabs( qa * qb / ( 4*CONST_PI*CONST_E0 * mr * vbar ) );

    return log( debyeLength / bcl );
}

/**
 * @brief Evaluate collision parameter
 *
 *\f$c_{ab} = \frac{n_b q_a^2q_b^2 \ln\Lambda_{ab}}{4\pi\epsilon_0^2}\f$
 *
 * where
 *
 * - \f$q_a\f$ is test particle charge [C]
 * - \f$q_b\f$ is plasma species charge [C]
 * - \f$n_b\f$ is plasma species density [m^-3]
 * - \f$\ln\Lambda_{ab}\f$ is Coulomb logarithm.
 */
double mccc_coefs_cab(int s1, int s2, Config* config) {
    double qa = config->species[s1].q;
    double qb = config->species[s1].q;
    return qa*qa * qb*qb * mccc_coefs_clog(s1, s2, config) / ( 8 * CONST_PI * CONST_E0*CONST_E0 );
}

}