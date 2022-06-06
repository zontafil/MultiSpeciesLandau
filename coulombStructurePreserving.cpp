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
        config->n0 = 1.;
        config->v0 = 1.;
        config->t0 = 1.;
        config->T0 = 1.;
        config->nu0 = 1.;
    }

    // init markers
    Particle2d* p_mesh = initOutputPrintMesh(config);
    Particle2d** p0 = new Particle2d*[config->nspecies];
    Particle2d** p1 = new Particle2d*[config->nspecies];
    for (int i=0; i<config->nspecies; i++) {
        Specie specie = config->species[i];
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

    // init eq. motion iterations helper variables
    VectorXd* eom_f = new VectorXd[config->nspecies];
    for (int s=0; s<config->nspecies; s++) {
        eom_f[s] = VectorXd(2*config->nmarkers);
    }

    // initial energy
    double E0 = K(p1, config), E;
    double** f_mesh = new double*[config->nspecies];
    Vector2d P0 = Momentum(p1, config), P;
    VectorXd* dSdV = new VectorXd[config->nspecies];
    for (int s=0; s<config->nspecies; s++) {
        f_mesh[s] = new double[config->_nmarkers_outputmesh];
        dSdV[s] = VectorXd(2*config->nmarkers);
        print_out(VERBOSE_NORMAL, "Specie %d Epsilon: %e\n", s, config->species[s].eps);
    }

    print_out(VERBOSE_NORMAL, "Markers: %d dT: %e\n", config->nmarkers, config->dt);
    print_out(VERBOSE_NORMAL, "Initial Energy [eV]: %e\n", E0);
    print_out(VERBOSE_NORMAL, "Initial Momentum: %e %e\n", P0[0], P0[1]);
    print_out(VERBOSE_NORMAL, "Thermalization time: %e s\n", thermalizationTime(config0));

    double t = 0;
    int nsteps = 0, ETA;
    time_t cpuTime0 = time(NULL);
    char ETAstring[100];
    config->time_dsdv = 0;
    config->time_eqmotion = 0;
    config->time_total = -getTime();
    while (t<=config0->t1) {

        ETA = (time(NULL) - cpuTime0) * (config0->t1/t - 1.);
        format_duration(ETA, ETAstring);
        print_out(VERBOSE_NORMAL, "Timestep: %d Time: %e s ETA: %s\n", nsteps, t, ETAstring);

        for (int s=0; s<config->nspecies; s++) {
            copy(p1[s], p1[s]+config->nmarkers, p0[s]);
        }

        if (nsteps%config->recordAtStep == 0) {
            printState(f_mesh, p_mesh, p1, config, config0, nsteps, E0, P0);
        }

        // precompute entropy gradient
        config->time_dsdv -= getTime();
        Kernel::computedSdv(dSdV, p0, config);
        config->time_dsdv += getTime();
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
        config->time_eqmotion -= getTime();
        for (int j=0; j<config->maxEOMIterations; j++) {
            print_out(VERBOSE_DEBUG, "Iteration %d ", j);
            if (config->useNewton) {
                if (pushForwardNewtonIteration(p0, p1, dSdV, config)) {
                    break;
                }
            } else {
                if (pushForward_dv(p0, p1, dSdV, eom_f, config)) {
                    break;
                }
            }
        }
        config->time_eqmotion += getTime();

        t += config0->dt;
        nsteps++;
    }
    config->time_total += getTime();
    config->time_total *= config->dt / config->t1;
    config->time_eqmotion *= config->dt / config->t1;
    config->time_dsdv *= config->dt / config->t1;
    FILE* benchOut = fopen("benchmarks.txt", "a+");
    
    fprintf(benchOut, "%d %f %f %f\n", config->nx, config->time_total, config->time_dsdv, config->time_eqmotion);

    printState(f_mesh, p_mesh, p1, config, config0, nsteps, E0, P0);

    delete[] eom_f;
    for (int s=0; s<config->nspecies; s++) {
        free(p0[s]);
        free(p1[s]);
        free(f_mesh[s]);
    }
    free(p_mesh);
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
    double nu0 = coefs_nu(0, 0, config);
    printf("NU0 %e\n", nu0);
    double t0 = v0 * v0 * v0 * CONST_ME * CONST_ME / (n0 * nu0);
    ret->dt /= t0;
    ret->t1 /= t0;

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
        }
        for (int i=0; i<ret->species[s].npeaks; i++) {
            ret->species[s].peaks[i] /= v0;
        }
    }

    ret->n0 = n0;
    ret->v0 = v0;
    ret->t0 = t0;
    ret->T0 = T0;
    ret->nu0 = nu0;
    
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
    ret->species = new Specie[config->nspecies];
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
    Particle2d* p_mesh,
    Particle2d** p1,
    Config* config,
    Config* config0,
    int t,
    double E0,
    Vector2d P0
) {
    // build distribution at mesh nodes and print to file
    mesh_distribution(f_mesh, p_mesh, p1, config);
    char filename[30];
    snprintf (filename, sizeof filename, "out/data/step_C_%08d.txt", t);
    FILE* fout = fopen(filename, "w+");

    // print system state and debug info to screen
    double E = K(p1, config);
    print_out(VERBOSE_NORMAL, "Energy: %.15e Error: %.15e\n", E, (E-E0)/E0);

    Vector2d P = Momentum(p1, config), V;
    print_out(VERBOSE_NORMAL, "Momentum: %.15e %.15e\n", P[0], P[1]);
    print_out(VERBOSE_NORMAL, "Momentum Error: %.15e %.15e\n", (P[0]-P0[0])/P0[0], (P[1]-P0[1])/P0[1]);

    double thermalAnalytic, thermalTime, Ta, Tb;
    thermalTime = thermalizationTime(config0);
    Ta = sqrt(config0->species[0].Tx*config0->species[0].Ty);
    Tb = sqrt(config0->species[1].Tx*config0->species[1].Ty);
    thermalAnalytic = (Ta+Tb)/2. + (Ta-Tb)/2.* exp(-2.*t*config0->dt/thermalTime);

    double dt = config->dt;
    if (config->normalize) {
        dt *= config->t0;
    }
    fprintf(fout, "%d %d %d %e\n", config->nspecies, config->nmarkers, config->_nmarkers_outputmesh, dt);
    fprintf(fout, "%e %e %e %e %e %e %e\n", E, (E-E0)/E0, P[0], P[1], (P[0]-P0[0]), (P[1]-P0[1]), thermalAnalytic);
    double T, Tx, Ty;
    for (int s=0; s<config->nspecies; s++) {
        E = Kspecie(p1, s, config);
        V = averageVelocitySpecie(p1, s, config);
        T = TemperatureSpecie(p1, s, config);
        Tx = TemperatureSpecieSingleAxis(p1, s, 0, config);
        Ty = TemperatureSpecieSingleAxis(p1, s, 1, config);
        if (config->normalize) {
            V *= config->v0;
        }
        fprintf(fout, "%e %e %e %e %e %e %s\n", E, V[0], V[1], T, Tx, Ty, config->species[s].name);
    }

    if (t % config->recordMeshAtStep == 0) {
        for (int s=0; s<config->nspecies; s++) {
            for (int i=0; i<config->_nmarkers_outputmesh; i++) {
                Vector2d z = p_mesh[i].z;
                if (config->normalize) {
                    z *= config->v0;
                }
                fprintf(fout, "%d %d %e %e %e\n", s, i, z(0), z(1), f_mesh[s][i]);
            }
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
    double Tx, Ty;

    if (config->normalize) {
        Tx = specie.Tx;
        Ty = specie.Ty;
    } else {
        Tx = specie.Tx * CONST_E; // compute T in J
        Ty = specie.Ty * CONST_E; // compute T in J
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
    Particle2d* p_mesh,
    Particle2d** p,
    Config* config
) {
    for (int s=0; s<config->nspecies; s++) {
        for (int i=0; i<config->_nmarkers_outputmesh; i++) {
            ret[s][i] = 0;
            double CONST_2EPS_M1 = 1./(2.*config->species[s].eps);
            double CONST_2PIEPS_M1 = 1./(CONST_2PI*config->species[s].eps);
            for (int j=0; j<config->nmarkers; j++) {
                ret[s][i] += exp(-(p_mesh[i].z - p[s][j].z).squaredNorm()*CONST_2EPS_M1)*p[s][j].weight;
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
                ret[idx].z[1] = double(i+0.5) / (config->ny) * (specie.ymax-specie.ymin) + specie.ymin;
            } else if (type == MESH_SHIFT) {
                double shift = 0;
                if (s == 1) {
                    shift = 0.5;
                }
                ret[idx].z[0] = double(j+shift) / (config->nx) * (specie.xmax-specie.xmin) + specie.xmin;
                ret[idx].z[1] = double(i+shift) / (config->ny) * (specie.ymax-specie.ymin) + specie.ymin;
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

/**
 * @brief Init the mesh used to print the global output distribution used for plotting
 * The mesh covers all specie meshes and has the resolution of the mesh with the highest resolution
 * 
 * @param config 
 * @return Particle2d* 
 */
Particle2d* initOutputPrintMesh(Config* config) {

    // find limits of the mesh, find the largest limit of the species with the lowest separation
    double xmin = config->species[0].xmin;
    double xmax = config->species[0].xmax;
    double ymin = config->species[0].ymin;
    double ymax = config->species[0].ymax;
    double dx = (xmax - xmin) / config->nx;
    double dy = (ymax - ymin) / config->ny;
    for (int s=1; s<config->nspecies; s++) {
        xmin = fmin(xmin, config->species[s].xmin);
        ymin = fmin(ymin, config->species[s].ymin);
        xmax = fmax(xmax, config->species[s].xmax);
        ymax = fmax(ymax, config->species[s].ymax);
        if (config->highResolutionMesh) {
            dx = fmin(dx, (config->species[s].xmax - config->species[s].xmin)/config->nx);
            dy = fmin(dy, (config->species[s].ymax - config->species[s].ymin)/config->ny);
        } else {
            dx = fmax(dx, (config->species[s].xmax - config->species[s].xmin)/config->nx);
            dy = fmax(dy, (config->species[s].ymax - config->species[s].ymin)/config->ny);
        }
    }

    int nx = (xmax - xmin) / dx;
    int ny = (ymax - ymin) / dy;
    int nmarkers = nx*ny;
    config->_nmarkers_outputmesh = nmarkers;

    printf("INIT OUTPUT MESH: %d %d %e %e %e\n", nx, ny, xmin, xmax, dx);

    Particle2d* ret;
    ret = new Particle2d[nmarkers];
    int idx;
    for (int i = 0; i<ny; i++) {
        for (int j = 0; j<nx; j++) {
            idx = i*nx + j;
            ret[idx].z[0] = double(j+0.5) / nx * (xmax-xmin) + xmin;
            ret[idx].z[1] = double(i+0.5) / ny * (ymax-ymin) + ymin;
        }
    }
    return ret;
}

int pushForward_dv(
    Particle2d** p0,
    Particle2d** p1,
    VectorXd* dSdV,
    VectorXd* f,
    Config* config
) {
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
 * @brief Compute specie temperature of a single axis: 0.5*m*<(v-<v>)^2>
 * 
 * @param p 
 * @param s 
 * @param axis 
 * @param config 
 * @return double 
 */
double TemperatureSpecieSingleAxis(
    Particle2d** p,
    int s,
    int axis,
    Config* config
) {
    double T = 0;
    double V = 0;

    double rho = nSpecie(p, s, config);

    // compute <v>
    for (int i=0; i<config->nmarkers; i++) {
        V += p[s][i].weight * p[s][i].z(axis);
    }
    V /= rho;

    // compute Energy
    for (int i=0; i<config->nmarkers; i++) {
        T += p[s][i].weight * (p[s][i].z(axis) - V) * (p[s][i].z(axis) - V);
    }
    T *= config->species[s].m / rho;
    if (config->normalize) {
        T *= config->T0; // compute T in real units [eV]
    } else {
        T /= CONST_E;
    }
    return T;
}

/**
 * @brief Compute the average velocity (flow) of a specie
 * 
 * @param p 
 * @param s 
 * @param config 
 * @return Vector2d 
 */
Vector2d averageVelocitySpecie(
    Particle2d** p,
    int s,
    Config* config
) {
    Vector2d V(0,0);
    double rho = nSpecie(p, s, config);
    for (int i=0; i<config->nmarkers; i++) {
        V += p[s][i].weight * p[s][i].z;
    }
    V /= rho;

    return V;
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
    double rho = nSpecie(p, s, config);
    Vector2d V = averageVelocitySpecie(p, s, config);

    // compute Energy
    for (int i=0; i<config->nmarkers; i++) {
        T += p[s][i].weight * (p[s][i].z - V).squaredNorm();
    }
    T *= 0.5 * config->species[s].m / rho;
    if (config->normalize) {
        T *= config->T0; // compute T in real units [eV]
    } else {
        T /= CONST_E;
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
 */
double mccc_coefs_clog(int s1, int s2, Config* config) {

    /* Evaluate Debye length */
    double sum = 0;
    for(int i = 0; i < config->nspecies; i++){
        double qb = config->species[i].q;
        sum += config->species[i].n * qb * qb / sqrt(config->species[i].Tx*config->species[i].Ty);
    }
    // printf("sum %e \n", sum);
    double debyeLength = sqrt(CONST_E0*CONST_E/sum);

    /* Evaluate classical impact parameter */
    double va = config->species[s1].peaks[0].norm();
    double qa = config->species[s1].q;
    double qb = config->species[s2].q;
    double Tb = sqrt(config->species[s2].Tx*config->species[s2].Ty); // [eV]
    double ma = config->species[s1].m;
    double mb = config->species[s2].m;
    double vbar = va * va + 2 * Tb * CONST_E / mb;
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
double coefs_nu(int s1, int s2, Config* config) {
    double qa = config->species[s1].q;
    double qb = config->species[s1].q;
    double clogab = mccc_coefs_clog(s1, s2, config);
    return qa*qa * qb*qb * mccc_coefs_clog(s1, s2, config) / ( 8 * CONST_PI * CONST_E0*CONST_E0 );

}

double thermalizationTime(Config* config) {
    if (config->nspecies != 2) {
        print_out(VERBOSE_NORMAL, "Thermalization time not available with nspecies != 2\n");
        return 0;
    }
    double ma = config->species[0].m;
    double mb = config->species[1].m;
    double Ta = sqrt(config->species[0].Tx*config->species[0].Ty);
    double Tb = sqrt(config->species[1].Tx*config->species[1].Ty);
    double TaK = Ta * CONST_E / CONST_KB;
    double TbK = Tb * CONST_E / CONST_KB;
    double mt = sqrt(ma*Tb + mb*Ta);
    double mu = 2.;
    double na = config->species[0].n * 1E-6; // n_a in cm^-3
    double l = mccc_coefs_clog(0, 1, config);
    return (mt*mt*mt) /
            (1.8E-19 * sqrt(ma*mb)* na * l);
}

void format_duration(int seconds, char* ret) {
    int minutes = (int) ((seconds / (60)) % 60);
    int hours   = (int) ((seconds / (60*60)) % 24);

    sprintf(ret, "%d Hours, %d Minutes, %d Seconds", hours, minutes, seconds%60);
}

/**
 * @brief milliseconds since epoch
 * 
 * @return double 
 */
double getTime() {
    struct timespec spec;
    clock_gettime(CLOCK_REALTIME, &spec);
    return (spec.tv_sec*1000 + spec.tv_nsec * 1E-6);
}

}