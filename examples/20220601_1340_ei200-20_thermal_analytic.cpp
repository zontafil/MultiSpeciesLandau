#include "coulombStructurePreserving.h"

Coulomb::Config buildConfig() {
    Coulomb::Config config;

    int MARKERS_PER_DIM = 20;

    double n = 10E19;

    config.dt = 1E-8;
    config.t1 = 1E-4;
    config.newtonTolerance = 1E-14; // minimum target for error of eq. of motions
    config.useNewton = 0;
    config.maxEOMIterations = 20;
    config.nx = MARKERS_PER_DIM;
    config.ny = MARKERS_PER_DIM;
    config.distributionType = MESH_SHIFT;
    config.recordAtStep = 100;

    // initial peaks config
    Vector2d u1;
    config.nspecies = 2;
    config.species = new Specie[config.nspecies];
    double L, h;
    double k = 2;

    // SPECIE 1
    config.species[0].Tx = 400.;
    config.species[0].Ty = 300.;
    config.species[0].m = CONST_ME;
    L = sqrt(CONST_E * sqrt(config.species[0].Tx*config.species[0].Ty)/config.species[0].m) * 5;
    h = L/MARKERS_PER_DIM;
    config.species[0].q = -CONST_E;
    config.species[0].n = n;
    config.species[0].npeaks = 1;
    config.species[0].peaks = new Vector2d[1];
    // u1 = Vector2d(L/1000., 0);
    u1 = Vector2d(0, 0);
    config.species[0].peaks[0] = u1;
    config.species[0].xmin = -L;
    config.species[0].xmax = L;
    config.species[0].ymin = -L;
    config.species[0].ymax = L;
    config.species[0].eps = k*0.64*pow(h, 1.98);
    sprintf(config.species[0].name, "Electrons");

    // SPECIE 2
    if (config.nspecies > 1) {
        config.species[1].m = 200.*CONST_ME;
        config.species[1].Tx = 250;
        config.species[1].Ty = 200;
        L = sqrt(CONST_E * sqrt(config.species[1].Tx*config.species[1].Ty)/config.species[1].m) * 5;
        h = L/MARKERS_PER_DIM;
        config.species[1].q = CONST_E;
        config.species[1].n = n;
        config.species[1].npeaks = 1;
        config.species[1].peaks = new Vector2d[1];
        u1 = Vector2d(0., 0);
        config.species[1].peaks[0] = u1;
        config.species[1].xmin = -L;
        config.species[1].xmax = L;
        config.species[1].ymin = -L;
        config.species[1].ymax = L;
        config.species[1].eps = k*0.64*pow(h, 1.98);
        sprintf(config.species[1].name, "Ions");
    }

    // config.species[0].eps = 100.*2.055463e+11;
    // config.species[1].eps = 100.*2.055463e+11;

    printf("L %e\n", L);
    // config.h = 2.*L/MARKERS_PER_DIM / config.nspecies;
    // config.h = 2.*L/MARKERS_PER_DIM;
    printf("h %e\n", h);
    // config.eps = 40*2.055463e+13;
    // config.eps = 10*pow(config.h, 2.);
    
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
    config.dx = 1E-4; // dx for finite difference derivative
    config.cudaThreadsPerBlock = 1024;
    config.normalize = 1;
    config.highResolutionMesh = 0;


    // init coulomb logarithm and nu parameters
    for (int s1 = 0; s1<config.nspecies; s1++) {
        config.species[s1].nu = new double[config.nspecies];
        for (int s2 = 0; s2<=s1; s2++) {
            config.species[s1].nu[s2] = coefs_nu(s1, s2, &config);
            config.species[s2].nu[s1] = coefs_nu(s1, s2, &config);
            // if (s2<s1) {
            //     config.species[s2].nu[s1] *= 10.;
            //     config.species[s1].nu[s2] *= 10.;
            // }
        }
    }

    return config;
}

int main() {
    Coulomb::Config config = buildConfig();
    Coulomb::Run(&config);
    return 0;
}