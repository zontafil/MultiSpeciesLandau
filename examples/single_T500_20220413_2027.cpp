#include "coulombStructurePreserving.h"

Coulomb::Config buildConfig() {
    Coulomb::Config config;

    double L = 10;
    int MARKERS_PER_DIM = 20;

    config.dx = 1E-4; // dx for finite difference derivative
    config.dt = 1/16.;
    config.n_timesteps = 200; // number of timesteps
    config.newtonTolerance = 1E-14; // minimum target for error of eq. of motions
    config.useNewton = 0;
    config.maxEOMIterations = 20;
    config.h = 2.*L/MARKERS_PER_DIM;
    config.eps = 0.64*pow(config.h, 1.98);
    config.xmin = -L;
    config.xmax = L;
    config.ymin = -L;
    config.ymax = L;
    config.distributionType = UNIFORM;
    config.nx = MARKERS_PER_DIM;
    config.ny = MARKERS_PER_DIM;
    config.distributionType = MESH;
    // config.recordAtStep = max(config.n_timesteps - 1, 1);
    config.recordAtStep = 1;

    // initial peaks config
    Vector2d u1 = Vector2d(-2, 1);
    Vector2d u2 = Vector2d(0, -1);
    config.nspecies = 1;
    config.species = new Specie[config.nspecies];
    config.species[0].m = 1;
    config.species[0].npeaks = 2;
    config.species[0].peaks = new Vector2d[2];
    config.species[0].peaks[0] = u1;
    config.species[0].peaks[1] = u2;
    config.species[0].nu = new double[config.nspecies];
    config.species[0].nu[0] = 1;
    sprintf(config.species[0].name, "Single specie");
    
    
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
    config.cudaThreadsPerBlock = 1024;

    return config;
}

int main() {
    Coulomb::Config config = buildConfig();
    Coulomb::Run(&config);
    return 0;
}