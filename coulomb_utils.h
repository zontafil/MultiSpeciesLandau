#ifndef COULOMB_UTILS_H
#define COULOMB_UTILS_H

typedef double real;  /**< Double precision float   */
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <limits>
#include <stdio.h>
#include <time.h>

#ifdef INTELLISENSE
#include <eigen3/Eigen/Dense> // global installation
#else
#include "Eigen/Dense" // local installation in same folder. FIX this
// Aalto workstation has Eigen as module but it's slightly outdated and has compile warnings
// CSC has no Eigen module at all
// #include <eigen3/Eigen/Dense> // global installation
#endif

using namespace std;
using namespace Eigen;

#define CONST_PI     3.1415926535897932384626
#define CONST_2PI    6.2831853071795862319959
#define PI2E0_5 2.50662827463
#define CONST_KB 1.3807e-23
#define CONST_E 1.602176565e-19
#define CONST_ME 9.1093837E-31
#define CONST_MD 3.3435837724E-27 // deuteron mass

/** @brief Electric constant [m^-3*kg^-1*s^4*A^2]                             */
#define CONST_E0 8.8542e-12

enum {
    VERBOSE_SILLY = 3,
    VERBOSE_DEBUG = 2,
    VERBOSE_NORMAL = 1,
    VERBOSE_MINIMAL = 0,
    VERBOSE_IO = 1
};

/**
 * @brief Print to standard output
 */
#define print_out(v,...) { if(VERBOSE_LEVEL >= (v)) printf(__VA_ARGS__); }

namespace Coulomb {
    typedef struct {
        double weight;
        Vector2d z;
    } Particle2d;

    typedef struct {
        int npeaks;
        Vector2d* peaks;
        double Tx;
        double Ty;
        double m;
        double* nu;
        double q;
        double n;
        double xmin;
        double xmax;
        double ymin;
        double ymax;
        char name[20];
        double eps;
    } Specie;

    typedef enum {
        MESH,
        UNIFORM,
        MESH_SHIFT,
        MESH_PEAK_CENTERED
    } DistributionType;

    typedef struct {
        int nx;
        int ny;
        int nmarkers;
        DistributionType distributionType;
        double dt;
        double dx;
        double t1;
        double newtonTolerance;
        int useNewton;
        int maxEOMIterations;
        double* kHermite;
        double* wHermite;
        int nHermite;
        int cudaThreadsPerBlock;

        // print to file parameters
        int recordAtStep; // print state every n-th step
        int recordMeshAtStep; // print mesh every n-th step
        int highResolutionMesh; // flag, use the highest resolution for building the output plot mesh

        // normalization constants
        int normalize; // normalize flag. Recommended to be set
        double n0;
        double v0;
        double t0;
        double T0;
        double nu0;

        // benchmark times
        double time_dsdv;
        double time_eqmotion;
        double time_total;

        int nspecies;
        Specie* species;

        // internal variables
        int _nmarkers_outputmesh;

    } Config;
}

#endif