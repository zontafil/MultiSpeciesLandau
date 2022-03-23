#ifndef COULOMB_UTILS_H
#define COULOMB_UTILS_H

typedef double real;  /**< Double precision float   */
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <limits>
#include <stdio.h>

#ifdef INTELLISENSE
#include <eigen3/Eigen/Dense> // global installation
#else
#include "Eigen/Dense" // local installation in same folder. FIX this
// Aalto workstation has Eigen as module but it's slightly outdated and has compile warnings
// CSC has no Eigen module at all
#endif

using namespace std;
using namespace Eigen;

#define CONST_PI     3.1415926535897932384626
#define CONST_2PI    6.2831853071795862319959
#define PI2E0_5 2.50662827463

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

typedef struct {
    double xmin;
    double xmax;
    int nx;
    double ymin;
    double ymax;
    int ny;
    int nmarkers;
    double dt;
    double dx;
    double nu;
    double m;
    double h;
    double eps;
    int n_timesteps;
    double newtonTolerance;
    int useNewton;
    int maxEOMIterations;
    double* kHermite;
    double* wHermite;
    int nHermite;
    int recordAtStep;
    int cudaThreadsPerBlock;
    Vector2d u1;
    Vector2d u2;
} Config;

typedef struct {
    double weight;
    Vector2d z;
} Particle2d;

#endif