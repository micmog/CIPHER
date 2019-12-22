/**
 * Type definitions header file
 */

#ifndef TYPEDEF_H
#define TYPEDEF_H

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>

/* Set EXTERN macro: */
#ifdef TYPEDEF_IMPORT
    #define EXTERN
#else
    #define EXTERN extern
#endif

/* Constants declarations */

#define LARGE 1.0e+6
#define TOL   1.0e-6
/* max active phases... depends on neighbourhood */
#define MAXCP 5
#define MAXAP 15
#define R_GAS_CONST 8.314462

/* field offsets */
PetscInt AS_SIZE,   PF_SIZE,   DP_SIZE,   CP_SIZE;
PetscInt AS_OFFSET, PF_OFFSET, DP_OFFSET, CP_OFFSET;

/* Types declarations */

/* model type */
typedef enum {
   QUADRATIC_CHEMENERGY,
   CALPHAD_CHEMENERGY,
   NONE_CHEMENERGY
} model_t;

/* interpolation type */
typedef enum {
   LINEAR_INTERPOLATION,
   QUADRATIC_INTERPOLATION,
   CUBIC_INTERPOLATION,
   NONE_INTERPOLATION
} interpolation_t;

/* Temperature series */
typedef struct TSeries {
    PetscInt  nTser;
    PetscInt  exp[10];
    PetscReal coeff[10];
    PetscReal logCoeff;
} TSeries;

/* RK coefficients */
typedef struct RK {
    PetscInt n, *i;
    TSeries *enthalpy;
} RK;

/* CALPHAD energy parameters container */
typedef struct CALPHAD {
    PetscReal ref;
    TSeries *unary;
    RK *binary, *ternary;
    PetscReal *mobilityc;
} CALPHAD;

/* Quadratic energy parameters container */
typedef struct QUAD {
    TSeries ref;
    TSeries *ceq;
    TSeries *unary, *binary;
    PetscReal *mobilityc;
} QUAD;

/* Energy container */
typedef union CHEMFE {
   QUAD    quad;
   CALPHAD calphad;
} CHEMFE;

/* Phase container */
typedef struct MATERIAL {
    model_t model;
    PetscReal molarvolume, statekineticcoeff;
    PetscReal *c0;
    CHEMFE energy;
} MATERIAL;

/* interface container */
typedef struct INTERFACE {
    PetscReal energy, mobility;
    PetscReal *potential, *mobilityc;
} INTERFACE;

/* AMR solparams */
typedef struct AMRPARAMS {
    PetscInt  initrefine, initcoarsen, maxnrefine, minnrefine, amrinterval, *initblocksize;
} AMRPARAMS;

/* solution solparams */
typedef struct SOLUTIONPARAMS {
    /* time parameters */
    PetscReal finaltime, timestep, mintimestep, maxtimestep;
    PetscInt  step;
    /* phase field parameters */
    PetscReal interfacewidth;
    /* temperature */
    PetscReal temperature;
    /* tolerances */
    PetscReal reltol, abstol;
    /* output parameters */
    PetscInt  outputfreq;
    char      outfile[128];
    char      petscoptions[PETSC_MAX_PATH_LEN];
    PetscViewer viewer;
} SOLUTIONPARAMS;

/* solution parameters */
typedef struct AppCtx {
    /* number of phases, materials, interfaces and components */
    PetscInt npf, ndp, ncp, ntp;
    PetscInt nf, nmat;
    char **componentname;
    /* grid resolution */
    PetscInt dim, *resolution;
    PetscReal *size;
    /* time step */
    PetscInt step;
    /* exception flag */
    PetscErrorCode rejectstage;
    /* aux grids and vecs */
    DM da_solution, da_solforest;
    DM da_output;
    PetscInt *localcells, nlocalcells, ninteriorcells;
    PetscInt *localfaces, nlocalfaces;
    PetscReal *cellgeom;
    /* phase material parameters */
    MATERIAL *material;
    uint16_t *phasevoxelmapping, *phasematerialmapping;
    interpolation_t interpolation;
    /* interface material parameters */
    INTERFACE *interface;
    uint16_t *interfacelist;
    /* solution solparams */
    SOLUTIONPARAMS solparams;
    /* AMR solparams */
    AMRPARAMS amrparams;
} AppCtx;

#undef TYPEDEF_IMPORT
#undef EXTERN
#endif
