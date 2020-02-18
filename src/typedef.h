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

/* nucleation model type */
typedef enum {
   CNT_NUCLEATION,
   CONST_NUCLEATION,
   NONE_NUCLEATION
} nucleation_model_t;

/* chemical FE model type */
typedef enum {
   QUADRATIC_CHEMENERGY,
   CALPHAD_CHEMENERGY,
   NONE_CHEMENERGY
} chemfe_model_t;

/* interpolation type */
typedef enum {
   LINEAR_INTERPOLATION,
   QUADRATIC_INTERPOLATION,
   CUBIC_INTERPOLATION,
   NONE_INTERPOLATION
} interpolation_t;

/* CNT nucleation model container */
typedef struct CNT_NUC {
    PetscReal D0, migration, gamma, shapefactor;
    PetscReal minsize, atomicvolume, lengthscale;
} CNT_NUC;

/* constant nucleation model container */
typedef struct CONST_NUC {
    PetscReal nucleation_rate;
} CONST_NUC;

/* nucleation model */
typedef union NUC {
    CNT_NUC    cnt;
    CONST_NUC  constant;
} NUC;

/* Nucleus parameters */
typedef struct NUCLEUS {
    nucleation_model_t nuc_model;
    uint16_t matrixlist[MAXAP+1];
    NUC nucleation; 
} NUCLEUS;

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

/* Mobility */
typedef struct MOBILITY {
    PetscReal m0;
    TSeries *unary;
    RK *binary;
} MOBILITY;

/* CALPHAD energy parameters container */
typedef struct CALPHAD {
    TSeries ref;
    TSeries *unary;
    RK *binary, *ternary;
    MOBILITY *mobilityc;
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
    chemfe_model_t chemfe_model;
    PetscReal molarvolume, chempot_ex_kineticcoeff;
    PetscReal *c0;
    CHEMFE energy;
} MATERIAL;

/* interface container */
typedef struct INTERFACE {
    PetscReal energy;
    PetscReal *potential, *mobilityc;
    MOBILITY *mobility;
} INTERFACE;

/* AMR solparams */
typedef struct AMRPARAMS {
    PetscInt  initrefine, initcoarsen, maxnrefine, minnrefine, amrinterval, *initblocksize;
} AMRPARAMS;

/* solution solparams */
typedef struct SOLUTIONPARAMS {
    /* load case parameters */
    PetscReal *time, *temperature, timestep, mintimestep, maxtimestep;
    PetscInt  nloadcases, currentloadcase;
    PetscInt  step;
    /* phase field parameters */
    PetscReal interfacewidth;
    /* interpolation */
    interpolation_t interpolation;
    /* tolerances */
    PetscReal reltol, abstol;
    /* random seed */
    PetscInt  randomseed;
    /* output parameters */
    PetscInt  outputfreq;
    char      outfile[128];
    char      petscoptions[PETSC_MAX_PATH_LEN];
    PetscViewer viewer;
} SOLUTIONPARAMS;

/* solution parameters */
typedef struct AppCtx {
    /* MPI rank */
    PetscMPIInt worldrank, worldsize;
    /* number of phases, materials, interfaces and components */
    PetscInt npf, ndp, ncp, ntp;
    PetscInt nf, nmat;
    char **componentname, **materialname, **interfacename, **nucleusname;
    /* grid resolution */
    PetscInt dim, *resolution;
    PetscReal *size;
    /* time step */
    PetscInt step;
    /* aux grids and vecs */
    DM da_solution, da_solforest;
    DM da_output;
    PetscInt *localcells, nlocalcells, ninteriorcells;
    PetscInt *localfaces, nlocalfaces;
    PetscReal *cellgeom;
    /* phase material parameters */
    MATERIAL *material;
    PetscInt *voxelphasemapping, *phasematerialmapping;
    PetscInt *voxelsitemapping, *sitenucleusmapping, *sitephasemapping;
    /* nucleation parameters */
    PetscInt nsites, nsites_local, siteoffset, nnuclei;
    NUCLEUS *nucleus;
    char *siteactivity_local, *siteactivity_global;
    PetscSF nucleation_sf;
    /* interface material parameters */
    INTERFACE *interface;
    PetscInt *interfacelist;
    /* solution solparams */
    SOLUTIONPARAMS solparams;
    /* AMR solparams */
    AMRPARAMS amrparams;
    /* outputs */
    PetscInt noutputs;
    char **outputname;
} AppCtx;

#undef TYPEDEF_IMPORT
#undef EXTERN
#endif
