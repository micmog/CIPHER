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
#define TOL   1.0e-8
/* max active phases... depends on neighbourhood */
#define MAXFP 4 
#define MAXIP 4*MAXFP
#define MAXCP 5
#define MAXAP 15

/* field offsets */
#define AS_OFFSET 0 
#define AS_SIZE   MAXFP 
#define PF_OFFSET AS_OFFSET+AS_SIZE
#define PF_SIZE   MAXAP 
#define DP_OFFSET PF_OFFSET+PF_SIZE 
#define DP_SIZE   MAXCP 
#define CP_OFFSET DP_OFFSET+DP_SIZE 
#define CP_SIZE   MAXAP*MAXCP 

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

/* phase field dofs */
typedef struct PFIELD {
   PetscReal   p[PF_SIZE];
} PFIELD;

/* field dofs */
typedef struct DFIELD {
   PetscReal   d[DP_SIZE];
} DFIELD;

/* composition dofs */
typedef struct STATE {
   PetscReal   c[MAXAP*MAXCP];
} STATE;

/* output dofs */
typedef struct O_DOFS {
   PetscReal   p;
   PetscReal   c[MAXCP];
} O_DOFS;

/* Float to int conversion */
typedef union F2I {
   PetscReal   f[MAXFP];
   uint16_t    i[MAXIP];
} F2I;

/* RK coefficients */
typedef struct RK {
    PetscInt n, *i;
    PetscReal *enthalpy;
} RK;

/* CALPHAD energy parameters container */
typedef struct CALPHAD {
    PetscReal ref, RT;
    PetscReal *unary;
    RK *binary, *ternary;
    PetscReal *mobilityc;
} CALPHAD;

/* Quadratic energy parameters container */
typedef struct QUAD {
    PetscReal ref;
    PetscReal *ceq;
    PetscReal *unary, *binary;
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
    PetscReal molarvolume;
    PetscReal *c0;
    CHEMFE energy;
} MATERIAL;

/* interface container */
typedef struct INTERFACE {
    PetscReal energy, mobility;
    PetscReal *potential, *mobilityc;
} INTERFACE;

/* solution params */
typedef struct SOLUTIONPARAMS {
    /* time parameters */
    PetscReal finaltime, timestep, mintimestep, maxtimestep;
    PetscInt  step;
    /* phase field parameters */
    PetscReal interfacewidth;
    /* AMR parameters */
    PetscInt  initrefine, maxnrefine, amrinterval;
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
    PetscInt np, nmat, nf, nc;
    char **componentname;
    /* grid resolution */
    PetscInt resolution[3];
    PetscReal size[3];
    /* time step */
    PetscInt step;
    /* exception flag */
    PetscErrorCode rejectstage;
    /* aux grids and vecs */
    DM da_solution, da_solforest;
    DM da_matstate, da_matforest;
    Vec matstate;
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
    /* solution params */
    SOLUTIONPARAMS params;
} AppCtx;

#undef TYPEDEF_IMPORT
#undef EXTERN
#endif
