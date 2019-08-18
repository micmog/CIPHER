/**
 * Type definitions header file
 */

#ifndef TYPEDEF_H
#define TYPEDEF_H

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <petscdm.h>
#include <petscdmda.h>
#include "roaring.h"

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
#define MAXTP 3
#define MAXCP 4
#define MAXAP 12

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
   CUBIC_INTERPOLATION,
   NONE_INTERPOLATION
} interpolation_t;

/* field dofs */
typedef struct F_DOFS {
   PetscScalar   p[MAXAP];
   PetscScalar   m[MAXCP];
} F_DOFS;

/* composition dofs */
typedef struct C_DOFS {
   PetscScalar   c[MAXAP*MAXCP];
} C_DOFS;

/* output dofs */
typedef struct O_DOFS {
   PetscScalar   ***o;
} O_DOFS;

/* Float to int conversion */
typedef union F2I {
   PetscScalar   f[MAXTP];
   uint16_t      i[MAXAP];
} F2I;

/* RK coefficients */
typedef struct RK {
    PetscInt n, *i;
    PetscScalar *enthalpy;
} RK;

/* CALPHAD energy parameters container */
typedef struct CALPHAD {
    PetscScalar ref, RT;
    PetscScalar *unary;
    RK *binary, *ternary;
    PetscScalar *mobilityc;
} CALPHAD;

/* Quadratic energy parameters container */
typedef struct QUAD {
    PetscScalar ref;
    PetscScalar *ceq;
    PetscScalar *unary, *binary;
    PetscScalar *mobilityc;
} QUAD;

/* Energy container */
typedef union CHEMFE {
   QUAD    quad;
   CALPHAD calphad;
} CHEMFE;

/* Phase container */
typedef struct MATERIAL {
    model_t model;
    PetscScalar molarvolume;
    PetscScalar *c0;
    CHEMFE energy;
} MATERIAL;

/* Phase container */
typedef struct INTERFACE {
    PetscScalar energy, mobility;
    PetscScalar *potential, *mobilityc;
} INTERFACE;

/* solution parameters */
typedef struct AppCtx {
    /* number of phases, materials, interfaces and components */
    PetscInt np, nmat, nf, nc;
    char **componentname;
    /* grid size and resolution */
    PetscReal len;
    PetscInt resolution[3];
    /* time step */
    PetscReal dt, dtmax, time;
    PetscInt step;
    /* tolerances */
    PetscReal ptol, ctol, error, error0, error1, kI;
    /* aux grids and vecs */
    DM daCmp, daOut, daIdx;
    Vec activephases, activeneighbourphases, rhs0, cvec;
    /* phase material parameters */
    MATERIAL *material;
    roaring_bitmap_t **phasevoxelmapping;
    uint16_t *phasematerialmapping;
    interpolation_t interpolation;
    /* interface material parameters */
    INTERFACE *interface;
    uint16_t *interfacelist;
    /* output */
    PetscInt outputfreq;
    char outfile[128];
} AppCtx;

#undef TYPEDEF_IMPORT
#undef EXTERN
#endif
