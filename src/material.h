/**
 * Material header file
 */

#ifndef MATERIAL_H
#define MATERIAL_H

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <petscdm.h>
#include <petscdmda.h>
#include "utility.h"
#include "typedef.h"

/* Set EXTERN macro: */
#ifdef MATERIAL_IMPORT
    #define EXTERN
#else
    #define EXTERN extern
#endif

/* Constants declarations */

/* Types declarations */

/* Global variables declarations */

/* Function prototypes */
EXTERN void Chemicalpotential  (PetscScalar *, const PetscScalar *, const PetscScalar *, const uint16_t *, const AppCtx *);
EXTERN void Chemenergy         (PetscScalar *, const PetscScalar *, const PetscScalar *, const uint16_t *, const AppCtx *);
EXTERN int  Composition        (PetscScalar *, const PetscScalar *,                      const uint16_t *, const AppCtx *);
EXTERN void CompositionTangent (PetscScalar *, const PetscScalar *, const PetscScalar *, const uint16_t *, const AppCtx *);
EXTERN void CompositionMobility(PetscScalar *, const PetscScalar *, const PetscScalar *, const uint16_t *, const AppCtx *);
EXTERN void material_init(const AppCtx *);

#undef MATERIAL_IMPORT
#undef EXTERN
#endif