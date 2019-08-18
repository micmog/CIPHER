/**
 * Utility functions header file
 */

#ifndef UTILITY_H
#define UTILITY_H

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <petscdm.h>
#include <petscdmda.h>
#include "typedef.h"

/* Set EXTERN macro: */
#ifdef UTILITY_IMPORT
    #define EXTERN
#else
    #define EXTERN extern
#endif

/* Function prototypes */
EXTERN char *Extract(const char *const, const char *const, const char *const);
EXTERN PetscScalar FastPow(PetscScalar, unsigned);
EXTERN PetscScalar FastSqrt(PetscScalar);
EXTERN PetscScalar FastExp(PetscScalar);
EXTERN PetscScalar FastLog(PetscScalar);
EXTERN void Invertmatrix(PetscScalar *, PetscScalar *);
EXTERN void EvalInterpolant(PetscScalar *, const PetscScalar *, const uint16_t);
EXTERN void MatMulInterpolantDerivative(PetscScalar *, const PetscScalar *,  const PetscScalar *, const uint16_t);
EXTERN void SimplexProjection(PetscScalar *, PetscScalar *, int);
EXTERN void utility_init(const AppCtx *);

#undef UTILITY_IMPORT
#undef EXTERN
#endif