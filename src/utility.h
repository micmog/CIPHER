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
EXTERN void F2IFUNC(uint16_t *,PetscScalar *);
EXTERN void I2FFUNC(PetscScalar *,uint16_t *);
EXTERN char *Extract(const char *const, const char *const, const char *const);
EXTERN PetscReal FastPow(PetscReal, unsigned);
EXTERN PetscReal FastSqrt(PetscReal);
EXTERN PetscReal FastExp(PetscReal);
EXTERN PetscReal FastLog(PetscReal);
EXTERN void SetIntersection(uint16_t *, uint16_t *, uint16_t *, uint16_t *, uint16_t *);
EXTERN void SetUnion(uint16_t *, uint16_t *, uint16_t *, uint16_t *, uint16_t *);
EXTERN void Invertmatrix(PetscReal *, PetscReal *, const uint16_t);
EXTERN void MatMatMult_CIPHER(PetscReal *, const PetscReal *, const PetscReal *, const uint16_t);
EXTERN void MatVecMult_CIPHER(PetscReal *, const PetscReal *, const PetscReal *, const uint16_t);
EXTERN void EvalInterpolant(PetscReal *, const PetscReal *, const uint16_t);
EXTERN void MatMulInterpolantDerivative(PetscReal *,  const PetscReal *, const uint16_t);
EXTERN void SimplexProjection(PetscReal *, PetscReal *, int);
EXTERN void utility_init(const AppCtx *);

#undef UTILITY_IMPORT
#undef EXTERN
#endif