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
EXTERN void Chemenergy(PetscReal *, const PetscReal *, const PetscReal *, const PetscReal, const uint16_t *, const AppCtx *);
EXTERN void ChemicalpotentialExplicit(PetscReal *, const PetscReal *, const PetscReal, const uint16_t, const AppCtx *);
EXTERN void ChemicalpotentialExplicitTangent(PetscReal *, const PetscReal *, const PetscReal, const uint16_t, const AppCtx *);
EXTERN void ChemicalpotentialImplicit(PetscReal *, const PetscReal *, const PetscReal, const uint16_t, const AppCtx *);
EXTERN void ChemicalpotentialImplicitTangent(PetscReal *, const PetscReal *, const PetscReal, const uint16_t, const AppCtx *);
EXTERN void Chemicalpotential(PetscReal *, const PetscReal *, const PetscReal, const uint16_t, const AppCtx *);
EXTERN void ChemicalpotentialTangent(PetscReal *, const PetscReal *, const PetscReal, const uint16_t, const AppCtx *);
EXTERN void Composition(PetscReal *, const PetscReal *, const PetscReal, const uint16_t, const AppCtx *);
EXTERN void CompositionTangent(PetscReal *, const PetscReal *, const PetscReal, const uint16_t, const AppCtx *);
EXTERN void CompositionMobility(PetscReal *, const PetscReal *, const uint16_t, const AppCtx *);
EXTERN void material_init(const AppCtx *);

#undef MATERIAL_IMPORT
#undef EXTERN
#endif