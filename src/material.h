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
EXTERN char NucleationEvent(const PetscReal, const PetscReal, const PetscReal, const PetscReal, const PetscReal, const PetscReal, const PetscInt, const AppCtx *);
EXTERN void NucleationBarrier(PetscReal *, PetscReal *, const PetscReal *, const PetscReal *, const PetscReal, const PetscReal *, const uint16_t, const uint16_t *, const AppCtx *);
EXTERN void Chemsource(PetscReal *, const PetscReal *, const PetscReal *, const uint16_t *, const AppCtx *);
EXTERN void Heatflux(PetscReal *, const PetscReal *, const uint16_t *, const AppCtx *);
EXTERN void Chemenergy(PetscReal *, const PetscReal *, const PetscReal *, const PetscReal, const uint16_t *, const AppCtx *);
EXTERN void SitepotentialExplicit(PetscReal *, const PetscReal *, const PetscReal, const uint16_t, const AppCtx *);
EXTERN void SitepotentialImplicit(PetscReal *, const PetscReal *, const PetscReal, const uint16_t, const AppCtx *);
EXTERN void Sitepotential(PetscReal *, const PetscReal *, const PetscReal, const uint16_t, const AppCtx *);
EXTERN void Sitefrac(PetscReal *, const PetscReal *, const PetscReal, const uint16_t, const AppCtx *);
EXTERN void SitefracTangent(PetscReal *, PetscReal *, PetscReal *, const PetscReal *, const PetscReal *, const PetscReal, const uint16_t, const AppCtx *);
EXTERN void CompositionMobilityComponent(PetscReal *, const PetscReal *, const PetscReal, const uint16_t, const AppCtx *);
EXTERN void CompositionMobilityVolumeRef(PetscReal *, const PetscReal *, const PetscReal *, const AppCtx *);
EXTERN void material_init(const AppCtx *);

#undef MATERIAL_IMPORT
#undef EXTERN
#endif