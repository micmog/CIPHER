/**
 * Initialisation header file
 */

#ifndef INIT_H
#define INIT_H

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <petscdm.h>
#include <petscdmda.h>
#include "material.h"
#include "utility.h"
#include "typedef.h"

/* Set EXTERN macro: */
#ifdef INIT_IMPORT
    #define EXTERN
#else
    #define EXTERN extern
#endif

/* Function prototypes */
EXTERN char *strtok_r(char *, const char *, char **);
EXTERN PetscErrorCode SetUpGeometry(AppCtx *);
EXTERN PetscErrorCode SetUpInterface(AppCtx *);
EXTERN PetscErrorCode SetUpProblem(DM,AppCtx *,Vec);

#undef INIT_IMPORT
#undef EXTERN
#endif