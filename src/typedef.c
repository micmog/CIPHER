/**
 * Utility functions module
 */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <petscdm.h>
#include <petscdmda.h>


/* Including my own header for checking by compiler */
#define TYPEDEF_IMPORT
#include "typedef.h"

/* global variable definitions */
PetscInt MAXSITES = 0;
PetscInt SP_SIZE = 0;
PetscInt SF_SIZE = 0;
PetscInt AS_SIZE = 0;
PetscInt PF_SIZE = 0;
PetscInt DP_SIZE = 0;
PetscInt EX_SIZE = 0;
PetscInt TM_SIZE = 0;
PetscInt AS_OFFSET = 0;
PetscInt PF_OFFSET = 0;
PetscInt DP_OFFSET = 0;
PetscInt EX_OFFSET = 0;
PetscInt TM_OFFSET = 0;
