/**
 * Utility functions module
 */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <petscdm.h>
#include <petscdmda.h>
#include "typedef.h"


/* Including my own header for checking by compiler */
#define UTILITY_IMPORT
#include "utility.h"

/* nc x nc matrix inversion function */
static void (*Invertmatrixf)(PetscScalar *, PetscScalar *);

/* interpolation shape functions */
static void (*Interpolant          )(PetscScalar *, const PetscScalar *, const uint16_t);
static void (*InterpolantDerivative)(PetscScalar *, const PetscScalar *, const uint16_t);

/*
 Extract string between tags
 */
char *Extract(const char *const string, const char *const left, const char *const right)
{
    char  *head, *tail, *result;
    size_t length;

    if ((string == NULL) || (left == NULL) || (right == NULL))
        return NULL;
    length = strlen(left);
    head   = strstr(string, left);
    if (head == NULL)
        return NULL;
    head += length;
    tail  = strstr(head, right);
    if (tail == NULL)
        return tail;
    length = tail - head;
    result = malloc(1 + length);
    if (result == NULL)
        return NULL;
    result[length] = '\0';

    memcpy(result, head, length);
    return result;
}

/*
 PetscScalar to power of unsigned int
 */
PetscScalar FastPow(PetscScalar base, unsigned expn)
{
    PetscScalar result = 1.0;
    while (expn)
    {
        if (expn & 1)
            result *= base;
        expn >>= 1;
        base *= base;
    }

    return result;
}

/*
 PetscScalar fast square root
 */
PetscScalar FastSqrt(PetscScalar number)
{
    if (number < TOL) {
      return 0.0;
    } else {
      PetscScalar y = number;
      PetscScalar x2 = y * 0.5;
      int64_t i = *(int64_t *) &y;
      // The magic number is for doubles is from https://cs.uwaterloo.ca/~m32rober/rsqrt.pdf
      i = 0x5fe6eb50c7b537a9 - (i >> 1);
      y = *(PetscScalar *) &i;
      // 1st iteration
      y = y * (1.5 - (x2 * y * y));        
      // 2nd iteration, this can be removed
      y  = y * ( 1.5 - ( x2 * y * y ) );   
      return 1.0/y;
    }
}
  
/*
 PetscScalar fast exponential
 */
PetscScalar FastExp(PetscScalar a) {
  union { PetscScalar d; long long x; } u;
  u.x = (long long)(6497320848556798LL * a + 0x3fef127e83d16f12LL);
  return u.d;
}

/*
 PetscScalar fast log_e
 */
PetscScalar FastLog(PetscScalar a) {
  union { PetscScalar d; long long x; } u = { a };
  return (u.x - 4606921278410026770) * 1.539095918623324e-16; /* 1 / 6497320848556798.0; */
}

static void Invert1x1(PetscScalar * dst, PetscScalar * src)
{
    /* Compute adjoint and multiply with reciprocal of determinant: */

    dst[0] = 1.0 / src[0];
}

static void Invert2x2(PetscScalar * dst, PetscScalar * src)
{
    PetscScalar det;

    /* Compute adjoint: */

    dst[0] = + src[3];
    dst[2] = - src[2];
    dst[3] = + src[0];

    /* Compute determinant: */

    det = src[0] * dst[0] + src[1] * dst[2];

    /* Multiply adjoint with reciprocal of determinant: */

    det = 1.0 / det;

    dst[0] *= det;
    dst[2] *= det;
    dst[3] *= det;
    dst[1] = dst[2];
}

static void Invert3x3(PetscScalar * dst, PetscScalar * src)
{
    PetscScalar det;

    /* Compute adjoint: */

    dst[0] = + src[4] * src[8] - src[5] * src[7];
    dst[3] = - src[3] * src[8] + src[5] * src[6];
    dst[4] = + src[0] * src[8] - src[2] * src[6];
    dst[6] = + src[3] * src[7] - src[4] * src[6];
    dst[7] = - src[0] * src[7] + src[1] * src[6];
    dst[8] = + src[0] * src[4] - src[1] * src[3];

    /* Compute determinant: */

    det = src[0] * dst[0] + src[1] * dst[3] + src[2] * dst[6];

    /* Multiply adjoint with reciprocal of determinant: */

    det = 1.0 / det;

    dst[0] *= det;
    dst[3] *= det;
    dst[4] *= det;
    dst[6] *= det;
    dst[7] *= det;
    dst[8] *= det;
    dst[1] = dst[3];
    dst[2] = dst[6];
    dst[5] = dst[7];
}

static void Invert4x4(PetscScalar * dst, PetscScalar * src)
{
    PetscScalar det;

    /* Compute adjoint: */

    dst[0] =
        + src[ 5] * src[10] * src[15]
        - src[ 5] * src[11] * src[14]
        - src[ 9] * src[ 6] * src[15]
        + src[ 9] * src[ 7] * src[14]
        + src[13] * src[ 6] * src[11]
        - src[13] * src[ 7] * src[10];

    dst[4] =
        - src[ 4] * src[10] * src[15]
        + src[ 4] * src[11] * src[14]
        + src[ 8] * src[ 6] * src[15]
        - src[ 8] * src[ 7] * src[14]
        - src[12] * src[ 6] * src[11]
        + src[12] * src[ 7] * src[10];

    dst[5] =
        + src[ 0] * src[10] * src[15]
        - src[ 0] * src[11] * src[14]
        - src[ 8] * src[ 2] * src[15]
        + src[ 8] * src[ 3] * src[14]
        + src[12] * src[ 2] * src[11]
        - src[12] * src[ 3] * src[10];

    dst[8] =
        + src[ 4] * src[ 9] * src[15]
        - src[ 4] * src[11] * src[13]
        - src[ 8] * src[ 5] * src[15]
        + src[ 8] * src[ 7] * src[13]
        + src[12] * src[ 5] * src[11]
        - src[12] * src[ 7] * src[ 9];

    dst[9] =
        - src[ 0] * src[ 9] * src[15]
        + src[ 0] * src[11] * src[13]
        + src[ 8] * src[ 1] * src[15]
        - src[ 8] * src[ 3] * src[13]
        - src[12] * src[ 1] * src[11]
        + src[12] * src[ 3] * src[ 9];

    dst[10] =
        + src[ 0] * src[ 5] * src[15]
        - src[ 0] * src[ 7] * src[13]
        - src[ 4] * src[ 1] * src[15]
        + src[ 4] * src[ 3] * src[13]
        + src[12] * src[ 1] * src[ 7]
        - src[12] * src[ 3] * src[ 5];

    dst[12] =
        - src[ 4] * src[ 9] * src[14]
        + src[ 4] * src[10] * src[13]
        + src[ 8] * src[ 5] * src[14]
        - src[ 8] * src[ 6] * src[13]
        - src[12] * src[ 5] * src[10]
        + src[12] * src[ 6] * src[ 9];

    dst[13] =
        + src[ 0] * src[ 9] * src[14]
        - src[ 0] * src[10] * src[13]
        - src[ 8] * src[ 1] * src[14]
        + src[ 8] * src[ 2] * src[13]
        + src[12] * src[ 1] * src[10]
        - src[12] * src[ 2] * src[ 9];

    dst[14] =
        - src[ 0] * src[ 5] * src[14]
        + src[ 0] * src[ 6] * src[13]
        + src[ 4] * src[ 1] * src[14]
        - src[ 4] * src[ 2] * src[13]
        - src[12] * src[ 1] * src[ 6]
        + src[12] * src[ 2] * src[ 5];

    dst[15] =
        + src[ 0] * src[ 5] * src[10]
        - src[ 0] * src[ 6] * src[ 9]
        - src[ 4] * src[ 1] * src[10]
        + src[ 4] * src[ 2] * src[ 9]
        + src[ 8] * src[ 1] * src[ 6]
        - src[ 8] * src[ 2] * src[ 5];

    /* Compute determinant: */

    det = + src[0] * dst[0] + src[1] * dst[4] + src[2] * dst[8] + src[3] * dst[12];

    /* Multiply adjoint with reciprocal of determinant: */

    det = 1.0 / det;

    dst[ 0] *= det;
    dst[ 4] *= det;
    dst[ 5] *= det;
    dst[ 8] *= det;
    dst[ 9] *= det;
    dst[10] *= det;
    dst[12] *= det;
    dst[13] *= det;
    dst[14] *= det;
    dst[15] *= det;
    dst[ 1] = dst[ 4];
    dst[ 2] = dst[ 8];
    dst[ 3] = dst[12];
    dst[ 6] = dst[ 9];
    dst[ 7] = dst[13];
    dst[11] = dst[14];
}

static void Invertnxn(PetscScalar * dst, PetscScalar * src)
{
    PetscInt n = sizeof(src)/sizeof(PetscScalar);
    memset(dst,0,n*n*sizeof(PetscScalar));
    for(PetscInt i = 0; i < n; i++) dst[i*n+i] = 1.0;

    for(PetscInt i = 0; i < n; i++) {
        for(PetscInt j = 0; j < n; j++) {
            if (i != j) {
                PetscScalar ratio = src[j*n+i]/src[i*n+i];
                for(PetscInt k = 0; k < n; k++) src[j*n+k] -= ratio * src[i*n+k];
                for(PetscInt k = 0; k < n; k++) dst[j*n+k] -= ratio * dst[i*n+k];
            }
        }
    }

    for(PetscInt i = 0; i < n; i++) {
        PetscScalar a = src[i*n+i];
        for(PetscInt j = 0; j < n; j++) dst[i*n+j] /= a;
    }
}

void Invertmatrix(PetscScalar * dst, PetscScalar * src)
{
    Invertmatrixf(dst,src);
}

/*
 Linear interpolation shape function
 */
static void Interpolant_linear(PetscScalar *interpolant, const PetscScalar *phasefrac, const uint16_t nphases)
{
    memcpy(interpolant,phasefrac,nphases*sizeof(PetscScalar));
}

/*
 Linear interpolation shape function derivative
 */
static void InterpolantDerivative_linear(PetscScalar *interpolantderivative, const PetscScalar *phasefrac, const uint16_t nphases)
{
    memset(interpolantderivative,0,nphases*nphases*sizeof(PetscScalar));
    for (PetscInt g = 0; g < nphases; g++) interpolantderivative[g*nphases+g] = 1.0;
}

/*
 Cubic interpolation shape function
 */
static void Interpolant_cubic(PetscScalar *interpolant, const PetscScalar *phasefrac, const uint16_t nphases)
{
    PetscScalar interpolantsum = 0.0;
    for (PetscInt g = 0; g < nphases; g++) {
        interpolant[g] = phasefrac[g]*phasefrac[g]*(3.0 - 2.0*phasefrac[g]);
        interpolantsum += interpolant[g];
    }
    interpolantsum = 1.0/interpolantsum;
    for (PetscInt g = 0; g < nphases; g++) interpolant[g] *= interpolantsum;
}

/*
 Cubic interpolation shape function derivative
 */
static void InterpolantDerivative_cubic(PetscScalar *interpolantderivative, const PetscScalar *phasefrac, const uint16_t nphases)
{
    memset(interpolantderivative,0,nphases*nphases*sizeof(PetscScalar));
    PetscScalar interpolant[nphases], interpolantd[nphases], interpolantsum = 0.0;
    for (PetscInt g = 0; g < nphases; g++) {
        interpolant [g] = phasefrac[g]*phasefrac[g]*(3.0 - 2.0*phasefrac[g]);
        interpolantd[g] = 6.0*phasefrac[g]*(1.0 - phasefrac[g]);
        interpolantsum += interpolant[g];
    }    
    interpolantsum = 1.0/interpolantsum;
    PetscScalar interpolantsum2 = interpolantsum*interpolantsum;
    for (PetscInt gk = 0; gk < nphases; gk++) {
        interpolantderivative[gk*nphases+gk] += interpolantd[gk]*interpolantsum;
        for (PetscInt gj = 0; gj < nphases; gj++) {
            interpolantderivative[gk*nphases+gj] -= interpolant[gk]*interpolantd[gj]*interpolantsum2;
        }
    }    
}

/*
 Interpolation shape function
 */
void Shapefunc(PetscScalar *interpolant, const PetscScalar *phasefrac, const uint16_t nphases)
{
    Interpolant(interpolant,phasefrac,nphases);
}

/*
 Interpolation shape function derivative
 */
void ShapefuncDerivative(PetscScalar *interpolantderivative, const PetscScalar *phasefrac, const uint16_t nphases)
{
    InterpolantDerivative(interpolantderivative,phasefrac,nphases);
}

/*
 SimplexProjection - Project given vector on to Gibbs simplex
 */
static int Comparison(const void *x, const void *y){   
   if (*(PetscScalar*)y > *(PetscScalar*)x) {
       return 1;
   } else if (*(PetscScalar*)y < *(PetscScalar*)x) {    
       return -1;
   } else {
       return 0;
   }
}
void SimplexProjection(PetscScalar *out, PetscScalar *in, int size)
{
    char bget=0;
    PetscScalar s[size], tmax, tempsum=0.0;
    memcpy(s,in,size*sizeof(PetscScalar));
    qsort(s,size,sizeof(PetscScalar),Comparison);
    for (int i=0; i<size-1; i++) {
        tempsum += s[i];
        tmax = (tempsum-1.0)/((double)(i+1));
        if (tmax >= s[i+1]) {
            bget = 1;
            break;
        }
    }
    if (!bget) tmax = (tempsum + s[size-1] - 1.0)/((PetscScalar)(size));
    for (int i=0; i<size; i++) out[i] = in[i]-tmax > 0.0 ? in[i]-tmax : 0.0;
}

/*
 Initialise material module
 */
void utility_init(const AppCtx *user)
{
    if        (user->nc-1 == 1) {
        Invertmatrixf = &Invert1x1;
    } else if (user->nc-1 == 2) {
        Invertmatrixf = &Invert2x2;
    } else if (user->nc-1 == 3) {
        Invertmatrixf = &Invert3x3;
    } else if (user->nc-1 == 4) {
        Invertmatrixf = &Invert4x4;
    } else {
        Invertmatrixf = &Invertnxn;
    }  
    if        (user->interpolation == LINEAR_INTERPOLATION) {
        Interpolant = &Interpolant_linear;
        InterpolantDerivative = &InterpolantDerivative_linear;
    } else if (user->interpolation == CUBIC_INTERPOLATION ) {
        Interpolant = &Interpolant_cubic;
        InterpolantDerivative = &InterpolantDerivative_cubic;
    }
}
