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

/* interpolation shape functions */
static void (*Interpolant          )(PetscReal *, const PetscReal *, const uint16_t);
static void (*InterpolantDerivative)(PetscReal *, const PetscReal *, const uint16_t);

/*
 Set conversion routines
 */
void F2IFUNC(uint16_t *IARRAY,PetscScalar *FARRAY) {
    IARRAY[0] = (uint16_t) round(FARRAY[0]);
    for (int MACROIDX = 1; MACROIDX <= IARRAY[0]; ++MACROIDX) {
        IARRAY[MACROIDX] = (uint16_t) round(FARRAY[MACROIDX]);
    }
}
void I2FFUNC(PetscScalar *FARRAY,uint16_t *IARRAY) {
    FARRAY[0] = (PetscScalar) IARRAY[0];
    for (int MACROIDX = 1; MACROIDX <= IARRAY[0]; ++MACROIDX) {
        FARRAY[MACROIDX] = (PetscScalar) IARRAY[MACROIDX];
    }
}

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
 PetscReal to power of unsigned int
 */
PetscReal FastPow(PetscReal base, unsigned expn)
{
    PetscReal result = 1.0;
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
 PetscReal fast square root
 */
PetscReal FastSqrt(PetscReal number)
{
    if (number < TOL) {
      return 0.0;
    } else {
      PetscReal y = number;
      PetscReal x2 = y * 0.5;
      int64_t i = *(int64_t *) &y;
      // The magic number is for doubles is from https://cs.uwaterloo.ca/~m32rober/rsqrt.pdf
      i = 0x5fe6eb50c7b537a9 - (i >> 1);
      y = *(PetscReal *) &i;
      // 1st iteration
      y = y * (1.5 - (x2 * y * y));        
      // 2nd iteration, this can be removed
      y  = y * ( 1.5 - ( x2 * y * y ) );   
      return 1.0/y;
    }
}
  
/*
 PetscReal fast exponential
 */
PetscReal FastExp(PetscReal a) {
  union { PetscReal d; long long x; } u;
  u.x = (long long)(6497320848556798LL * a + 0x3fef127e83d16f12LL);
  return u.d;
}

/*
 PetscReal fast log_e
 */
PetscReal FastLog(PetscReal a) {
  union { PetscReal d; long long x; } u = { a };
  return (u.x - 4606921278410026770) * 1.539095918623324e-16; /* 1 / 6497320848556798.0; */
}

/*
 SetIntersection - C <-- A \cap B, A(CInA(i)) = C(i), B(CInB(i)) = C(i)
 */
void SetIntersection(uint16_t *OutC, uint16_t *CInA, uint16_t *CInB, uint16_t *InA, uint16_t *InB)
{
    uint16_t i=0, j=0;
    OutC[0] = 0;
    while (i < InA[0] && j < InB[0])
    {
        if (InA[i+1] == InB[j+1])
        {
            CInA[OutC[0]] = i; 
            CInB[OutC[0]] = j; 
            i++; j++; (OutC[0])++; 
            OutC[OutC[0]] = InA[i]; 
        }
        else if (InA[i+1] < InB[j+1])
        {
            i++;
        }
        else
        {
            j++;
        }
    }
}

/*
 SetUnion - C <-- A U B, C(AInC(i)) = A(i), C(AInB(i)) = B(i) 
 */
void SetUnion(uint16_t *OutC, uint16_t *AInC, uint16_t *BInC, uint16_t *InA, uint16_t *InB)
{
    uint16_t i=0, j=0;
    OutC[0] = 0;
    while (i < InA[0] && j < InB[0])
    {
        if      (InA[i+1] == InB[j+1])
        {
            AInC[i] = OutC[0]; 
            BInC[j] = OutC[0]; 
            i++; j++; (OutC[0])++;
            OutC[OutC[0]] = InA[i];
        }
        else if (InA[i+1] <  InB[j+1])
        {
            AInC[i] = OutC[0]; 
            i++;      (OutC[0])++;
            OutC[OutC[0]] = InA[i]; 
        }
        else
        {
            BInC[j] = OutC[0]; 
            j++;      (OutC[0])++; 
            OutC[OutC[0]] = InB[j]; 
        }
    }
    
    while (i < InA[0]) {
        AInC[i] = OutC[0];
        OutC[++(OutC[0])] = InA[++i];
    }
    while (j < InB[0]) {
        BInC[j] = OutC[0];
        OutC[++(OutC[0])] = InB[++j];
    }
}

static void Invert1x1(PetscReal * dst, PetscReal * src, const uint16_t n)
{
    /* Compute adjoint and multiply with reciprocal of determinant: */

    dst[0] = 1.0 / src[0];
}

static void Invert2x2(PetscReal * dst, PetscReal * src, const uint16_t n)
{
    PetscReal det;

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

static void Invert3x3(PetscReal * dst, PetscReal * src, const uint16_t n)
{
    PetscReal det;

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

static void Invert4x4(PetscReal * dst, PetscReal * src, const uint16_t n)
{
    PetscReal det;

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

static void Invertnxn(PetscReal * dst, PetscReal * src, const uint16_t n)
{
    memset(dst,0,n*n*sizeof(PetscReal));
    for(PetscInt i = 0; i < n; i++) dst[i*n+i] = 1.0;

    for(PetscInt i = 0; i < n; i++) {
        for(PetscInt j = 0; j < n; j++) {
            if (i != j) {
                PetscReal ratio = src[j*n+i]/src[i*n+i];
                for(PetscInt k = 0; k < n; k++) src[j*n+k] -= ratio * src[i*n+k];
                for(PetscInt k = 0; k < n; k++) dst[j*n+k] -= ratio * dst[i*n+k];
            }
        }
    }

    for(PetscInt i = 0; i < n; i++) {
        PetscReal a = src[i*n+i];
        for(PetscInt j = 0; j < n; j++) dst[i*n+j] /= a;
    }
}

void Invertmatrix(PetscReal * dst, PetscReal * src, const uint16_t n)
{
    if        (n == 1) {
        Invert1x1(dst,src,n);
    } else if (n == 2) {
        Invert2x2(dst,src,n);
    } else if (n == 3) {
        Invert3x3(dst,src,n);
    } else if (n == 4) {
        Invert4x4(dst,src,n);
    } else {
        Invertnxn(dst,src,n);
    }  
}

void MatMatMult_CIPHER(PetscReal * dst, const PetscReal * srcA, const PetscReal * srcB, const uint16_t n)
{
    PetscInt i,j,k;
    memset(dst,0,n*n*sizeof(PetscReal));
    for (k=0; k<n; k++) {
        for (j=0; j<n; j++) {
            for (i=0; i<n; i++) {
                dst[k*n+j] += srcA[k*n+i]*srcB[i*n+j];
            }
        }
    }        
}

void MatVecMult_CIPHER(PetscReal * dst, const PetscReal * srcA, const PetscReal * srcB, const uint16_t n)
{
    PetscInt j,k;
    memset(dst,0,n*sizeof(PetscReal));
    for (k=0; k<n; k++) {
        for (j=0; j<n; j++) {
            dst[k] += srcA[k*n+j]*srcB[j];
        }
    }        
}

/*
 Linear interpolation shape function
 */
static void Interpolant_linear(PetscReal *interpolant, const PetscReal *phasefrac, const uint16_t nphases)
{
    memcpy(interpolant,phasefrac,nphases*sizeof(PetscReal));
}

/*
 Linear interpolation shape function derivative
 */
static void InterpolantDerivative_linear(PetscReal *in, const PetscReal *phasefrac, const uint16_t nphases)
{}

/*
 Quadratic interpolation shape function
 */
static void Interpolant_quad(PetscReal *interpolant, const PetscReal *phasefrac, const uint16_t nphases)
{
    if (nphases == 1) {
        interpolant[0] = 1.0;
        return;
    }
    PetscReal interpolantsum = 0.0;
    for (PetscInt g = 0; g < nphases; g++) {
        interpolant[g] = phasefrac[g]*phasefrac[g];
        interpolantsum += interpolant[g];
    }
    for (PetscInt g = 0; g < nphases; g++) interpolant[g] /= interpolantsum;
}

/*
 Quadratic interpolation shape function derivative
 */
static void InterpolantDerivative_quad(PetscReal *in, const PetscReal *phasefrac, const uint16_t nphases)
{
    if (nphases == 1) {
        in[0] = 0.0;
        return;
    }
    PetscReal out[nphases];
    PetscReal avgin = 0.0, interpolantsum = 0.0;
    for (PetscInt g = 0; g < nphases; g++) {
        PetscReal interpolant = phasefrac[g]*phasefrac[g];
        interpolantsum += interpolant;
        avgin += interpolant*in[g];
    }    
    for (PetscInt g = 0; g < nphases; g++)
        out[g] = 2.0*phasefrac[g]*(in[g] - avgin/interpolantsum)/interpolantsum;
    
    memcpy(in,out,nphases*sizeof(PetscReal));
}

/*
 Cubic interpolation shape function
 */
static void Interpolant_cubic(PetscReal *interpolant, const PetscReal *phasefrac, const uint16_t nphases)
{
    if (nphases == 1) {
        interpolant[0] = 1.0;
        return;
    }
    PetscReal interpolantsum = 0.0;
    for (PetscInt g = 0; g < nphases; g++) {
        interpolant[g] = phasefrac[g]*phasefrac[g]*(3.0 - 2.0*phasefrac[g]);
        interpolantsum += interpolant[g];
    }
    for (PetscInt g = 0; g < nphases; g++) interpolant[g] /= interpolantsum;
}

/*
 Cubic interpolation shape function derivative
 */
static void InterpolantDerivative_cubic(PetscReal *in, const PetscReal *phasefrac, const uint16_t nphases)
{
    if (nphases == 1) {
        in[0] = 0.0;
        return;
    }
    PetscReal out[nphases];
    PetscReal avgin = 0.0, interpolantsum = 0.0;
    for (PetscInt g = 0; g < nphases; g++) {
        PetscReal interpolant = phasefrac[g]*phasefrac[g]*(3.0 - 2.0*phasefrac[g]);
        interpolantsum += interpolant;
        avgin += interpolant*in[g];
    }    
    for (PetscInt g = 0; g < nphases; g++)
        out[g] = 6.0*phasefrac[g]*(1.0 - phasefrac[g])*(in[g] - avgin/interpolantsum)/interpolantsum;
    
    memcpy(in,out,nphases*sizeof(PetscReal));
}

/*
 Interpolation shape function
 */
void EvalInterpolant(PetscReal *interpolant, const PetscReal *phasefrac, const uint16_t nphases)
{
    Interpolant(interpolant,phasefrac,nphases);
}

/*
 Interpolation shape function derivative
 */
void MatMulInterpolantDerivative(PetscReal *in, const PetscReal *phasefrac, const uint16_t nphases)
{
    InterpolantDerivative(in,phasefrac,nphases);
}

/*
 SimplexProjection - Project given vector on to Gibbs simplex
 */
static int Comparison(const void *x, const void *y){   
   if (*(PetscReal*)y > *(PetscReal*)x) {
       return 1;
   } else if (*(PetscReal*)y < *(PetscReal*)x) {    
       return -1;
   } else {
       return 0;
   }
}
void SimplexProjection(PetscReal *out, PetscReal *in, int size)
{
    char bget=0;
    PetscReal s[size], tmax, tempsum=0.0;
    memcpy(s,in,size*sizeof(PetscReal));
    qsort(s,size,sizeof(PetscReal),Comparison);
    for (int i=0; i<size-1; i++) {
        tempsum += s[i];
        tmax = (tempsum-1.0)/((double)(i+1));
        if (tmax >= s[i+1]) {
            bget = 1;
            break;
        }
    }
    if (!bget) tmax = (tempsum + s[size-1] - 1.0)/((PetscReal)(size));
    for (int i=0; i<size; i++) out[i] = in[i]-tmax > 0.0 ? in[i]-tmax : 0.0;
}

/*
 Initialise material module
 */
void utility_init(const AppCtx *user)
{
    if        (user->interpolation == LINEAR_INTERPOLATION) {
        Interpolant = &Interpolant_linear;
        InterpolantDerivative = &InterpolantDerivative_linear;
    } else if (user->interpolation == QUADRATIC_INTERPOLATION ) {
        Interpolant = &Interpolant_quad;
        InterpolantDerivative = &InterpolantDerivative_quad;
    } else if (user->interpolation == CUBIC_INTERPOLATION ) {
        Interpolant = &Interpolant_cubic;
        InterpolantDerivative = &InterpolantDerivative_cubic;
    }
}
