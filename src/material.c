
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <petscdm.h>
#include <petscdmda.h>
#include "utility.h"
#include "typedef.h"

/* Including my own header for checking by compiler */
#define MATERIAL_IMPORT
#include "material.h"

/* constitutive functions */
typedef struct MATFUNC {
    void (*Chemenergy) (PetscReal *, const PetscReal *, const PetscReal, const PetscReal, const CHEMFE, const PetscInt);
    void (*ChemicalpotentialExplicit) (PetscReal *, const PetscReal *, const PetscReal, const CHEMFE, const PetscInt);
    void (*ChemicalpotentialExplicitTangent) (PetscReal *, const PetscReal *, const PetscReal, const CHEMFE, const PetscInt);
    void (*ChemicalpotentialImplicit) (PetscReal *, const PetscReal *, const PetscReal, const PetscReal, const CHEMFE, const PetscInt);
    void (*ChemicalpotentialImplicitTangent) (PetscReal *, const PetscReal *, const PetscReal, const PetscReal, const CHEMFE, const PetscInt);
    void (*Composition) (PetscReal *, const PetscReal *, const PetscReal, const PetscReal, const CHEMFE, const PetscInt);
    void (*CompositionTangent) (PetscReal *, const PetscReal *, const PetscReal, const PetscReal, const CHEMFE, const PetscInt);
    void (*CompositionMobility) (PetscReal *, const PetscReal *, const PetscReal, const CHEMFE, const PetscInt);
} MATFUNC;

static MATFUNC *Matfunc;

/*
 Chemenergy - CALPHAD model for chemical energy
 */
void Temperature(PetscReal *temperature, const PetscReal time, const AppCtx *user)
{
    if (time <= user->solparams.time_list[0]) {
        *temperature = user->solparams.temperature_list[0];
        return;
    }
    if (time >= user->solparams.time_list[user->solparams.n_temperatures-1]) {
        *temperature = user->solparams.temperature_list[user->solparams.n_temperatures-1];
        return;
    }
    for (PetscInt interval=0; interval<user->solparams.n_temperatures-1; interval++) {
        if (user->solparams.time_list[interval] <= time && user->solparams.time_list[interval+1] > time) {
            *temperature = user->solparams.temperature_list[interval] 
                         + (user->solparams.temperature_list[interval+1] - user->solparams.temperature_list[interval])
                         * (time                                         - user->solparams.time_list       [interval])
                         / (user->solparams.time_list       [interval+1] - user->solparams.time_list       [interval]);
        }
    }
}

/*
 Chemenergy - CALPHAD model for chemical energy
 */
static void Chemenergy_calphad(PetscReal *chemenergy, const PetscReal *composition, const PetscReal temperature, const PetscReal temperature0, const CHEMFE energy, const PetscInt numcomps)
{
    const CALPHAD *currentcalphad = &energy.calphad;
    const RK *currentbinary = &currentcalphad->binary[0];
    const RK *currentternary = &currentcalphad->ternary[0];
    *chemenergy = currentcalphad->ref + currentcalphad->ref_T*(temperature-temperature0);
    for (PetscInt ck=0; ck<numcomps; ck++) {
        *chemenergy += (currentcalphad->unary[ck]+currentcalphad->unary_T[ck]*(temperature-temperature0))*composition[ck] 
                     + (currentcalphad->R*temperature)*composition[ck]*log(composition[ck]);
        for (PetscInt cj=ck+1; cj<numcomps; cj++,currentbinary++) {
            for (PetscInt rko=0; rko<currentbinary->n; rko++)
                *chemenergy += (currentbinary->enthalpy[rko]+currentbinary->enthalpy_T[rko]*(temperature-temperature0))*composition[ck]*composition[cj]*FastPow(composition[ck]-composition[cj],rko);
            for (PetscInt ci=cj+1; ci<numcomps; ci++,currentternary++) {
                for (PetscInt rko=0; rko<currentternary->n; rko++)
                    *chemenergy += (currentternary->enthalpy[rko] + currentternary->enthalpy_T[rko]*(temperature-temperature0))
                                 * composition[ck]*composition[cj]*composition[ci]
                                 * (composition[currentternary->i[rko]] - (1.0-composition[ck]-composition[cj]-composition[ci])/3.0);
            }
        }
    }
}

/*
 Chemenergy - quadratic model for chemical energy
 */
static void Chemenergy_quad(PetscReal *chemenergy, const PetscReal *composition, const PetscReal temperature, const PetscReal temperature0, const CHEMFE energy, const PetscInt numcomps)
{
    PetscReal ceq;
    const QUAD *currentquad = &energy.quad;
    *chemenergy = currentquad->ref + currentquad->ref_T*(temperature-temperature0);
    for (PetscInt c=0; c<numcomps; c++) {
        ceq = currentquad->ceq[c] + currentquad->ceq_T[c]*(temperature-temperature0);
        *chemenergy += (currentquad->unary [c] + currentquad->unary_T [c]*(temperature-temperature0))*(composition[c] - ceq)
                     + (currentquad->binary[c] + currentquad->binary_T[c]*(temperature-temperature0))*(composition[c] - ceq)
                                                                                                     *(composition[c] - ceq);
    }
}

/*
 Chemenergy - none model for chemical energy
 */
static void Chemenergy_none(PetscReal *chemenergy, const PetscReal *composition, const PetscReal temperature, const PetscReal temperature0, const CHEMFE energy, const PetscInt numcomps)
{
    *chemenergy = 0.0;
}

/*
 Chemenergy - CALPHAD model for chemical energy
 */
void Chemenergy(PetscReal *chemenergy, const PetscReal *composition, const PetscReal *chempot, const PetscReal temperature, const uint16_t *phaseID, const AppCtx *user)
{
    for (PetscInt g=0; g<phaseID[0]; g++) {
        const MATERIAL *currentmaterial = &user->material[user->phasematerialmapping[phaseID[g+1]]];
        Matfunc[user->phasematerialmapping[phaseID[g+1]]].Chemenergy(&chemenergy[g],&composition[g*user->ncp],temperature,currentmaterial->temperature0,currentmaterial->energy,user->ncp);
        for (PetscInt c=0; c<user->ndp; c++) chemenergy[g] -= chempot[c]*composition[g*user->ncp+c];
        chemenergy[g] /= currentmaterial->molarvolume;    
    }
}

/*
 Chemicalpotential - CALPHAD model for chemical potential
 */
static void ChemicalpotentialExplicit_calphad(PetscReal *chempot, const PetscReal *composition, const PetscReal deltaT, const CHEMFE energy, const PetscInt numcomps)
{
    memset(chempot,0,numcomps*sizeof(PetscReal));
    PetscReal  binaryfactora[numcomps][numcomps]          ,  binaryfactorb[numcomps][numcomps];
    PetscReal ternaryfactora[numcomps][numcomps][numcomps], ternaryfactorb[numcomps][numcomps][numcomps];
    const CALPHAD *currentcalphad = &energy.calphad;
    const RK *currentbinary = &currentcalphad->binary[0];
    const RK *currentternary = &currentcalphad->ternary[0];
    for (PetscInt ck=0; ck<numcomps; ck++) {  
        for (PetscInt cj=ck+1; cj<numcomps; cj++,currentbinary++) {
            binaryfactora[ck][cj] = 0.0;
            for (PetscInt rko=0; rko<currentbinary->n; rko++) {
                binaryfactora[ck][cj] += (currentbinary->enthalpy[rko]+currentbinary->enthalpy_T[rko]*deltaT)*FastPow(composition[ck] - composition[cj],rko);
            }
            binaryfactorb[ck][cj] = 0.0;
            for (PetscInt rko=1; rko<currentbinary->n; rko++) {
                binaryfactorb[ck][cj] += ((PetscReal) (rko))*(currentbinary->enthalpy[rko]+currentbinary->enthalpy_T[rko]*deltaT)*((PetscReal) rko)*FastPow(composition[ck] - composition[cj],rko-1);
            }
            for (PetscInt ci=cj+1; ci<numcomps; ci++,currentternary++) {
                ternaryfactora[ck][cj][ci] = 0.0; ternaryfactorb[ck][cj][ci] = 0.0;
                for (PetscInt rko=0; rko<currentternary->n; rko++) {
                    ternaryfactora[ck][cj][ci] += (currentternary->enthalpy[rko] + currentternary->enthalpy_T[rko]*deltaT)
                                                * (  composition[currentternary->i[rko]] 
                                                   + (1.0 - composition[ck] - composition[cj] - composition[ci])/3.0);
                    ternaryfactorb[ck][cj][ci] += (currentternary->enthalpy[rko] + currentternary->enthalpy_T[rko]*deltaT);
                }
            }
        }
    }
    currentternary = &currentcalphad->ternary[0];
    for (PetscInt ck=0; ck<numcomps; ck++) {  
        for (PetscInt cj=ck+1; cj<numcomps; cj++) {
            chempot[ck] += composition[cj]*(binaryfactora[ck][cj] + composition[ck]*binaryfactorb[ck][cj]);
            chempot[cj] += composition[ck]*(binaryfactora[ck][cj] - composition[cj]*binaryfactorb[ck][cj]);
            for (PetscInt ci=cj+1; ci<numcomps; ci++,currentternary++) {
                chempot[ck] += composition[cj]*composition[ci]*(ternaryfactora[ck][cj][ci] - composition[ck]*ternaryfactorb[ck][cj][ci]/3.0);
                chempot[cj] += composition[ci]*composition[ck]*(ternaryfactora[ck][cj][ci] - composition[cj]*ternaryfactorb[ck][cj][ci]/3.0);
                chempot[ci] += composition[ck]*composition[cj]*(ternaryfactora[ck][cj][ci] - composition[ci]*ternaryfactorb[ck][cj][ci]/3.0);
                for (PetscInt rko=0; rko<currentternary->n; rko++) {
                    chempot[currentternary->i[rko]] += composition[ck]*composition[cj]*composition[ci]
                                                     * (currentternary->enthalpy[rko] + currentternary->enthalpy_T[rko]*deltaT);
                }
            }
        }
    }
}

/*
 Chemicalpotential - Quadratic model for chemical potential
 */
static void ChemicalpotentialExplicit_quad(PetscReal *chempot, const PetscReal *composition, const PetscReal deltaT, const CHEMFE energy, const PetscInt numcomps)
{
    memset(chempot,0,numcomps*sizeof(PetscReal));
}

/*
 Chemicalpotential - None model for chemical potential
 */
static void ChemicalpotentialExplicit_none(PetscReal *chempot, const PetscReal *composition, const PetscReal deltaT, const CHEMFE energy, const PetscInt numcomps)
{
    memset(chempot,0,numcomps*sizeof(PetscReal));
}

/*
 Chemicalpotential
 */
void ChemicalpotentialExplicit(PetscReal *chempot, const PetscReal *composition, const PetscReal temperature, const uint16_t phaseID, const AppCtx *user)
{
    PetscReal chempotk[user->ncp];
    const MATERIAL *currentmaterial = &user->material[user->phasematerialmapping[phaseID]];
    Matfunc[user->phasematerialmapping[phaseID]].ChemicalpotentialExplicit(chempotk,composition,temperature-currentmaterial->temperature0,currentmaterial->energy,user->ncp);
    memset(chempot,0,(user->ndp)*sizeof(PetscReal));
    for (PetscInt ck=0; ck<user->ndp; ck++) {  
        chempot[ck] = (chempotk[ck] - chempotk[user->ndp]);
    }                              
}

/*
 Chemicalpotential - CALPHAD model for chemical potential
 */
static void ChemicalpotentialExplicitTangent_calphad(PetscReal *chempottangent, const PetscReal *composition, const PetscReal deltaT, const CHEMFE energy, const PetscInt numcomps)
{
    PetscReal chempottangent_f[numcomps*numcomps];
    PetscReal val;
    PetscReal binaryfactora[numcomps][numcomps],  binaryfactorb[numcomps][numcomps],  binaryfactorc[numcomps][numcomps];
    PetscReal ternaryfactora[numcomps][numcomps][numcomps], ternaryfactorb[numcomps][numcomps][numcomps];
    const CALPHAD *currentcalphad = &energy.calphad;
    const RK *currentbinary = &currentcalphad->binary[0];
    const RK *currentternary = &currentcalphad->ternary[0];
    for (PetscInt ck=0; ck<numcomps; ck++) {  
        for (PetscInt cj=ck+1; cj<numcomps; cj++,currentbinary++) {
            binaryfactora[ck][cj] = 0.0;
            for (PetscInt rko=0; rko<currentbinary->n; rko++) {
                binaryfactora[ck][cj] += (currentbinary->enthalpy[rko]+currentbinary->enthalpy_T[rko]*deltaT)*FastPow(composition[ck] - composition[cj],rko);
            }
            binaryfactorb[ck][cj] = 0.0;
            for (PetscInt rko=1; rko<currentbinary->n; rko++) {
                binaryfactorb[ck][cj] += ((PetscReal) (rko))*(currentbinary->enthalpy[rko]+currentbinary->enthalpy_T[rko]*deltaT)*((PetscReal) rko)*FastPow(composition[ck] - composition[cj],rko-1);
            }
            binaryfactorc[ck][cj] = 0.0;
            for (PetscInt rko=2; rko<currentbinary->n; rko++) {
                binaryfactorc[ck][cj] += ((PetscReal) (rko))*((PetscReal) (rko - 1))*(currentbinary->enthalpy[rko]+currentbinary->enthalpy_T[rko]*deltaT)*((PetscReal) rko)*FastPow(composition[ck] - composition[cj],rko-2);
            }
            for (PetscInt ci=cj+1; ci<numcomps; ci++,currentternary++) {
                ternaryfactora[ck][cj][ci] = 0.0; ternaryfactorb[ck][cj][ci] = 0.0;
                for (PetscInt rko=0; rko<currentternary->n; rko++) {
                    ternaryfactora[ck][cj][ci] += (currentternary->enthalpy[rko]+currentternary->enthalpy_T[rko]*deltaT)
                                                * (  composition[currentternary->i[rko]] 
                                                   + (1.0 - composition[ck] - composition[cj] - composition[ci])/3.0);
                    ternaryfactorb[ck][cj][ci] += (currentternary->enthalpy[rko]+currentternary->enthalpy_T[rko]*deltaT);
                }
            }
        }
    }        
    memset(chempottangent_f,0,numcomps*numcomps*sizeof(PetscReal));
    currentternary = &currentcalphad->ternary[0];
    for (PetscInt ck=0; ck<numcomps; ck++) {  
        for (PetscInt cj=ck+1; cj<numcomps; cj++) {
            val = binaryfactora[ck][cj] + (composition[ck] - composition[cj])*binaryfactorb[ck][cj] - composition[ck]*composition[cj]*binaryfactorc[ck][cj];
            chempottangent_f[ck*numcomps+cj] += val;
            chempottangent_f[cj*numcomps+ck] += val;
            for (PetscInt ci=cj+1; ci<numcomps; ci++,currentternary++) {
                val = composition[ci]*ternaryfactora[ck][cj][ci] - (composition[ci]*composition[ck] + composition[ci]*composition[cj])*ternaryfactorb[ck][cj][ci];
                chempottangent_f[ck*numcomps+cj] += val;
                chempottangent_f[cj*numcomps+ck] += val;
                val = composition[ck]*ternaryfactora[ck][cj][ci] - (composition[ck]*composition[cj] + composition[ck]*composition[ci])*ternaryfactorb[ck][cj][ci];
                chempottangent_f[cj*numcomps+ci] += val;
                chempottangent_f[ci*numcomps+cj] += val;
                val = composition[cj]*ternaryfactora[ck][cj][ci] - (composition[cj]*composition[ck] + composition[cj]*composition[ci])*ternaryfactorb[ck][cj][ci];
                chempottangent_f[ck*numcomps+ci] += val;
                chempottangent_f[ci*numcomps+ck] += val;
                for (PetscInt rko=0; rko<currentternary->n; rko++) {
                    val = composition[ck]*composition[cj]*(currentternary->enthalpy[rko]+currentternary->enthalpy_T[rko]*deltaT);
                    chempottangent_f[currentternary->i[rko]*numcomps+ci] += val;
                    chempottangent_f[ci*numcomps+currentternary->i[rko]] += val;
                    val = composition[cj]*composition[ci]*(currentternary->enthalpy[rko]+currentternary->enthalpy_T[rko]*deltaT);
                    chempottangent_f[currentternary->i[rko]*numcomps+ck] += val;
                    chempottangent_f[ck*numcomps+currentternary->i[rko]] += val;
                    val = composition[ci]*composition[ck]*(currentternary->enthalpy[rko]+currentternary->enthalpy_T[rko]*deltaT);
                    chempottangent_f[currentternary->i[rko]*numcomps+cj] += val;
                    chempottangent_f[cj*numcomps+currentternary->i[rko]] += val;
                }
            }
        }
    }
    memset(chempottangent,0,(numcomps-1)*(numcomps-1)*sizeof(PetscReal));
    for (PetscInt ck=0; ck<numcomps-1; ck++) {  
        for (PetscInt cj=0; cj<numcomps-1; cj++) {
            chempottangent[ck*(numcomps-1)+cj] = chempottangent_f[ ck         *numcomps+ cj         ]
                                               - chempottangent_f[ ck         *numcomps+(numcomps-1)]
                                               - chempottangent_f[(numcomps-1)*numcomps+ cj         ]
                                               + chempottangent_f[(numcomps-1)*numcomps+(numcomps-1)];
        }
    }
}

/*
 Chemicalpotential - Quadratic model for chemical potential
 */
static void ChemicalpotentialExplicitTangent_quad(PetscReal *chempottangent, const PetscReal *composition, const PetscReal deltaT, const CHEMFE energy, const PetscInt numcomps)
{
    memset(chempottangent,0,(numcomps-1)*(numcomps-1)*sizeof(PetscReal));
}

/*
 Chemicalpotential - None model for chemical potential
 */
static void ChemicalpotentialExplicitTangent_none(PetscReal *chempottangent, const PetscReal *composition, const PetscReal deltaT, const CHEMFE energy, const PetscInt numcomps)
{
    memset(chempottangent,0,(numcomps-1)*(numcomps-1)*sizeof(PetscReal));
}

/*
 Composition tangent wrt chemical potential
 */
void ChemicalpotentialExplicitTangent(PetscReal *chempottangent, const PetscReal *composition, const PetscReal temperature, const uint16_t phaseID, const AppCtx *user)
{
    const MATERIAL *currentmaterial = &user->material[user->phasematerialmapping[phaseID]];
    Matfunc[user->phasematerialmapping[phaseID]].ChemicalpotentialExplicitTangent(chempottangent,composition,temperature-currentmaterial->temperature0, currentmaterial->energy,user->ncp);
}

/*
 Chemicalpotential - CALPHAD model for chemical potential
 */
static void ChemicalpotentialImplicit_calphad(PetscReal *chempot, const PetscReal *composition, const PetscReal temperature, const PetscReal temperature0, const CHEMFE energy, const PetscInt numcomps)
{
    memset(chempot,0,numcomps*sizeof(PetscReal));
    const CALPHAD *currentcalphad = &energy.calphad;
    for (PetscInt ck=0; ck<numcomps; ck++) {  
        chempot[ck] = currentcalphad->unary[ck] + currentcalphad->unary_T[ck]*(temperature - temperature0) 
                    + currentcalphad->R*temperature*log(composition[ck]);
    }
}

/*
 Chemicalpotential - Quadratic model for chemical potential
 */
static void ChemicalpotentialImplicit_quad(PetscReal *chempot, const PetscReal *composition, const PetscReal temperature, const PetscReal temperature0, const CHEMFE energy, const PetscInt numcomps)
{
    memset(chempot,0,numcomps*sizeof(PetscReal));
    const QUAD *currentquad = &energy.quad;
    for (PetscInt c=0; c<numcomps-1; c++)
        chempot[c] =     (currentquad->unary [c] + currentquad->unary_T [c]*(temperature - temperature0)) 
                   + 2.0*(currentquad->binary[c] + currentquad->binary_T[c]*(temperature - temperature0))*(composition[c] - (currentquad->ceq[c]+currentquad->ceq_T[c]*(temperature-temperature0)));
}

/*
 Chemicalpotential - None model for chemical potential
 */
static void ChemicalpotentialImplicit_none(PetscReal *chempot, const PetscReal *composition, const PetscReal temperature, const PetscReal temperature0, const CHEMFE energy, const PetscInt numcomps)
{
    memset(chempot,0,numcomps*sizeof(PetscReal));
}

/*
 Chemicalpotential
 */
void ChemicalpotentialImplicit(PetscReal *chempot, const PetscReal *composition, const PetscReal temperature, const uint16_t phaseID, const AppCtx *user)
{
    PetscReal chempotk[user->ncp];
    const MATERIAL *currentmaterial = &user->material[user->phasematerialmapping[phaseID]];
    Matfunc[user->phasematerialmapping[phaseID]].ChemicalpotentialImplicit(chempotk,composition, temperature, currentmaterial->temperature0, currentmaterial->energy,user->ncp);
    memset(chempot,0,(user->ndp)*sizeof(PetscReal));
    for (PetscInt ck=0; ck<user->ndp; ck++) {  
        chempot[ck] = (chempotk[ck] - chempotk[user->ndp]);
    }                              
}

/*
 Chemicalpotential - CALPHAD model for chemical potential
 */
static void ChemicalpotentialImplicitTangent_calphad(PetscReal *chempottangent, const PetscReal *composition, const PetscReal temperature, const PetscReal temperature0, const CHEMFE energy, const PetscInt numcomps)
{
    memset(chempottangent,0,(numcomps-1)*(numcomps-1)*sizeof(PetscReal));
    const CALPHAD *currentcalphad = &energy.calphad;
    for (PetscInt ck=0; ck<numcomps-1; ck++) {  
        chempottangent[ck*(numcomps-1)+ck] = (currentcalphad->R*temperature)/composition[ck];
        for (PetscInt cj=0; cj<numcomps-1; cj++) {  
            chempottangent[ck*(numcomps-1)+cj] += (currentcalphad->R*temperature)/composition[numcomps-1];
        }    
    }
}

/*
 Chemicalpotential - Quadratic model for chemical potential
 */
static void ChemicalpotentialImplicitTangent_quad(PetscReal *chempottangent, const PetscReal *composition, const PetscReal temperature, const PetscReal temperature0, const CHEMFE energy, const PetscInt numcomps)
{
    memset(chempottangent,0,(numcomps-1)*(numcomps-1)*sizeof(PetscReal));
    const QUAD *currentquad = &energy.quad;
    for (PetscInt ck=0; ck<numcomps-1; ck++) {  
        chempottangent[ck*(numcomps-1)+ck] = 2.0*(currentquad->binary[ck]+currentquad->binary_T[ck]*(temperature - temperature0));
        for (PetscInt cj=0; cj<numcomps-1; cj++) {
            chempottangent[ck*(numcomps-1)+cj] += 2.0*(currentquad->binary[numcomps-1]+currentquad->binary_T[numcomps-1]*(temperature - temperature0));
        }
    }
}

/*
 Chemicalpotential - None model for chemical potential
 */
static void ChemicalpotentialImplicitTangent_none(PetscReal *chempottangent, const PetscReal *composition, const PetscReal temperature, const PetscReal temperature0, const CHEMFE energy, const PetscInt numcomps)
{
    memset(chempottangent,0,(numcomps-1)*(numcomps-1)*sizeof(PetscReal));
}

/*
 Composition tangent wrt chemical potential
 */
void ChemicalpotentialImplicitTangent(PetscReal *chempottangent, const PetscReal *composition, const PetscReal temperature, const uint16_t phaseID, const AppCtx *user)
{
    const MATERIAL *currentmaterial = &user->material[user->phasematerialmapping[phaseID]];
    Matfunc[user->phasematerialmapping[phaseID]].ChemicalpotentialImplicitTangent(chempottangent,composition, temperature, currentmaterial->temperature0, currentmaterial->energy,user->ncp);
}

/*
 Chemicalpotential
 */
void Chemicalpotential(PetscReal *chempot, const PetscReal *composition, const PetscReal temperature, const uint16_t phaseID, const AppCtx *user)
{
    PetscReal chempotim[user->ndp], chempotex[user->ndp];
    ChemicalpotentialImplicit(chempotim,composition, temperature, phaseID,user);
    ChemicalpotentialExplicit(chempotex,composition, temperature, phaseID,user);
    for (PetscInt ck=0; ck<user->ndp; ck++) {  
        chempot[ck] = (chempotim[ck] + chempotex[ck]);
    }                              
}

/*
 Composition tangent wrt chemical potential
 */
void ChemicalpotentialTangent(PetscReal *chempottangent, const PetscReal *composition, const PetscReal temperature, const uint16_t phaseID, const AppCtx *user)
{
    PetscReal chempotim_tangent[user->ndp*user->ndp], chempotex_tangent[user->ndp*user->ndp];
    ChemicalpotentialImplicitTangent(chempotim_tangent,composition, temperature, phaseID,user);
    ChemicalpotentialExplicitTangent(chempotex_tangent,composition, temperature, phaseID,user);
    for (PetscInt ck=0; ck<user->ndp*user->ndp; ck++) {  
        chempottangent[ck] = (chempotim_tangent[ck] + chempotex_tangent[ck]);
    }                              
}

/*
 Composition - Semi-implicit concentration per phase
 */
static void Composition_calphad(PetscReal *composition, const PetscReal *chempot_im, const PetscReal temperature, const PetscReal temperature0, const CHEMFE energy, const PetscInt numcomps)
{
    const CALPHAD *currentcalphad = &energy.calphad;
    PetscReal sumexp = 0.0;
    for (PetscInt c=0; c<numcomps-1; c++) {
        composition[c] = exp((   
                                chempot_im[c] 
                              - (currentcalphad->unary[c]         +currentcalphad->unary_T[c]         *(temperature - temperature0)) 
                              + (currentcalphad->unary[numcomps-1]+currentcalphad->unary_T[numcomps-1]*(temperature - temperature0))
                             )/(currentcalphad->R*temperature));
        sumexp += composition[c];
    }
    composition[numcomps-1] = 1.0;
    for (PetscInt c=0; c<numcomps-1; c++) {
        composition[c] /= (1.0 + sumexp);   
        composition[numcomps-1] -= composition[c];
    } 
}

/*
 Composition - Semi-implicit concentration per phase
 */
static void Composition_quad(PetscReal *composition, const PetscReal *chempot_im, const PetscReal temperature, const PetscReal temperature0, const CHEMFE energy, const PetscInt numcomps)
{
    const QUAD *currentquad = &energy.quad;
    composition[numcomps-1] = 1.0;
    for (PetscInt c=0; c<numcomps-1; c++) {
        composition[c] = (currentquad->ceq[c]+currentquad->ceq_T[c]*(temperature-temperature0))
                       + 0.5*(chempot_im[c] - (currentquad->unary[c]+currentquad->unary_T[c]*(temperature - temperature0)))/(currentquad->binary[c]+currentquad->binary_T[c]*(temperature - temperature0));
        composition[numcomps-1] -= composition[c];
    }
}

/*
 Composition - const concentration per phase
 */
static void Composition_none(PetscReal *composition, const PetscReal *chempot_im, const PetscReal temperature, const PetscReal temperature0, const CHEMFE energy, const PetscInt numcomps)
{
    memset(composition,0,numcomps*sizeof(PetscReal));
    composition[numcomps-1] = 1.0;
}

/*
 Composition - Semi-implicit concentration per phase
 */
void Composition(PetscReal *composition, const PetscReal *chempot_im, const PetscReal temperature, const uint16_t phaseID, const AppCtx *user)
{
    const MATERIAL *currentmaterial = &user->material[user->phasematerialmapping[phaseID]];
    Matfunc[user->phasematerialmapping[phaseID]].Composition(composition,chempot_im, temperature, currentmaterial->temperature0, currentmaterial->energy,user->ncp);
}

/*
 Composition tangent wrt chemical potential
 */
static void CompositionTangent_calphad(PetscReal *compositiontangent, const PetscReal *composition, const PetscReal temperature, const PetscReal temperature0, const CHEMFE energy, const PetscInt numcomps)
{
    memset(compositiontangent,0,(numcomps-1)*(numcomps-1)*sizeof(PetscReal));
    const CALPHAD *currentcalphad = &energy.calphad;
    for (PetscInt cj=0; cj<numcomps-1; cj++) {
        compositiontangent[cj*(numcomps-1)+cj] = composition[cj]/(currentcalphad->R*temperature);
        for (PetscInt ci=0; ci<numcomps-1; ci++) {
            compositiontangent[cj*(numcomps-1)+ci] -= composition[cj]*composition[ci]/(currentcalphad->R*temperature);
        }
    }        
}

/*
 Composition tangent wrt chemical potential
 */
static void CompositionTangent_quad(PetscReal *compositiontangent, const PetscReal *composition, const PetscReal temperature, const PetscReal temperature0, const CHEMFE energy, const PetscInt numcomps)
{
    memset(compositiontangent,0,(numcomps-1)*(numcomps-1)*sizeof(PetscReal));
    const QUAD *currentquad = &energy.quad;
    for (PetscInt cj=0; cj<numcomps-1; cj++) {
        compositiontangent[cj*(numcomps-1)+cj] = 0.5/(currentquad->binary[cj]-currentquad->binary_T[cj]*(temperature - temperature0));
    }        
}

/*
 Composition tangent wrt chemical potential
 */
static void CompositionTangent_none(PetscReal *compositiontangent, const PetscReal *composition, const PetscReal temperature, const PetscReal temperature0, const CHEMFE energy, const PetscInt numcomps)
{
    memset(compositiontangent,0,(numcomps-1)*(numcomps-1)*sizeof(PetscReal));
}

/*
 Composition tangent wrt chemical potential
 */
void CompositionTangent(PetscReal *compositiontangent, const PetscReal *composition, const PetscReal temperature, const uint16_t phaseID, const AppCtx *user)
{
    memset(compositiontangent,0,(user->ndp)*(user->ndp)*sizeof(PetscReal));
    const MATERIAL *currentmaterial = &user->material[user->phasematerialmapping[phaseID]];
    Matfunc[user->phasematerialmapping[phaseID]].CompositionTangent(compositiontangent,composition,temperature, currentmaterial->temperature0, currentmaterial->energy,user->ncp);
}

/*
 Composition mobility
 */
static void CompositionMobility_calphad(PetscReal *mobilityc, const PetscReal *composition, const PetscReal deltaT, const CHEMFE energy, const PetscInt numcomps)
{
    memset(mobilityc,0,(numcomps-1)*(numcomps-1)*sizeof(PetscReal));
    const CALPHAD *currentcalphad = &energy.calphad;
    PetscScalar summobility = 0.0, val;
    for (PetscInt ck=0; ck<numcomps; ck++) {
        summobility += composition[ck]*(currentcalphad->mobilityc[ck]+currentcalphad->mobilityc_T[ck]*deltaT);
    }        
    for (PetscInt ck=0; ck<numcomps-1; ck++) {
        mobilityc[ck*(numcomps-1)+ck] = composition[ck]*((currentcalphad->mobilityc[ck]+currentcalphad->mobilityc_T[ck]*deltaT) + composition[ck]*(summobility - 2.0*(currentcalphad->mobilityc[ck]+currentcalphad->mobilityc_T[ck]*deltaT)));
        for (PetscInt cj=ck+1; cj<numcomps-1; cj++) {
            val = composition[ck]*composition[cj]*(summobility - (currentcalphad->mobilityc[cj]+currentcalphad->mobilityc_T[cj]*deltaT) - (currentcalphad->mobilityc[ck] + currentcalphad->mobilityc_T[ck]*deltaT));
            mobilityc[ck*(numcomps-1)+cj] = val; 
            mobilityc[cj*(numcomps-1)+ck] = val; 
        }        
    }        
}

/*
 Composition mobility
 */
static void CompositionMobility_quad(PetscReal *mobilityc, const PetscReal *composition, const PetscReal deltaT, const CHEMFE energy, const PetscInt numcomps)
{
    memset(mobilityc,0,(numcomps-1)*(numcomps-1)*sizeof(PetscReal));
    const QUAD *currentquad = &energy.quad;
    PetscScalar summobility = 0.0, val;
    for (PetscInt ck=0; ck<numcomps; ck++) {
        summobility += composition[ck]*(currentquad->mobilityc[ck]+currentquad->mobilityc_T[ck]*deltaT);
    }        
    for (PetscInt ck=0; ck<numcomps-1; ck++) {
        mobilityc[ck*(numcomps-1)+ck] = composition[ck]*((currentquad->mobilityc[ck]+currentquad->mobilityc_T[ck]*deltaT) + composition[ck]*(summobility - 2.0*(currentquad->mobilityc[ck]+currentquad->mobilityc_T[ck]*deltaT)));
        for (PetscInt cj=ck+1; cj<numcomps-1; cj++) {
            val = composition[ck]*composition[cj]*(summobility - (currentquad->mobilityc[cj]+currentquad->mobilityc_T[cj]*deltaT) - (currentquad->mobilityc[ck]+currentquad->mobilityc_T[ck]*deltaT));
            mobilityc[ck*(numcomps-1)+cj] = val; 
            mobilityc[cj*(numcomps-1)+ck] = val; 
        }        
    }        
}

/*
 Composition mobility
 */
static void CompositionMobility_none(PetscReal *mobilityc, const PetscReal *composition, const PetscReal deltaT, const CHEMFE energy, const PetscInt numcomps)
{
    memset(mobilityc,0,(numcomps-1)*(numcomps-1)*sizeof(PetscReal));
}

/*
 Composition mobility
 */
void CompositionMobility(PetscReal *mobilityc, const PetscReal *composition, const PetscReal temperature, const uint16_t phaseID, const AppCtx *user)
{
    memset(mobilityc,0,(user->ndp)*(user->ndp)*sizeof(PetscReal));
    const MATERIAL *currentmaterial = &user->material[user->phasematerialmapping[phaseID]];
    Matfunc[user->phasematerialmapping[phaseID]].CompositionMobility(mobilityc,composition, temperature-currentmaterial->temperature0, currentmaterial->energy,user->ncp);
}

/*
 Initialise material module
 */
void material_init(const AppCtx *user)
{
    Matfunc = (MATFUNC *) malloc(user->nmat*sizeof(struct MATFUNC));
    MATERIAL *currentmaterial = &user->material[0];
    for (PetscInt mat=0; mat<user->nmat; mat++, currentmaterial++) {
        if        (currentmaterial->model == QUADRATIC_CHEMENERGY) {
            Matfunc[mat].Chemenergy = &Chemenergy_quad;
            Matfunc[mat].ChemicalpotentialExplicit = &ChemicalpotentialExplicit_quad;
            Matfunc[mat].ChemicalpotentialExplicitTangent = &ChemicalpotentialExplicitTangent_quad;
            Matfunc[mat].ChemicalpotentialImplicit = &ChemicalpotentialImplicit_quad;
            Matfunc[mat].ChemicalpotentialImplicitTangent = &ChemicalpotentialImplicitTangent_quad;
            Matfunc[mat].Composition = &Composition_quad;
            Matfunc[mat].CompositionTangent = &CompositionTangent_quad;
            Matfunc[mat].CompositionMobility = &CompositionMobility_quad;
        } else if (currentmaterial->model == CALPHAD_CHEMENERGY  ) {
            Matfunc[mat].Chemenergy = &Chemenergy_calphad;
            Matfunc[mat].ChemicalpotentialExplicit = &ChemicalpotentialExplicit_calphad;
            Matfunc[mat].ChemicalpotentialExplicitTangent = &ChemicalpotentialExplicitTangent_calphad;
            Matfunc[mat].ChemicalpotentialImplicit = &ChemicalpotentialImplicit_calphad;
            Matfunc[mat].ChemicalpotentialImplicitTangent = &ChemicalpotentialImplicitTangent_calphad;
            Matfunc[mat].Composition = &Composition_calphad;
            Matfunc[mat].CompositionTangent = &CompositionTangent_calphad;
            Matfunc[mat].CompositionMobility = &CompositionMobility_calphad;
        } else if (currentmaterial->model == NONE_CHEMENERGY     ) {
            Matfunc[mat].Chemenergy = &Chemenergy_none;
            Matfunc[mat].ChemicalpotentialExplicit = &ChemicalpotentialExplicit_none;
            Matfunc[mat].ChemicalpotentialExplicitTangent = &ChemicalpotentialExplicitTangent_none;
            Matfunc[mat].ChemicalpotentialImplicit = &ChemicalpotentialImplicit_none;
            Matfunc[mat].ChemicalpotentialImplicitTangent = &ChemicalpotentialImplicitTangent_none;
            Matfunc[mat].Composition = &Composition_none;
            Matfunc[mat].CompositionTangent = &CompositionTangent_none;
            Matfunc[mat].CompositionMobility = &CompositionMobility_none;
        }
    }    
}
