
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
    void (*Chemicalpotential  )(PetscScalar *, const PetscScalar *, const CHEMFE, const PetscInt);
    void (*Chemenergy         )(PetscScalar *, const PetscScalar *, const CHEMFE, const PetscInt);
    int  (*Composition        )(PetscScalar *, const PetscScalar *, const CHEMFE, const PetscInt);
    void (*CompositionTangent )(PetscScalar *, const PetscScalar *, const CHEMFE, const PetscInt);
    void (*CompositionMobility)(PetscScalar *, const PetscScalar *, const CHEMFE, const PetscInt);
} MATFUNC;

static MATFUNC *Matfunc;

/*
 Chemicalpotential - CALPHAD model for chemical potential
 */
static void Chemicalpotential_calphad(PetscScalar *chempot, const PetscScalar *composition, const CHEMFE energy, const PetscInt numcomps)
{
    PetscScalar  binaryfactora[numcomps][numcomps]          ,  binaryfactorb[numcomps][numcomps];
    PetscScalar ternaryfactora[numcomps][numcomps][numcomps], ternaryfactorb[numcomps][numcomps][numcomps];
    const CALPHAD *currentcalphad = &energy.calphad;
    const RK *currentbinary = &currentcalphad->binary[0];
    const RK *currentternary = &currentcalphad->ternary[0];
    for (PetscInt ck=0; ck<numcomps; ck++) {  
        for (PetscInt cj=ck+1; cj<numcomps; cj++,currentbinary++) {
            binaryfactora[ck][cj] = 0.0;
            for (PetscInt rko=0; rko<currentbinary->n; rko++) {
                binaryfactora[ck][cj] += currentbinary->enthalpy[rko]*FastPow(composition[ck] - composition[cj],rko);
            }
            binaryfactorb[ck][cj] = 0.0;
            for (PetscInt rko=1; rko<currentbinary->n; rko++) {
                binaryfactorb[ck][cj] += currentbinary->enthalpy[rko]*((PetscScalar) rko)*FastPow(composition[ck] - composition[cj],rko-1);
            }
            for (PetscInt ci=cj+1; ci<numcomps; ci++,currentternary++) {
                ternaryfactora[ck][cj][ci] = 0.0; ternaryfactorb[ck][cj][ci] = 0.0;
                for (PetscInt rko=0; rko<currentternary->n; rko++) {
                    ternaryfactora[ck][cj][ci] += currentternary->enthalpy[rko]
                                                * (  composition[currentternary->i[rko]] 
                                                   + (1.0 - composition[ck] - composition[cj] - composition[ci])/3.0);
                    ternaryfactorb[ck][cj][ci] += currentternary->enthalpy[rko];
                }
            }
        }
    }        
    currentternary = &currentcalphad->ternary[0];
    for (PetscInt ck=0; ck<numcomps; ck++) {  
        chempot[ck] += currentcalphad->unary[ck] + currentcalphad->RT*log(composition[ck]+TOL);
        for (PetscInt cj=ck+1; cj<numcomps; cj++) {
            chempot[ck] += composition[cj]*(binaryfactora[ck][cj] + composition[ck]*binaryfactorb[ck][cj]);
            chempot[cj] += composition[ck]*(binaryfactora[ck][cj] - composition[cj]*binaryfactorb[ck][cj]);
            for (PetscInt ci=cj+1; ci<numcomps; ci++,currentternary++) {
                chempot[ck] += composition[cj]*composition[ci]*(ternaryfactora[ck][cj][ci] - composition[ck]*ternaryfactorb[ck][cj][ci]/3.0);
                chempot[cj] += composition[ci]*composition[ck]*(ternaryfactora[ck][cj][ci] - composition[cj]*ternaryfactorb[ck][cj][ci]/3.0);
                chempot[ci] += composition[ck]*composition[cj]*(ternaryfactora[ck][cj][ci] - composition[ci]*ternaryfactorb[ck][cj][ci]/3.0);
                for (PetscInt rko=0; rko<currentternary->n; rko++) {
                    chempot[currentternary->i[rko]] += composition[ck]*composition[cj]*composition[ci]
                                                     * currentternary->enthalpy[rko];
                }
            }
        }
    }
}

/*
 Chemicalpotential - Quadratic model for chemical potential
 */
static void Chemicalpotential_quad(PetscScalar *chempot, const PetscScalar *composition, const CHEMFE energy, const PetscInt numcomps)
{
    const QUAD *currentquad = &energy.quad;
    for (PetscInt c=0; c<numcomps; c++)
        chempot[c] = currentquad->unary[c] + 2.0*currentquad->binary[c]*(composition[c] - currentquad->ceq[c]);
}

/*
 Chemicalpotential - None model for chemical potential
 */
static void Chemicalpotential_none(PetscScalar *chempot, const PetscScalar *composition, const CHEMFE energy, const PetscInt numcomps)
{
    memset(chempot,0,(numcomps-1)*sizeof(PetscScalar));
}

/*
 Chemicalpotential
 */
void Chemicalpotential(PetscScalar *chempot, const PetscScalar *composition, const PetscScalar *phasefrac, const uint16_t *phaseID, const AppCtx *user)
{
    memset(chempot,0,(user->nc-1)*sizeof(PetscScalar));
    for (PetscInt g=0; g<phaseID[0]; g++) {
        PetscScalar chempotk[user->nc];
        memset(chempotk,0,user->nc*sizeof(PetscScalar));
        const MATERIAL *currentmaterial = &user->material[user->phasematerialmapping[phaseID[g+1]]];
        Matfunc[user->phasematerialmapping[phaseID[g+1]]].Chemicalpotential(chempotk,&composition[g*user->nc],currentmaterial->energy,user->nc);
        for (PetscInt ck=0; ck<user->nc-1; ck++) {  
            chempot[ck] += phasefrac[g]*(chempotk[ck] - chempotk[user->nc-1]);
        }                              
    }
}

/*
 Chemenergy - CALPHAD model for chemical energy
 */
static void Chemenergy_calphad(PetscScalar *chemenergy, const PetscScalar *composition, const CHEMFE energy, const PetscInt numcomps)
{
    for (PetscInt c=0; c<numcomps; c++) assert(composition[c] + TOL > 0.0);
    const CALPHAD *currentcalphad = &energy.calphad;
    const RK *currentbinary = &currentcalphad->binary[0];
    const RK *currentternary = &currentcalphad->ternary[0];
    *chemenergy = currentcalphad->ref;
    for (PetscInt ck=0; ck<numcomps; ck++) {
        *chemenergy += currentcalphad->unary[ck]*composition[ck] + currentcalphad->RT*composition[ck]*log(composition[ck]+TOL);
        for (PetscInt cj=ck+1; cj<numcomps; cj++,currentbinary++) {
            for (PetscInt rko=0; rko<currentbinary->n; rko++)
                *chemenergy += currentbinary->enthalpy[rko]*composition[ck]*composition[cj]*FastPow(composition[ck]-composition[cj],rko);
            for (PetscInt ci=cj+1; ci<numcomps; ci++,currentternary++) {
                for (PetscInt rko=0; rko<currentternary->n; rko++)
                    *chemenergy += currentternary->enthalpy[rko]
                                 * composition[ck]*composition[cj]*composition[ci]
                                 * (composition[currentternary->i[rko]] - (1.0-composition[ck]-composition[cj]-composition[ci])/3.0);
            }
        }
    }
}

/*
 Chemenergy - quadratic model for chemical energy
 */
static void Chemenergy_quad(PetscScalar *chemenergy, const PetscScalar *composition, const CHEMFE energy, const PetscInt numcomps)
{
    const QUAD *currentquad = &energy.quad;
    *chemenergy = currentquad->ref;
    for (PetscInt c=0; c<numcomps; c++)
        *chemenergy += currentquad->unary [c]*(composition[c] - currentquad->ceq[c])
                     + currentquad->binary[c]*(composition[c] - currentquad->ceq[c])
                                             *(composition[c] - currentquad->ceq[c]);
}

/*
 Chemenergy - none model for chemical energy
 */
static void Chemenergy_none(PetscScalar *chemenergy, const PetscScalar *composition, const CHEMFE energy, const PetscInt numcomps)
{
    *chemenergy = 0.0;
}

/*
 Chemenergy - CALPHAD model for chemical energy
 */
void Chemenergy(PetscScalar *chemenergy, const PetscScalar *composition, const PetscScalar *chempot, const uint16_t *phaseID, const AppCtx *user)
{
    for (PetscInt g=0; g<phaseID[0]; g++) {
        const MATERIAL *currentmaterial = &user->material[user->phasematerialmapping[phaseID[g+1]]];
        Matfunc[user->phasematerialmapping[phaseID[g+1]]].Chemenergy(&chemenergy[g],&composition[g*user->nc],currentmaterial->energy,user->nc);
        for (PetscInt c=0; c<user->nc-1; c++) chemenergy[g] -= chempot[c]*composition[g*user->nc+c];
        chemenergy[g] /= currentmaterial->molarvolume;    
    }
}

/*
 Composition - Semi-implicit concentration per phase
 */
static int Composition_calphad(PetscScalar *composition, const PetscScalar *chempot, const CHEMFE energy, const PetscInt numcomps)
{
    PetscScalar  binaryfactora[numcomps][numcomps]          ,  binaryfactorb[numcomps][numcomps];
    PetscScalar ternaryfactora[numcomps][numcomps][numcomps], ternaryfactorb[numcomps][numcomps][numcomps];
    const CALPHAD *currentcalphad = &energy.calphad;
    const RK *currentbinary = &currentcalphad->binary[0];
    const RK *currentternary = &currentcalphad->ternary[0];
    for (PetscInt ck=0; ck<numcomps; ck++) {  
        for (PetscInt cj=ck+1; cj<numcomps; cj++,currentbinary++) {
            binaryfactora[ck][cj] = 0.0;
            for (PetscInt rko=0; rko<currentbinary->n; rko++) {
                binaryfactora[ck][cj] += currentbinary->enthalpy[rko]*FastPow(composition[ck] - composition[cj],rko);
            }
            binaryfactorb[ck][cj] = 0.0;
            for (PetscInt rko=1; rko<currentbinary->n; rko++) {
                binaryfactorb[ck][cj] += 
                    currentbinary->enthalpy[rko]*((PetscScalar) rko)*FastPow(composition[ck] - composition[cj],rko-1);
            }
            for (PetscInt ci=cj+1; ci<numcomps; ci++,currentternary++) {
                ternaryfactora[ck][cj][ci] = 0.0; ternaryfactorb[ck][cj][ci] = 0.0;
                for (PetscInt rko=0; rko<currentternary->n; rko++) {
                    ternaryfactora[ck][cj][ci] += 
                        currentternary->enthalpy[rko]*(composition[currentternary->i[rko]] + (1.0 - composition[ck] - composition[cj] - composition[ci])/3.0);
                    ternaryfactorb[ck][cj][ci] += currentternary->enthalpy[rko];
                }
            }
        }
    }        
    PetscScalar chempote[numcomps];
    memset(chempote,0,numcomps*sizeof(PetscScalar));
    currentternary = &currentcalphad->ternary[0];
    for (PetscInt ck=0; ck<numcomps; ck++) {  
        chempote[ck] += currentcalphad->unary[ck];
        for (PetscInt cj=ck+1; cj<numcomps; cj++) {
            chempote[ck] += composition[cj]*(binaryfactora[ck][cj] + composition[ck]*binaryfactorb[ck][cj]);
            chempote[cj] += composition[ck]*(binaryfactora[ck][cj] - composition[cj]*binaryfactorb[ck][cj]);
            for (PetscInt ci=cj+1; ci<numcomps; ci++,currentternary++) {
                chempote[ck] += composition[cj]*composition[ci]*(ternaryfactora[ck][cj][ci] - composition[ck]*ternaryfactorb[ck][cj][ci]/3.0);
                chempote[cj] += composition[ci]*composition[ck]*(ternaryfactora[ck][cj][ci] - composition[cj]*ternaryfactorb[ck][cj][ci]/3.0);
                chempote[ci] += composition[ck]*composition[cj]*(ternaryfactora[ck][cj][ci] - composition[ci]*ternaryfactorb[ck][cj][ci]/3.0);
                for (PetscInt rko=0; rko<currentternary->n; rko++) {
                    chempote[currentternary->i[rko]] += composition[ck]*composition[cj]*composition[ci]
                                                      * currentternary->enthalpy[rko];
                }
            }
        }
    }        
    PetscScalar sumexp = 0.0;
    for (PetscInt c=0; c<numcomps-1; c++) {
        composition[c] = exp((chempot[c] - chempote[c] + chempote[numcomps-1])/currentcalphad->RT);
        sumexp += composition[c];
    }
    composition[numcomps-1] = 1.0;
    for (PetscInt c=0; c<numcomps-1; c++) {
        composition[c] /= (1.0 + sumexp);   
        composition[numcomps-1] -= composition[c];
    } 
    int ierr = 1;
    for (PetscInt c=0; c<numcomps; c++) {
        ierr = ierr && (composition[c] + TOL > 0.0);
    }
    return ierr;    
}

/*
 Composition - Semi-implicit concentration per phase
 */
static int Composition_quad(PetscScalar *composition, const PetscScalar *chempot, const CHEMFE energy, const PetscInt numcomps)
{
    const QUAD *currentquad = &energy.quad;
    PetscScalar rhs[numcomps];
    PetscScalar sument = 0.0, sumc = 0.0;
    for (PetscInt c=0; c<numcomps-1; c++) {
        rhs[c] = 0.5*chempot[c] - 0.5*currentquad->unary[c] + 0.5*currentquad->unary[numcomps-1]
               + currentquad->binary[c]*currentquad->ceq[c] - currentquad->binary[numcomps-1]*currentquad->ceq[numcomps-1]
               + currentquad->binary[numcomps-1];
        composition[c] = rhs[c]/currentquad->binary[c];
        sumc += composition[c]; sument += 1.0/currentquad->binary[c];
    }
    sument = currentquad->binary[numcomps-1]/(1.0 + sument*currentquad->binary[numcomps-1]);            
    composition[numcomps-1] = 1.0;
    for (PetscInt c=0; c<numcomps-1; c++) {
        composition[c] -= sument*sumc/currentquad->binary[c];
        composition[numcomps-1] -= composition[c];
    }
    return 1;
}

/*
 Composition - const concentration per phase
 */
static int Composition_none(PetscScalar *composition, const PetscScalar *chempot, const CHEMFE energy, const PetscInt numcomps)
{
    composition[0] = 1.0;
    return 1;
}

/*
 Composition - Semi-implicit concentration per phase
 */
int Composition(PetscScalar *composition, const PetscScalar *chempot, const uint16_t *phaseID, const AppCtx *user)
{
    int err;
    for (PetscInt g=0; g<phaseID[0]; g++) {
        const MATERIAL *currentmaterial = &user->material[user->phasematerialmapping[phaseID[g+1]]];
        err = Matfunc[user->phasematerialmapping[phaseID[g+1]]].Composition(&composition[g*user->nc],chempot,currentmaterial->energy,user->nc);
    }
    return err;
}

/*
 Composition tangent wrt chemical potential
 */
static void CompositionTangent_calphad(PetscScalar *compositiontangent, const PetscScalar *composition, const CHEMFE energy, const PetscInt numcomps)
{
    memset(compositiontangent,0,(numcomps-1)*(numcomps-1)*sizeof(PetscScalar));
    const CALPHAD *currentcalphad = &energy.calphad;
    for (PetscInt cj=0; cj<numcomps-1; cj++) {
        compositiontangent[cj*(numcomps-1)+cj] = composition[cj]/currentcalphad->RT;
        for (PetscInt ci=0; ci<numcomps-1; ci++) {
            compositiontangent[cj*(numcomps-1)+ci] -= composition[cj]*composition[ci]/currentcalphad->RT;
        }
    }        
}

/*
 Composition tangent wrt chemical potential
 */
static void CompositionTangent_quad(PetscScalar *compositiontangent, const PetscScalar *composition, const CHEMFE energy, const PetscInt numcomps)
{
    memset(compositiontangent,0,(numcomps-1)*(numcomps-1)*sizeof(PetscScalar));
    const QUAD *currentquad = &energy.quad;
    PetscScalar sument = 0.0;
    for (PetscInt c=0; c<numcomps-1; c++) sument += 1.0/currentquad->binary[c];
    sument = currentquad->binary[numcomps-1]/(1.0 + sument*currentquad->binary[numcomps-1]);            
    for (PetscInt cj=0; cj<numcomps-1; cj++) {
        compositiontangent[cj*(numcomps-1)+cj] = 0.5/currentquad->binary[cj];
        for (PetscInt ci=0; ci<numcomps-1; ci++) {
            compositiontangent[cj*(numcomps-1)+ci] -= 0.5*sument/currentquad->binary[ci]/currentquad->binary[cj];
        }
    }        
}

/*
 Composition tangent wrt chemical potential
 */
static void CompositionTangent_none(PetscScalar *compositiontangent, const PetscScalar *composition, const CHEMFE energy, const PetscInt numcomps)
{
    memset(compositiontangent,0,(numcomps-1)*(numcomps-1)*sizeof(PetscScalar));
}

/*
 Composition tangent wrt chemical potential
 */
void CompositionTangent(PetscScalar *compositiontangent, const PetscScalar *composition, const PetscScalar *phasefrac, const uint16_t *phaseID, const AppCtx *user)
{
    memset(compositiontangent,0,(user->nc-1)*(user->nc-1)*sizeof(PetscScalar));
    for (PetscInt g=0; g<phaseID[0]; g++) {
        const MATERIAL *currentmaterial = &user->material[user->phasematerialmapping[phaseID[g+1]]];
        PetscScalar compositiontangentg[(user->nc-1)*(user->nc-1)];
        Matfunc[user->phasematerialmapping[phaseID[g+1]]].CompositionTangent(compositiontangentg,&composition[g*user->nc],currentmaterial->energy,user->nc);
        for (PetscInt cj=0; cj<user->nc-1; cj++) {
            for (PetscInt ci=0; ci<user->nc-1; ci++) {
                compositiontangent[cj*(user->nc-1) + ci] += phasefrac[g]*compositiontangentg[cj*(user->nc-1) + ci];
            }
        }
    }
}

/*
 Composition mobility
 */
static void CompositionMobility_calphad(PetscScalar *mobilityc, const PetscScalar *composition, const CHEMFE energy, const PetscInt numcomps)
{
    memset(mobilityc,0,(numcomps-1)*sizeof(PetscScalar));
    const CALPHAD *currentcalphad = &energy.calphad;
    for (PetscInt c=0; c<numcomps-1; c++) {
        mobilityc[c] = currentcalphad->mobilityc[c];
    }        
}

/*
 Composition mobility
 */
static void CompositionMobility_quad(PetscScalar *mobilityc, const PetscScalar *composition, const CHEMFE energy, const PetscInt numcomps)
{
    memset(mobilityc,0,(numcomps-1)*sizeof(PetscScalar));
    const QUAD *currentquad = &energy.quad;
    for (PetscInt c=0; c<numcomps-1; c++) {
        mobilityc[c] = currentquad->mobilityc[c];
    }        
}

/*
 Composition mobility
 */
static void CompositionMobility_none(PetscScalar *mobilityc, const PetscScalar *composition, const CHEMFE energy, const PetscInt numcomps)
{
    memset(mobilityc,0,(numcomps-1)*sizeof(PetscScalar));
}

/*
 Composition mobility
 */
void CompositionMobility(PetscScalar *mobilityc, const PetscScalar *composition, const PetscScalar *phasefrac, const uint16_t *phaseID, const AppCtx *user)
{
    memset(mobilityc,0,(user->nc-1)*sizeof(PetscScalar));
    for (PetscInt g=0; g<phaseID[0]; g++) {
        const MATERIAL *currentmaterial = &user->material[user->phasematerialmapping[phaseID[g+1]]];
        PetscScalar mobilitycg[user->nc-1];
        Matfunc[user->phasematerialmapping[phaseID[g+1]]].CompositionMobility(mobilitycg,&composition[g*user->nc],currentmaterial->energy,user->nc);
        for (PetscInt c=0; c<user->nc-1; c++) {
            mobilityc[c] += phasefrac[g]*mobilitycg[c];
        }
    }
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
            Matfunc[mat].Chemicalpotential = &Chemicalpotential_quad;
            Matfunc[mat].Chemenergy = &Chemenergy_quad;
            Matfunc[mat].Composition = &Composition_quad;
            Matfunc[mat].CompositionTangent = &CompositionTangent_quad;
            Matfunc[mat].CompositionMobility = &CompositionMobility_quad;
        } else if (currentmaterial->model == CALPHAD_CHEMENERGY  ) {
            Matfunc[mat].Chemicalpotential = &Chemicalpotential_calphad;
            Matfunc[mat].Chemenergy = &Chemenergy_calphad;
            Matfunc[mat].Composition = &Composition_calphad;
            Matfunc[mat].CompositionTangent = &CompositionTangent_calphad;
            Matfunc[mat].CompositionMobility = &CompositionMobility_calphad;
        } else if (currentmaterial->model == NONE_CHEMENERGY     ) {
            Matfunc[mat].Chemicalpotential = &Chemicalpotential_none;
            Matfunc[mat].Chemenergy = &Chemenergy_none;
            Matfunc[mat].Composition = &Composition_none;
            Matfunc[mat].CompositionTangent = &CompositionTangent_none;
            Matfunc[mat].CompositionMobility = &CompositionMobility_none;
        }
    }    
}
