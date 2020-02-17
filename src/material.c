
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

#define ANGSTROM 1.0e-10
#define KBANGST2 1.38064852E-3

/* constitutive functions */
typedef struct MATFUNC {
    void (*Chemenergy)            (PetscReal *, const PetscReal *, const PetscReal, const CHEMFE, const PetscInt);
    void (*SitepotentialExplicit) (PetscReal *, const PetscReal *, const PetscReal, const CHEMFE, const PetscInt);
    void (*SitepotentialImplicit) (PetscReal *, const PetscReal *, const PetscReal, const CHEMFE, const PetscInt);
    void (*Sitefrac)              (PetscReal *, const PetscReal *, const PetscReal, const CHEMFE, const PetscInt);
    void (*SitefracTangent)       (PetscReal *, const PetscReal *, const PetscReal, const CHEMFE, const PetscInt);
    void (*CompositionMobility)   (PetscReal *, const PetscReal *, const PetscReal, const CHEMFE, const PetscInt);
} MATFUNC;

static MATFUNC *Matfunc;

/*
 Nucleation model for new phases
 */
char Nucleation(const PetscReal current_time, const PetscReal current_timestep, 
                const PetscReal temperature, const PetscReal volume, const PetscReal gv, 
                const PetscInt siteID, const AppCtx *user)
{
    NUCLEUS *currentnucleus = &user->nucleus[user->sitenucleusmapping[siteID]];
    PetscReal KbT = KBANGST2*temperature;
    PetscReal radius_c = 2.0*currentnucleus->gamma*volume/gv/ANGSTROM;
    if (   radius_c < currentnucleus->lengthscale*pow(volume,1.0/user->dim)
        && radius_c > currentnucleus->minsize) {
        PetscReal zeldovich = currentnucleus->atomicvolume
                            * sqrt(currentnucleus->gamma/KbT)
                            / (2.0*PETSC_PI*radius_c*radius_c);
        PetscReal beta = 4.0*PETSC_PI
                       * currentnucleus->D0*exp(-currentnucleus->migration/temperature)
                       * radius_c/currentnucleus->atomicvolume;
        PetscReal incubation_time = 1.0/(2.0*zeldovich*zeldovich*beta);
        PetscReal nucleation_probability = exp(- (4.0*PETSC_PI*currentnucleus->gamma*radius_c*radius_c)
                                               * currentnucleus->shapefactor
                                               / (3.0*KbT));
        PetscReal site_probability = current_timestep * zeldovich * beta
                                   * nucleation_probability * exp(-incubation_time/current_time);
        PetscReal random_number = (rand()/(double)RAND_MAX);
        if (random_number < site_probability) return 1;
    }
    return 0;
}

/*
 Chemenergy - CALPHAD 2SL model for chemical energy
 */
static void Chemenergy_calphad2sl(PetscReal *chemenergy, const PetscReal *sitefrac, const PetscReal temperature, const CHEMFE energy, const PetscInt numcomps)
{
    const CALPHAD2SL *currentcalphad2sl = &energy.calphad2sl;
    const RK *currentbinary, *currentternary;
    const PetscReal *sitefrac_p = &sitefrac[0], *sitefrac_q = &sitefrac[numcomps];
    
    *chemenergy = 0.0;
    for (PetscInt c=0; c<numcomps; c++) {
        *chemenergy +=  R_GAS_CONST*temperature*(  currentcalphad2sl->p*sitefrac_p[c]*log(sitefrac_p[c]) 
                                                 + currentcalphad2sl->q*sitefrac_q[c]*log(sitefrac_q[c]));
    }

    for (PetscInt cp=0; cp<numcomps; cp++) {
        for (PetscInt cq=0; cq<numcomps; cq++) {
            *chemenergy += SumTSeries(temperature,currentcalphad2sl->unary[cp*numcomps+cq])*sitefrac_p[cp]*sitefrac_q[cq];
        }
    }    

    currentbinary = &currentcalphad2sl->binaryp[0];
    for (PetscInt cp_i=0; cp_i<numcomps; cp_i++) {
        for (PetscInt cp_j=cp_i+1; cp_j<numcomps; cp_j++) {
            for (PetscInt cq=0; cq<numcomps; cq++, currentbinary++) {
                for (PetscInt rko=0; rko < currentbinary->n; rko++) {
                    *chemenergy += SumTSeries(temperature,currentbinary->enthalpy[rko])
                                 * sitefrac_p[cp_i]*sitefrac_p[cp_j]*sitefrac_q[cq]
                                 * FastPow(sitefrac_p[cp_i] - sitefrac_p[cp_j],rko);
                }
            }
        }
    }            
    currentbinary = &currentcalphad2sl->binaryq[0];
    for (PetscInt cp=0; cp<numcomps; cp++) {
        for (PetscInt cq_i=0; cq_i<numcomps; cq_i++) {
            for (PetscInt cq_j=cq_i+1; cq_j<numcomps; cq_j++, currentbinary++) {
                for (PetscInt rko=0; rko < currentbinary->n; rko++) {
                    *chemenergy += SumTSeries(temperature,currentbinary->enthalpy[rko])
                                 * sitefrac_p[cp]*sitefrac_q[cq_i]*sitefrac_q[cq_j]
                                 * FastPow(sitefrac_q[cq_i] - sitefrac_q[cq_j],rko);
                }
            }
        }
    }            

    currentternary = &currentcalphad2sl->ternaryp[0];
    for (PetscInt cp_i=0; cp_i<numcomps; cp_i++) {
        for (PetscInt cp_j=cp_i+1; cp_j<numcomps; cp_j++) {
            for (PetscInt cp_k=cp_j+1; cp_k<numcomps; cp_k++) {
                for (PetscInt cq=0; cq<numcomps; cq++, currentternary++) {
                    for (PetscInt rko=0; rko < currentternary->n; rko++) {
                        *chemenergy += SumTSeries(temperature,currentternary->enthalpy[rko])
                                     * sitefrac_p[cp_i]*sitefrac_p[cp_j]*sitefrac_p[cp_k]*sitefrac_q[cq];
                    }
                }
            }
        }
    }        
    currentternary = &currentcalphad2sl->ternaryq[0];
    for (PetscInt cp=0; cp<numcomps; cp++) {
        for (PetscInt cq_i=0; cq_i<numcomps; cq_i++) {
            for (PetscInt cq_j=cq_i+1; cq_j<numcomps; cq_j++) {
                for (PetscInt cq_k=cq_j+1; cq_k<numcomps; cq_k++, currentternary++) {
                    for (PetscInt rko=0; rko < currentternary->n; rko++) {
                        *chemenergy += SumTSeries(temperature,currentternary->enthalpy[rko])
                                     * sitefrac_q[cq_i]*sitefrac_q[cq_j]*sitefrac_q[cq_k]*sitefrac_p[cp];
                    }
                }
            }
        }
    }        
}

/*
 Chemenergy - CALPHADDIS model for chemical energy
 */
static void Chemenergy_calphaddis(PetscReal *chemenergy, const PetscReal *sitefrac, const PetscReal temperature, const CHEMFE energy, const PetscInt numcomps)
{
    const CALPHADDIS *currentcalphaddis = &energy.calphaddis;
    const RK *currentbinary = &currentcalphaddis->binary[0];
    const RK *currentternary = &currentcalphaddis->ternary[0];
    *chemenergy = SumTSeries(temperature,currentcalphaddis->ref);
    for (PetscInt ck=0; ck<numcomps; ck++) {
        *chemenergy += SumTSeries(temperature,currentcalphaddis->unary[ck])*sitefrac[ck] + R_GAS_CONST*temperature*sitefrac[ck]*log(sitefrac[ck]);
        for (PetscInt cj=ck+1; cj<numcomps; cj++,currentbinary++) {
            for (PetscInt rko=0; rko<currentbinary->n; rko++){
                *chemenergy += SumTSeries(temperature,currentbinary->enthalpy[rko])*sitefrac[ck]*sitefrac[cj]*FastPow(sitefrac[ck]-sitefrac[cj],rko);
            }
            for (PetscInt ci=cj+1; ci<numcomps; ci++,currentternary++) {
                for (PetscInt rko=0; rko<currentternary->n; rko++){
                    *chemenergy += SumTSeries(temperature,currentternary->enthalpy[rko])
                                 * sitefrac[ck]*sitefrac[cj]*sitefrac[ci]
                                 * (sitefrac[currentternary->i[rko]] + (1.0-sitefrac[ck]-sitefrac[cj]-sitefrac[ci])/3.0);
                }
            }
        }
    }
}

/*
 Chemenergy - quadratic model for chemical energy
 */
static void Chemenergy_quad(PetscReal *chemenergy, const PetscReal *sitefrac, const PetscReal temperature, const CHEMFE energy, const PetscInt numcomps)
{
    const QUAD *currentquad = &energy.quad;
    *chemenergy = SumTSeries(temperature,currentquad->ref);
    for (PetscInt c=0; c<numcomps; c++)
        *chemenergy += SumTSeries(temperature,currentquad->unary [c])
                     * (sitefrac[c] - SumTSeries(temperature,currentquad->ceq[c]))
                     + SumTSeries(temperature,currentquad->binary[c])
                     * (sitefrac[c] - SumTSeries(temperature,currentquad->ceq[c]))
                     * (sitefrac[c] - SumTSeries(temperature,currentquad->ceq[c]));
}

/*
 Chemenergy - none model for chemical energy
 */
static void Chemenergy_none(PetscReal *chemenergy, const PetscReal *sitefrac, const PetscReal temperature, const CHEMFE energy, const PetscInt numcomps)
{
    *chemenergy = 0.0;
}

/*
 Chemenergy - Chemical energy
 */
void Chemenergy(PetscReal *chemenergy, const PetscReal *sitefrac, const PetscReal *chempot, const PetscReal temperature, const uint16_t *phaseID, const AppCtx *user)
{
    for (PetscInt g=0; g<phaseID[0]; g++) {
        const MATERIAL *currentmaterial = &user->material[user->phasematerialmapping[phaseID[g+1]]];
        Matfunc[user->phasematerialmapping[phaseID[g+1]]].Chemenergy(&chemenergy[g],&sitefrac[g*SF_SIZE],temperature,currentmaterial->energy,user->ncp);
        chemenergy[g] /= currentmaterial->molarvolume;    
        for (PetscInt site=0; site<currentmaterial->nsites; site++) {  
            for (PetscInt c=0; c<user->ndp; c++) chemenergy[g] -= chempot[c]*currentmaterial->stochiometry[site]*sitefrac[g*SF_SIZE+site*user->ncp+c];
        }                             
    }
}

/*
 Chemenergy - CALPHAD 2SL model for chemical potential
 */
static void SitepotentialExplicit_calphad2sl(PetscReal *sitepot, const PetscReal *sitefrac, const PetscReal temperature, const CHEMFE energy, const PetscInt numcomps)
{
    const CALPHAD2SL *currentcalphad2sl = &energy.calphad2sl;
    const RK *currentbinary, *currentternary;
    const PetscReal *sitefrac_p = &sitefrac[0], *sitefrac_q = &sitefrac[numcomps];
    PetscReal *sitepot_p = &sitepot[0], *sitepot_q = &sitepot[numcomps];
    PetscReal workval_a, workval_b, workval_c;
    
    memset(sitepot_p,0,numcomps*sizeof(PetscReal));
    memset(sitepot_q,0,numcomps*sizeof(PetscReal));
    for (PetscInt cp=0; cp<numcomps; cp++) {
        for (PetscInt cq=0; cq<numcomps; cq++) {
            workval_a = SumTSeries(temperature,currentcalphad2sl->unary[cp*numcomps+cq]);
            sitepot_p[cp] += workval_a*sitefrac_q[cq];
            sitepot_q[cq] += workval_a*sitefrac_p[cp];
        }
    }    

    currentbinary = &currentcalphad2sl->binaryp[0];
    for (PetscInt cp_i=0; cp_i<numcomps; cp_i++) {
        for (PetscInt cp_j=cp_i+1; cp_j<numcomps; cp_j++) {
            for (PetscInt cq=0; cq<numcomps; cq++, currentbinary++) {
                for (PetscInt rko=0; rko < currentbinary->n; rko++) {
                    workval_a = SumTSeries(temperature,currentbinary->enthalpy[rko]);
                    workval_b = FastPow(sitefrac_p[cp_i] - sitefrac_p[cp_j],rko);
                    workval_c = (rko == 0) ? 0.0 : ((PetscReal) rko)*FastPow(sitefrac_p[cp_i] - sitefrac_p[cp_j],rko-1);
                    sitepot_p[cp_i] += workval_a*sitefrac_p[cp_j]*sitefrac_q[cq]*(workval_b + workval_c*sitefrac_p[cp_i]);
                    sitepot_p[cp_j] += workval_a*sitefrac_p[cp_i]*sitefrac_q[cq]*(workval_b - workval_c*sitefrac_p[cp_j]);
                    sitepot_q[cq  ] += workval_a*sitefrac_p[cp_i]               * workval_b            *sitefrac_p[cp_j] ;
                }
            }
        }
    }            
    currentbinary = &currentcalphad2sl->binaryq[0];
    for (PetscInt cp=0; cp<numcomps; cp++) {
        for (PetscInt cq_i=0; cq_i<numcomps; cq_i++) {
            for (PetscInt cq_j=cq_i+1; cq_j<numcomps; cq_j++, currentbinary++) {
                for (PetscInt rko=0; rko < currentbinary->n; rko++) {
                    workval_a = SumTSeries(temperature,currentbinary->enthalpy[rko]);
                    workval_b = FastPow(sitefrac_q[cq_i] - sitefrac_q[cq_j],rko);
                    workval_c = (rko == 0) ? 0.0 : ((PetscReal) rko)*FastPow(sitefrac_q[cq_i] - sitefrac_q[cq_j],rko-1);
                    sitepot_q[cq_i] += workval_a*sitefrac_q[cq_j]*sitefrac_p[cp]*(workval_b + workval_c*sitefrac_q[cq_i]);
                    sitepot_q[cq_j] += workval_a*sitefrac_q[cq_i]*sitefrac_p[cp]*(workval_b - workval_c*sitefrac_q[cq_j]);
                    sitepot_p[cp  ] += workval_a*sitefrac_q[cq_i]               * workval_b            *sitefrac_q[cq_j] ;
                }
            }
        }
    }            

    currentternary = &currentcalphad2sl->ternaryp[0];
    for (PetscInt cp_i=0; cp_i<numcomps; cp_i++) {
        for (PetscInt cp_j=cp_i+1; cp_j<numcomps; cp_j++) {
            for (PetscInt cp_k=cp_j+1; cp_k<numcomps; cp_k++) {
                for (PetscInt cq=0; cq<numcomps; cq++, currentternary++) {
                    for (PetscInt rko=0; rko < currentternary->n; rko++) {
                        workval_a = SumTSeries(temperature,currentternary->enthalpy[rko]);
                        sitepot_p[cp_i] += workval_a*sitefrac_p[cp_j]*sitefrac_p[cp_k]*sitefrac_q[cq  ];
                        sitepot_p[cp_j] += workval_a*sitefrac_p[cp_k]*sitefrac_p[cp_i]*sitefrac_q[cq  ];
                        sitepot_p[cp_k] += workval_a*sitefrac_p[cp_i]*sitefrac_p[cp_j]*sitefrac_q[cq  ];
                        sitepot_q[cq  ] += workval_a*sitefrac_p[cp_i]*sitefrac_p[cp_j]*sitefrac_p[cp_k];
                    }
                }
            }
        }
    }        
    currentternary = &currentcalphad2sl->ternaryq[0];
    for (PetscInt cp=0; cp<numcomps; cp++) {
        for (PetscInt cq_i=0; cq_i<numcomps; cq_i++) {
            for (PetscInt cq_j=cq_i+1; cq_j<numcomps; cq_j++) {
                for (PetscInt cq_k=cq_j+1; cq_k<numcomps; cq_k++, currentternary++) {
                    for (PetscInt rko=0; rko < currentternary->n; rko++) {
                        workval_a = SumTSeries(temperature,currentternary->enthalpy[rko]);
                        sitepot_q[cq_i] += workval_a*sitefrac_q[cq_j]*sitefrac_q[cq_k]*sitefrac_p[cp  ];
                        sitepot_q[cq_j] += workval_a*sitefrac_q[cq_k]*sitefrac_q[cq_i]*sitefrac_p[cp  ];
                        sitepot_q[cq_k] += workval_a*sitefrac_q[cq_i]*sitefrac_q[cq_j]*sitefrac_p[cp  ];
                        sitepot_p[cp  ] += workval_a*sitefrac_q[cq_i]*sitefrac_q[cq_j]*sitefrac_q[cq_k];
                    }
                }
            }
        }
    }        
}

/*
 Chemicalpotential - CALPHADDIS model for chemical potential
 */
static void SitepotentialExplicit_calphaddis(PetscReal *sitepot, const PetscReal *sitefrac, const PetscReal temperature, const CHEMFE energy, const PetscInt numcomps)
{
    const CALPHADDIS *currentcalphaddis = &energy.calphaddis;
    const RK *currentbinary, *currentternary;
    PetscReal workval_a, workval_b, workval_c;
    
    memset(sitepot,0,numcomps*sizeof(PetscReal));
    currentbinary = &currentcalphaddis->binary[0];
    for (PetscInt ci=0; ci<numcomps; ci++) {
        for (PetscInt cj=ci+1; cj<numcomps; cj++, currentbinary++) {
            for (PetscInt rko=0; rko < currentbinary->n; rko++) {
                workval_a = SumTSeries(temperature,currentbinary->enthalpy[rko]);
                workval_b = FastPow(sitefrac[ci] - sitefrac[cj],rko);
                workval_c = (rko == 0) ? 0.0 : ((PetscReal) rko)*FastPow(sitefrac[ci] - sitefrac[cj],rko-1);
                sitepot[ci] += workval_a*sitefrac[cj]*(workval_b + workval_c*sitefrac[ci]);
                sitepot[cj] += workval_a*sitefrac[ci]*(workval_b - workval_c*sitefrac[cj]);
            }
        }
    }            
    
    currentternary = &currentcalphaddis->ternary[0];
    for (PetscInt ci=0; ci<numcomps; ci++) {
        for (PetscInt cj=ci+1; cj<numcomps; cj++) {
            for (PetscInt ck=cj+1; ck<numcomps; ck++, currentternary++) {
                for (PetscInt rko=0; rko < currentternary->n; rko++) {
                    workval_a = SumTSeries(temperature,currentternary->enthalpy[rko]);
                    workval_b = sitefrac[currentternary->i[rko]] + (1.0-sitefrac[ck]-sitefrac[cj]-sitefrac[ci])/3.0;
                    sitepot[ci] += workval_a*sitefrac[cj]*sitefrac[ck]*(workval_b - sitefrac[ci]/3.0);
                    sitepot[cj] += workval_a*sitefrac[ck]*sitefrac[ci]*(workval_b - sitefrac[cj]/3.0);
                    sitepot[ck] += workval_a*sitefrac[ci]*sitefrac[cj]*(workval_b - sitefrac[ck]/3.0);
                    sitepot[currentternary->i[rko]] += workval_a*sitefrac[ci]*sitefrac[cj]*sitefrac[ck];
                }
            }
        }
    }        
}

/*
 Chemicalpotential - Quadratic model for chemical potential
 */
static void SitepotentialExplicit_quad(PetscReal *sitepot, const PetscReal *sitefrac, const PetscReal temperature, const CHEMFE energy, const PetscInt numcomps)
{
    memset(sitepot,0,numcomps*sizeof(PetscReal));
}

/*
 Chemicalpotential - None model for chemical potential
 */
static void SitepotentialExplicit_none(PetscReal *sitepot, const PetscReal *sitefrac, const PetscReal temperature, const CHEMFE energy, const PetscInt numcomps)
{
    memset(sitepot,0,numcomps*sizeof(PetscReal));
}

/*
 Chemicalpotential
 */
void SitepotentialExplicit(PetscReal *sitepot, const PetscReal *sitefrac, const PetscReal temperature, const uint16_t phaseID, const AppCtx *user)
{
    PetscReal sitepot_full[SF_SIZE];
    const MATERIAL *currentmaterial = &user->material[user->phasematerialmapping[phaseID]];
    Matfunc[user->phasematerialmapping[phaseID]].SitepotentialExplicit(sitepot_full,sitefrac,temperature,currentmaterial->energy,user->ncp);
    for (PetscInt site=0; site<currentmaterial->nsites; site++) {  
        for (PetscInt c=0; c<user->ndp; c++) {  
            sitepot[site*user->ndp+c] = (sitepot_full[site*user->ncp+c] - sitepot_full[site*user->ncp+user->ndp])/currentmaterial->molarvolume;
        } 
    }                             
}

/*
 Chemenergy - CALPHAD 2SL model for chemical potential
 */
static void SitepotentialImplicit_calphad2sl(PetscReal *sitepot, const PetscReal *sitefrac, const PetscReal temperature, const CHEMFE energy, const PetscInt numcomps)
{
    const CALPHAD2SL *currentcalphad2sl = &energy.calphad2sl;
    const PetscReal *sitefrac_p = &sitefrac[0], *sitefrac_q = &sitefrac[numcomps];
    PetscReal *sitepot_p = &sitepot[0], *sitepot_q = &sitepot[numcomps];
    
    for (PetscInt c=0; c<numcomps; c++) {
        sitepot_p[c] = R_GAS_CONST*temperature*currentcalphad2sl->p*log(sitefrac_p[c]);
        sitepot_q[c] = R_GAS_CONST*temperature*currentcalphad2sl->q*log(sitefrac_q[c]);
    }    
}

/*
 Chemicalpotential - CALPHADDIS model for chemical potential
 */
static void SitepotentialImplicit_calphaddis(PetscReal *sitepot, const PetscReal *sitefrac, const PetscReal temperature, const CHEMFE energy, const PetscInt numcomps)
{
    const CALPHADDIS *currentcalphaddis = &energy.calphaddis;
    for (PetscInt ck=0; ck<numcomps; ck++) {  
        sitepot[ck] = SumTSeries(temperature,currentcalphaddis->unary[ck]) + R_GAS_CONST*temperature*log(sitefrac[ck]);
    }
}

/*
 Chemicalpotential - Quadratic model for chemical potential
 */
static void SitepotentialImplicit_quad(PetscReal *sitepot, const PetscReal *sitefrac, const PetscReal temperature, const CHEMFE energy, const PetscInt numcomps)
{
    const QUAD *currentquad = &energy.quad;
    for (PetscInt c=0; c<numcomps-1; c++)
        sitepot[c] =     SumTSeries(temperature,currentquad->unary [c]) 
                   + 2.0*SumTSeries(temperature,currentquad->binary[c])
                   * (sitefrac[c] - SumTSeries(temperature,currentquad->ceq[c]));
}

/*
 Chemicalpotential - None model for chemical potential
 */
static void SitepotentialImplicit_none(PetscReal *sitepot, const PetscReal *sitefrac, const PetscReal temperature, const CHEMFE energy, const PetscInt numcomps)
{
    memset(sitepot,0,numcomps*sizeof(PetscReal));
}

/*
 Chemicalpotential
 */
void SitepotentialImplicit(PetscReal *sitepot, const PetscReal *sitefrac, const PetscReal temperature, const uint16_t phaseID, const AppCtx *user)
{
    PetscReal sitepot_full[SF_SIZE];
    const MATERIAL *currentmaterial = &user->material[user->phasematerialmapping[phaseID]];
    Matfunc[user->phasematerialmapping[phaseID]].SitepotentialImplicit(sitepot_full,sitefrac,temperature,currentmaterial->energy,user->ncp);
    for (PetscInt site=0; site<currentmaterial->nsites; site++) {  
        for (PetscInt c=0; c<user->ndp; c++) {  
            sitepot[site*user->ndp+c] = (sitepot_full[site*user->ncp+c] - sitepot_full[site*user->ncp+user->ndp])/currentmaterial->molarvolume;
        } 
    }                             
}

/*
 Chemicalpotential
 */
void Sitepotential(PetscReal *sitepot, const PetscReal *sitefrac, const PetscReal temperature, const uint16_t phaseID, const AppCtx *user)
{
    PetscReal sitepot_im[SP_SIZE], sitepot_ex[SP_SIZE];
    const MATERIAL *currentmaterial = &user->material[user->phasematerialmapping[phaseID]];
    SitepotentialImplicit(sitepot_im,sitefrac,temperature,phaseID,user);
    SitepotentialExplicit(sitepot_ex,sitefrac,temperature,phaseID,user);
    for (PetscInt site=0; site<currentmaterial->nsites; site++) {  
        for (PetscInt c=0; c<user->ndp; c++) {  
            sitepot[site*user->ndp+c] = sitepot_im[site*user->ndp+c] + sitepot_ex[site*user->ndp+c];
        } 
    }                             
}

/*
 Composition - Semi-implicit concentration per phase
 */
static void Sitefrac_calphad2sl(PetscReal *sitefrac, const PetscReal *sitepot_im, const PetscReal temperature, const CHEMFE energy, const PetscInt numcomps)
{
    const CALPHAD2SL *currentcalphad2sl = &energy.calphad2sl;
    PetscReal *sitefrac_p = &sitefrac[0], *sitefrac_q = &sitefrac[numcomps];
    const PetscReal *sitepot_p = &sitepot_im[0], *sitepot_q = &sitepot_im[numcomps-1];
    PetscReal RT = R_GAS_CONST*temperature, sum_p = 1.0, sum_q = 1.0;
    for (PetscInt c=0; c<numcomps-1; c++) {
        sitefrac_p[c] = exp(sitepot_p[c]/currentcalphad2sl->p/RT); sum_p += sitefrac_p[c];
        sitefrac_q[c] = exp(sitepot_q[c]/currentcalphad2sl->q/RT); sum_q += sitefrac_q[c];
    }
    sitefrac_p[numcomps-1] = 1.0; sitefrac_q[numcomps-1] = 1.0;
    for (PetscInt c=0; c<numcomps-1; c++) {
        sitefrac_p[c] /= sum_p; sitefrac_p[numcomps-1] -= sitefrac_p[c];
        sitefrac_q[c] /= sum_q; sitefrac_q[numcomps-1] -= sitefrac_q[c];
    } 
}

/*
 Composition - Semi-implicit concentration per phase
 */
static void Sitefrac_calphaddis(PetscReal *sitefrac, const PetscReal *sitepot_im, const PetscReal temperature, const CHEMFE energy, const PetscInt numcomps)
{
    const CALPHADDIS *currentcalphaddis = &energy.calphaddis;
    PetscReal RT = R_GAS_CONST*temperature, sum = 1.0;
    for (PetscInt c=0; c<numcomps-1; c++) {
        sitefrac[c] = exp((   sitepot_im[c] 
                            - SumTSeries(temperature,currentcalphaddis->unary[c         ]) 
                            + SumTSeries(temperature,currentcalphaddis->unary[numcomps-1]))/RT);
        sum += sitefrac[c];
    }
    sitefrac[numcomps-1] = 1.0;
    for (PetscInt c=0; c<numcomps-1; c++) {
        sitefrac[c] /= sum; sitefrac[numcomps-1] -= sitefrac[c];
    } 
}

/*
 Composition - Semi-implicit concentration per phase
 */
static void Sitefrac_quad(PetscReal *sitefrac, const PetscReal *sitepot_im, const PetscReal temperature, const CHEMFE energy, const PetscInt numcomps)
{
    const QUAD *currentquad = &energy.quad;
    sitefrac[numcomps-1] = 1.0;
    for (PetscInt c=0; c<numcomps-1; c++) {
        sitefrac[c] = SumTSeries(temperature,currentquad->ceq[c]) 
                    + 0.5*(sitepot_im[c] - SumTSeries(temperature,currentquad->unary[c]))
                    / SumTSeries(temperature,currentquad->binary[c]);
        sitefrac[numcomps-1] -= sitefrac[c];
    }
}

/*
 Composition - const concentration per phase
 */
static void Sitefrac_none(PetscReal *sitefrac, const PetscReal *sitepot_im, const PetscReal temperature, const CHEMFE energy, const PetscInt numcomps)
{
    memset(sitefrac,0,numcomps*sizeof(PetscReal));
    sitefrac[numcomps-1] = 1.0;
}

/*
 Composition - Semi-implicit concentration per phase
 */
void Sitefrac(PetscReal *sitefrac, const PetscReal *sitepot_im, const PetscReal temperature, const uint16_t phaseID, const AppCtx *user)
{
    const MATERIAL *currentmaterial = &user->material[user->phasematerialmapping[phaseID]];
    PetscReal sitepot_im_scaled[SP_SIZE];
    for (PetscInt site=0; site<currentmaterial->nsites; site++) {  
        for (PetscInt c=0; c<user->ndp; c++) {  
            sitepot_im_scaled[site*user->ndp+c] = currentmaterial->molarvolume*sitepot_im[site*user->ndp+c];
        } 
    }                             
    Matfunc[user->phasematerialmapping[phaseID]].Sitefrac(sitefrac,sitepot_im_scaled,temperature,currentmaterial->energy,user->ncp);
}

/*
 Composition tangent wrt chemical potential
 */
static void SitefracTangent_calphad2sl(PetscReal *sitefractangent, const PetscReal *sitefrac, const PetscReal temperature, const CHEMFE energy, const PetscInt numcomps)
{
    const CALPHAD2SL *currentcalphad2sl = &energy.calphad2sl;
    PetscReal *sitefractangent_p = &sitefractangent[0];
    PetscReal *sitefractangent_q = &sitefractangent[(numcomps-1)*(numcomps-1)];
    const PetscReal *sitefrac_p = &sitefrac[0], *sitefrac_q = &sitefrac[numcomps];
    
    PetscReal RT = R_GAS_CONST*temperature;
    memset(sitefractangent_p,0,(numcomps-1)*(numcomps-1)*sizeof(PetscReal));
    memset(sitefractangent_q,0,(numcomps-1)*(numcomps-1)*sizeof(PetscReal));
    for (PetscInt cj=0; cj<numcomps-1; cj++) {
        sitefractangent_p[cj*(numcomps-1)+cj] += sitefrac_p[cj]/currentcalphad2sl->p/RT;
        sitefractangent_q[cj*(numcomps-1)+cj] += sitefrac_q[cj]/currentcalphad2sl->q/RT;
        for (PetscInt ci=0; ci<numcomps-1; ci++) {
            sitefractangent_p[cj*(numcomps-1)+ci] -= sitefrac_p[cj]*sitefrac_p[ci]/currentcalphad2sl->p/RT;
            sitefractangent_q[cj*(numcomps-1)+ci] -= sitefrac_q[cj]*sitefrac_q[ci]/currentcalphad2sl->q/RT;
        }
    }        
}

/*
 Composition tangent wrt chemical potential
 */
static void SitefracTangent_calphaddis(PetscReal *sitefractangent, const PetscReal *sitefrac, const PetscReal temperature, const CHEMFE energy, const PetscInt numcomps)
{
    PetscReal RT = R_GAS_CONST*temperature;
    memset(sitefractangent,0,(numcomps-1)*(numcomps-1)*sizeof(PetscReal));
    for (PetscInt cj=0; cj<numcomps-1; cj++) {
        sitefractangent[cj*(numcomps-1)+cj] += sitefrac[cj]/RT;
        for (PetscInt ci=0; ci<numcomps-1; ci++) {
            sitefractangent[cj*(numcomps-1)+ci] -= sitefrac[cj]*sitefrac[ci]/RT;
        }
    }        
}

/*
 Composition tangent wrt chemical potential
 */
static void SitefracTangent_quad(PetscReal *sitefractangent, const PetscReal *sitefrac, const PetscReal temperature, const CHEMFE energy, const PetscInt numcomps)
{
    memset(sitefractangent,0,(numcomps-1)*(numcomps-1)*sizeof(PetscReal));
    const QUAD *currentquad = &energy.quad;
    for (PetscInt cj=0; cj<numcomps-1; cj++) {
        sitefractangent[cj*(numcomps-1)+cj] = 0.5/SumTSeries(temperature,currentquad->binary[cj]);
    }        
}

/*
 Composition tangent wrt chemical potential
 */
static void SitefracTangent_none(PetscReal *sitefractangent, const PetscReal *sitefrac, const PetscReal temperature, const CHEMFE energy, const PetscInt numcomps)
{
    memset(sitefractangent,0,(numcomps-1)*(numcomps-1)*sizeof(PetscReal));
}

/*
 Composition tangent wrt chemical potential
 */
void SitefracTangent(PetscReal *sitefractangent, const PetscReal *sitefrac, const PetscReal temperature, const uint16_t phaseID, const AppCtx *user)
{
    const MATERIAL *currentmaterial = &user->material[user->phasematerialmapping[phaseID]];
    Matfunc[user->phasematerialmapping[phaseID]].SitefracTangent(sitefractangent,sitefrac,temperature,currentmaterial->energy,user->ncp);
    for (PetscInt site=0; site<currentmaterial->nsites; site++) {  
        for (PetscInt c=0; c<user->ndp*user->ndp; c++) {  
            sitefractangent[site*user->ndp*user->ndp+c] *= currentmaterial->molarvolume;
        } 
    }                             
}

/*
 Composition mobility
 */
static void CompositionMobility_calphad2sl(PetscReal *mobilityc, const PetscReal *sitefrac, const PetscReal temperature, const CHEMFE energy, const PetscInt numcomps)
{
    memset(mobilityc,0,numcomps*sizeof(PetscReal));
    const CALPHAD2SL *currentcalphad2sl = &energy.calphad2sl;
    for (PetscInt ck=0; ck<numcomps; ck++) {
        mobilityc[ck] = currentcalphad2sl->mobilityc[ck];
    }        
}

/*
 Composition mobility
 */
static void CompositionMobility_calphaddis(PetscReal *mobilityc, const PetscReal *sitefrac, const PetscReal temperature, const CHEMFE energy, const PetscInt numcomps)
{
    PetscReal migration[numcomps];
    const CALPHADDIS *currentcalphaddis = &energy.calphaddis;
    MOBILITY *currentmobility = &currentcalphaddis->mobilityc[0];
    memset(mobilityc,0,numcomps*sizeof(PetscReal));
    for (PetscInt ck=0; ck<numcomps; ck++,currentmobility++) {
        RK *currentbinary = &currentmobility->binary[0];
        migration[ck] = 0.0;
        for (PetscInt cj=0; cj<numcomps; cj++) {
            migration[ck] += sitefrac[cj]*SumTSeries(temperature,currentmobility->unary[cj]);
            for (PetscInt ci=cj+1; ci<numcomps; ci++,currentbinary++) {
                for (PetscInt nrk=0; nrk < currentbinary->n; nrk++) {
                    migration[ck] += sitefrac[cj]*sitefrac[ci]
                                   * SumTSeries(temperature,currentbinary->enthalpy[nrk])
                                   * FastPow(sitefrac[cj]-sitefrac[ci],nrk);
                }
            }
        }
        mobilityc[ck] = sitefrac[ck]
                      * currentmobility->m0/R_GAS_CONST/temperature
                      * exp(migration[ck]/R_GAS_CONST/temperature);
    }        
}

/*
 Composition mobility
 */
static void CompositionMobility_quad(PetscReal *mobilityc, const PetscReal *sitefrac, const PetscReal temperature, const CHEMFE energy, const PetscInt numcomps)
{
    memset(mobilityc,0,numcomps*sizeof(PetscReal));
    const QUAD *currentquad = &energy.quad;
    for (PetscInt ck=0; ck<numcomps; ck++) {
        mobilityc[ck] = currentquad->mobilityc[ck];
    }        
}

/*
 Composition mobility
 */
static void CompositionMobility_none(PetscReal *mobilityc, const PetscReal *sitefrac, const PetscReal temperature, const CHEMFE energy, const PetscInt numcomps)
{
    memset(mobilityc,0,numcomps*sizeof(PetscReal));
}

/*
 Composition mobility
 */
void CompositionMobilityLatticeRef(PetscReal *mobilityc_latticeref, const PetscReal *sitefrac, const PetscReal temperature, const uint16_t phaseID, const AppCtx *user)
{
    const MATERIAL *currentmaterial = &user->material[user->phasematerialmapping[phaseID]];
    Matfunc[user->phasematerialmapping[phaseID]].CompositionMobility(mobilityc_latticeref,sitefrac,temperature,currentmaterial->energy,user->ncp);
    for (PetscInt c=0; c<user->ncp; c++) mobilityc_latticeref[c] *= currentmaterial->molarvolume;
}

/*
 Composition mobility
 */
void CompositionMobilityVolumeRef(PetscReal *mobilityc_volumeref, const PetscReal *mobilityc_latticeref, const AppCtx *user)
{
    memset(mobilityc_volumeref,0,user->ndp*user->ndp*sizeof(PetscReal));
    for (PetscInt ck=0; ck<user->ndp; ck++) {
        mobilityc_volumeref[ck*user->ndp+ck] += mobilityc_latticeref[ck];
//         for (PetscInt cj=0; cj<user->ndp; cj++) {
//             mobilityc_volumeref[ck*user->ndp+cj] -= mobilityc_latticeref[cj]/user->ncp; 
//         }        
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
            Matfunc[mat].Chemenergy = &Chemenergy_quad;
            Matfunc[mat].SitepotentialExplicit = &SitepotentialExplicit_quad;
            Matfunc[mat].SitepotentialImplicit = &SitepotentialImplicit_quad;
            Matfunc[mat].Sitefrac = &Sitefrac_quad;
            Matfunc[mat].SitefracTangent = &SitefracTangent_quad;
            Matfunc[mat].CompositionMobility = &CompositionMobility_quad;
        } else if (currentmaterial->model == CALPHADDIS_CHEMENERGY  ) {
            Matfunc[mat].Chemenergy = &Chemenergy_calphaddis;
            Matfunc[mat].SitepotentialExplicit = &SitepotentialExplicit_calphaddis;
            Matfunc[mat].SitepotentialImplicit = &SitepotentialImplicit_calphaddis;
            Matfunc[mat].Sitefrac = &Sitefrac_calphaddis;
            Matfunc[mat].SitefracTangent = &SitefracTangent_calphaddis;
            Matfunc[mat].CompositionMobility = &CompositionMobility_calphaddis;
        } else if (currentmaterial->model == CALPHAD2SL_CHEMENERGY  ) {
            Matfunc[mat].Chemenergy = &Chemenergy_calphad2sl;
            Matfunc[mat].SitepotentialExplicit = &SitepotentialExplicit_calphad2sl;
            Matfunc[mat].SitepotentialImplicit = &SitepotentialImplicit_calphad2sl;
            Matfunc[mat].Sitefrac = &Sitefrac_calphad2sl;
            Matfunc[mat].SitefracTangent = &SitefracTangent_calphad2sl;
            Matfunc[mat].CompositionMobility = &CompositionMobility_calphad2sl;
        } else if (currentmaterial->model == NONE_CHEMENERGY     ) {
            Matfunc[mat].Chemenergy = &Chemenergy_none;
            Matfunc[mat].SitepotentialExplicit = &SitepotentialExplicit_none;
            Matfunc[mat].SitepotentialImplicit = &SitepotentialImplicit_none;
            Matfunc[mat].Sitefrac = &Sitefrac_none;
            Matfunc[mat].SitefracTangent = &SitefracTangent_none;
            Matfunc[mat].CompositionMobility = &CompositionMobility_none;
        }
    }    
}
