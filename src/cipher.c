
static char help[] = "Solves multi phase field equations \n";

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <petscsf.h>
#include <petscdm.h>
#include <petscdmplex.h>
#include <petscdmforest.h>
#include <petscds.h>
#include <petscts.h>
#include <petscviewerhdf5.h>
#include "init.h"
#include "material.h"
#include "utility.h"
#include "typedef.h"

/* User-defined routines */
extern PetscErrorCode RHSFunction(TS,PetscReal,Vec,Vec,void *);
extern PetscErrorCode PostStep(TS);

/*
 RHSFunctionLocal - Evaluates RHS function, F(x)
 */

PetscErrorCode RHSFunction(TS ts,PetscReal ftime,Vec X,Vec F,void *ptr)
{
    PetscErrorCode    ierr;
    AppCtx            *user = (AppCtx*) ptr;
    Vec               localX, laplacian;
    PetscScalar       *fdof, *rhs, *lap, *offset;
    PetscScalar       *pcell, *pcellL, *pcellR, *ncell;
    PetscScalar       *plapl, *plaplL, *plaplR, *pfrhs, fluxp[PF_SIZE];
    PetscReal         *chempot, *chempotL, *chempotR;
    PetscReal         *composition, *compositionL, *compositionR;
    PetscReal         *mobilitycv, *mobilitycvL, *mobilitycvR;
    PetscReal         *sitefrac, *sitepot_ex, sitepot_im[SP_SIZE], chempot_interface[DP_SIZE];
    PetscScalar       *dlapl, *dlaplL, *dlaplR, *dprhs, *cprhs, fluxd[DP_SIZE];
    PetscScalar       *tlapl, *tlaplL, *tlaplR, *tmrhs, fluxt;
    uint16_t          slist[AS_SIZE], slistL[AS_SIZE], slistR[AS_SIZE], nlist[AS_SIZE], setintersection[AS_SIZE];
    const PetscInt    *scells, *cone;
    PetscReal         ffactor, cfactor, deltaL, deltaR, volL, volR;
    PetscInt          localcell, cell, localface, face; 
    PetscReal         interpolant[PF_SIZE], caplsource[PF_SIZE];
    PetscReal         chemsource[PF_SIZE], rhs_unconstrained[PF_SIZE];
    PetscBool         active[PF_SIZE];
    PetscReal         dmdc[user->ndp*user->ndp], dcdm[user->ndp*user->ndp];
    uint16_t          setintersect[AS_SIZE], injectionL[AS_SIZE], injectionR[AS_SIZE];
    PetscReal         nactivephases, rhsval, triplejunctionenergy;
    PetscInt          g, gi, gj, gk, c, s, interfacekj;
    MATERIAL          *currentmaterial;
    INTERFACE         *currentinterface;
    PetscReal         work_vec_PF[PF_SIZE], work_vec_DP[DP_SIZE], work_vec_MB[DP_SIZE*DP_SIZE];
    PetscReal         mobility_elem[user->ncp], composition_avg[user->ncp], compsource[user->ncp];
    PetscReal         work_vec_SP[SP_SIZE], work_vec_CP[PF_SIZE*DP_SIZE], work_vec_CTT[PF_SIZE*SP_SIZE];
    PetscReal         work_vec_CTPot[PF_SIZE*SP_SIZE*DP_SIZE], work_vec_CTEx[PF_SIZE*SP_SIZE*DP_SIZE];
    PetscReal         sitefrac_global[PF_SIZE*SF_SIZE*user->ninteriorcells]; 
    PetscReal         composition_global[user->ncp*user->ninteriorcells]; 
    PetscReal         mobilitycv_global[user->ncp*user->ninteriorcells];
    PetscReal         specific_heat, tconductivity[user->ninteriorcells];
    PetscReal         *temperature, *temperatureL, *temperatureR, temperatureavg;
    PetscReal         interface_mobility[PF_SIZE*PF_SIZE], interface_energy[PF_SIZE*PF_SIZE], interface_width[PF_SIZE*PF_SIZE]; 
    PetscReal         *gradient_matrix, dot_product;
    PetscInt          conesize, nsupp, supp, dim, dir, nleastsq;
    Vec               gradient_global[user->dim], gradient_local[user->dim];
    PetscScalar       *grad[user->dim], *grad_cell[user->dim], unitvec[user->dim];
    PetscScalar       *grad_cellL[user->dim], *grad_cellR[user->dim], grad_cellavg[user->dim][PF_SIZE];

    /* Gather FVM residuals */
    ierr = DMGetLocalVector(user->da_solution,&localX); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(user->da_solution,X,INSERT_VALUES,localX); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->da_solution,X,INSERT_VALUES,localX); CHKERRQ(ierr);
    ierr = VecGetArray(localX, &fdof);
    ierr = VecZeroEntries(F); CHKERRQ(ierr);
    ierr = VecGetArray(F, &rhs);
    ierr = DMGetLocalVector(user->da_solution,&laplacian); CHKERRQ(ierr);
    ierr = VecZeroEntries(laplacian); CHKERRQ(ierr);
    ierr = VecGetArray(laplacian,&lap); CHKERRQ(ierr);

    /* Precalculate gradient quantities */
    if (user->gradient_calculation) {
        for (dim=0; dim<user->dim; dim++) {
            ierr = DMGetGlobalVector(user->da_solution,&gradient_global[dim]); CHKERRQ(ierr);
            ierr = VecZeroEntries(gradient_global[dim]); CHKERRQ(ierr);
            ierr = VecGetArray(gradient_global[dim],&grad[dim]); CHKERRQ(ierr);
        }    
        for (localcell = 0; localcell < user->nlocalcells; ++localcell) {
            cell = user->localcells[localcell];

            /* get fields */
            offset = NULL;
            ierr = DMPlexPointLocalRef(user->da_solution, cell, fdof, &offset);
            F2IFUNC(slist,&offset[AS_OFFSET]);
            pcell = &offset[PF_OFFSET];
            for (dim=0; dim<user->dim; dim++) {
                offset = NULL;
                ierr = DMPlexPointGlobalRef(user->da_solution, cell, grad[dim], &offset);
                grad_cell[dim] = &offset[PF_OFFSET];
            }    

            /* calculate phase interpolants */
            EvalInterpolant(interpolant,pcell,slist[0]);

            /* pre-calculate calculate phase field gradients */
            gradient_matrix = &user->gradient_matrix[24*user->dim*localcell];
            for (g=0; g<slist[0]; g++) {
                for (dim=0; dim<user->dim; dim++) {
                    grad_cell[dim][g] = 0.0;
                    for (nleastsq=0; nleastsq<user->gradient_nleastsq[localcell]; nleastsq++)
                        grad_cell[dim][g] -= gradient_matrix[dim*user->gradient_nleastsq[localcell]+nleastsq]*pcell[g];
                }
            }        
            ierr = DMPlexGetConeSize(user->da_solution,cell,&conesize);
            ierr = DMPlexGetCone(user->da_solution,cell,&cone);
            nleastsq=0;
            for (face=0; face<conesize; face++) {
                ierr = DMPlexGetSupportSize(user->da_solution, cone[face], &nsupp);
                ierr = DMPlexGetSupport(user->da_solution, cone[face], &scells);
                for (supp=0; supp<nsupp; supp++) {
                    if (scells[supp] != cell && scells[supp] < user->ninteriorcells) {
                        offset = NULL;
                        ierr = DMPlexPointLocalRef(user->da_solution, scells[supp], fdof, &offset);
                        F2IFUNC(nlist,&offset[AS_OFFSET]);
                        ncell = &offset[PF_OFFSET];
                        SetIntersection(setintersection,slistL,slistR,slist,nlist);
                        for (g=0; g<setintersection[0]; g++) {
                            for (dim=0; dim<user->dim; dim++) {
                                grad_cell[dim][slistL[g]] += gradient_matrix[dim*user->gradient_nleastsq[localcell]+nleastsq]
                                                           * ncell[slistR[g]];
                            }
                        }        
                        nleastsq++;
                    }
                }
            }
        }
        for (dim=0; dim<user->dim; dim++) {
            ierr = VecRestoreArray(gradient_global[dim],&grad[dim]); CHKERRQ(ierr);
            ierr = DMGetLocalVector(user->da_solution,&gradient_local[dim]); CHKERRQ(ierr);
            ierr = DMGlobalToLocalBegin(user->da_solution,gradient_global[dim],INSERT_VALUES,gradient_local[dim]); CHKERRQ(ierr);
            ierr = DMGlobalToLocalEnd(user->da_solution,gradient_global[dim],INSERT_VALUES,gradient_local[dim]); CHKERRQ(ierr);
            ierr = DMRestoreGlobalVector(user->da_solution,&gradient_global[dim]); CHKERRQ(ierr);
            ierr = VecGetArray(gradient_local[dim],&grad[dim]); CHKERRQ(ierr);
        }    
    }

    /* Precalculate cell quantities */
    if (user->ndp) {
        for (cell = 0; cell < user->ninteriorcells; ++cell) {
            /* get fields */
            offset = NULL;
            ierr = DMPlexPointLocalRef(user->da_solution, cell, fdof, &offset); CHKERRQ(ierr);
            F2IFUNC(slist,&offset[AS_OFFSET]);
            pcell = &offset[PF_OFFSET];
            chempot = &offset[DP_OFFSET];
            sitepot_ex = &offset[EX_OFFSET];
            temperature =  &offset[TM_OFFSET];
            sitefrac = &sitefrac_global[cell*PF_SIZE*SF_SIZE];
            composition = &composition_global[cell*user->ncp];
            mobilitycv = &mobilitycv_global[cell*user->ncp];

            /* calculate phase interpolants */
            EvalInterpolant(interpolant,pcell,slist[0]);

            /* calculate cell quantities */
            memset(mobilitycv,0,user->ncp*sizeof(PetscReal));
            memset(composition,0,user->ncp*sizeof(PetscReal));
            memset(chempot_interface,0,user->ndp*sizeof(PetscReal));
            for (gk=0; gk<slist[0]; gk++) {
                for (gj=gk+1; gj<slist[0]; gj++) {
                    interfacekj = user->interfacelist[slist[gk+1]*user->npf+slist[gj+1]];
                    currentinterface = &user->interface[interfacekj];
                    if (currentinterface->potential) {
                        for (c=0; c<user->ndp; c++) chempot_interface[c] += interpolant[gk]*interpolant[gj]
                                                                          * currentinterface->potential[c];
                    }    
                }
            }    
            tconductivity[cell] = 0.0;
            for (g =0; g<slist[0];  g++) {
                currentmaterial = &user->material[user->phasematerialmapping[slist[g+1]]];
                for (s=0; s<currentmaterial->nsites; s++) {
                    for (c=0; c<user->ndp; c++) {
                        sitepot_im[s*user->ndp+c] = currentmaterial->stochiometry[s]*(chempot[c] - chempot_interface[c]) 
                                                  - sitepot_ex[g*SP_SIZE+s*user->ndp+c];
                    }
                }
                Sitefrac(&sitefrac[g*SF_SIZE],sitepot_im,(*temperature),slist[g+1],user);
                CompositionMobilityComponent(mobility_elem,&sitefrac[g*SF_SIZE],(*temperature),slist[g+1],user);
                for (c=0; c<user->ncp; c++) mobilitycv[c] += interpolant[g]*mobility_elem[c];
                for (s=0; s<currentmaterial->nsites; s++) {
                    for (c=0; c<user->ncp; c++) {
                        composition[c] += interpolant[g]
                                        * currentmaterial->stochiometry[s]
                                        * sitefrac[g*SF_SIZE+s*user->ncp+c];
                    }
                }
                if (currentmaterial->thermal.tconductivity.nTser) 
                  tconductivity[cell] += interpolant[g]*SumTSeries(*temperature,currentmaterial->thermal.tconductivity);
            }
            for (gk=0; gk<slist[0]; gk++) {
                for (gj=gk+1; gj<slist[0]; gj++) {
                    interfacekj = user->interfacelist[slist[gk+1]*user->npf+slist[gj+1]];
                    currentinterface = &user->interface[interfacekj];
                    if (currentinterface->mobilityc) {
                        for (c=0; c<user->ncp; c++) {
                            mobilitycv[c] *= exp(currentinterface->mobilityc[c]/R_GAS_CONST/(*temperature));
                        }
                    }    
                }
            }    
        }
    } else {
        for (cell = 0; cell < user->ninteriorcells; ++cell) {
            /* get fields */
            offset = NULL;
            ierr = DMPlexPointLocalRef(user->da_solution, cell, fdof, &offset); CHKERRQ(ierr);
            F2IFUNC(slist,&offset[AS_OFFSET]);
            pcell = &offset[PF_OFFSET];

            /* calculate phase interpolants */
            EvalInterpolant(interpolant,pcell,slist[0]);

            /* calculate cell quantities */
            tconductivity[cell] = 0.0;
            for (g =0; g<slist[0];  g++) {
                currentmaterial = &user->material[user->phasematerialmapping[slist[g+1]]];
                if (currentmaterial->thermal.tconductivity.nTser) 
                  tconductivity[cell] += interpolant[g]*SumTSeries(*temperature,currentmaterial->thermal.tconductivity);
            }
        }    
    }

    /* Loop over faces and compute flux contribution to the RHS */
    for (localface = 0; localface < user->nlocalfaces; ++localface) {
        face = user->localfaces[localface];

        /* skip boundary faces */
        DMPlexGetSupport(user->da_solution, face, &scells);
        if (scells[0] >= user->ninteriorcells || scells[1] >= user->ninteriorcells) continue;

        /* get neighbouring cell data */
        offset = NULL;
        ierr = DMPlexPointLocalRef(user->da_solution, scells[0], fdof, &offset); CHKERRQ(ierr);
        F2IFUNC(slistL,&offset[AS_OFFSET]);
        pcellL = &offset[PF_OFFSET];
        chempotL = &offset[DP_OFFSET];
        temperatureL = &offset[TM_OFFSET];
        mobilitycvL = &mobilitycv_global[scells[0]*user->ncp];
        compositionL = &composition_global[scells[0]*user->ncp];
        offset = NULL;
        ierr = DMPlexPointLocalRef(user->da_solution, scells[1], fdof, &offset); CHKERRQ(ierr);
        F2IFUNC(slistR,&offset[AS_OFFSET]);
        pcellR = &offset[PF_OFFSET];
        chempotR = &offset[DP_OFFSET];
        temperatureR = &offset[TM_OFFSET];
        mobilitycvR = &mobilitycv_global[scells[1]*user->ncp];
        compositionR = &composition_global[scells[1]*user->ncp];
        offset = NULL;
        ierr = DMPlexPointLocalRef(user->da_solution, scells[0], lap,  &offset); CHKERRQ(ierr);
        plaplL = &offset[PF_OFFSET];
        dlaplL = &offset[DP_OFFSET];
        tlaplL = &offset[TM_OFFSET];
        offset = NULL;
        ierr = DMPlexPointLocalRef(user->da_solution, scells[1], lap,  &offset); CHKERRQ(ierr);
        plaplR = &offset[PF_OFFSET];
        dlaplR = &offset[DP_OFFSET];
        tlaplR = &offset[TM_OFFSET];
        temperatureavg = ((*temperatureL) + (*temperatureL))/2.0;

        /* get geometric data */
        deltaL = user->cellgeom[scells[0]]; deltaR = user->cellgeom[scells[1]];
        volL = FastPow(deltaL,user->dim); volR = FastPow(deltaR,user->dim); 
        ffactor = 2.0*(deltaL < deltaR ? volL/deltaL : volR/deltaR)/(deltaL + deltaR);
        
        if (slistL[0] > 1 || slistL[1] > 1) {
            /* get common active phases */
            SetIntersection(setintersect,injectionL,injectionR,slistL,slistR);

            /* get fluxes on cell faces */
            memset(fluxp,0,setintersect[0]*sizeof(PetscReal));
            for (gk=0; gk<setintersect[0]; gk++) {
                fluxp[gk] = (pcellR[injectionR[gk]] - pcellL[injectionL[gk]]);
                plaplL[injectionL[gk]] += ffactor*fluxp[gk]/volL;
                plaplR[injectionR[gk]] -= ffactor*fluxp[gk]/volR;
            }
        }    
        
        if (user->ndp) {
            for (c=0; c<user->ncp; c++) {
                if (fabs(mobilitycvL[c]) > 1e-32 && fabs(mobilitycvR[c]) > 1e-32) {
                    mobility_elem[c] = 2.0*mobilitycvL[c]*mobilitycvR[c]/(mobilitycvL[c] + mobilitycvR[c]);
                } else {
                    mobility_elem[c] = 0.0;
                }    
                composition_avg[c] = 2.0*compositionL[c]*compositionR[c]/(compositionL[c] + compositionR[c]);
            }
            CompositionMobilityVolumeRef(work_vec_MB,mobility_elem,composition_avg,user);
            for (c=0; c<user->ndp; c++) {
                work_vec_DP[c] = (chempotR[c] - chempotL[c]);
            }
            MatVecMult_CIPHER(fluxd,work_vec_MB,work_vec_DP,user->ndp);
            for (c=0; c<user->ndp; c++) {
                dlaplL[c] += ffactor*fluxd[c]/volL;
                dlaplR[c] -= ffactor*fluxd[c]/volR;
            }
        }
        
        if (fabs(tconductivity[scells[0]]) > 1e-32 && fabs(tconductivity[scells[1]]) > 1e-32) {
            fluxt = 2.0*tconductivity[scells[0]]*tconductivity[scells[1]]
                  / (tconductivity[scells[0]] + tconductivity[scells[1]])
                  * ((*temperatureR) - (*temperatureL));
            *tlaplL += ffactor*fluxt/volL;
            *tlaplR -= ffactor*fluxt/volR;
        }   
    }    

    /* Loop over boundary faces and add boundary conditions */
    {
        DMLabel bclabel;
        IS bcIS;
        PetscInt bcface, nbcfaces, pstart, pend;
        const PetscInt *bcfaces;
        BOUNDARYCONDITIONS *currentboundary = &user->bcs[0];
        
        ierr = DMGetLabel(user->da_solution, "boundary", &bclabel); CHKERRQ(ierr);
        ierr = DMPlexGetHeightStratum(user->da_solution, 1, &pstart, &pend); CHKERRQ(ierr);
        for (PetscInt bc = 0; bc<user->nbcs; bc++, currentboundary++) {
            DMLabelGetStratumIS(bclabel, currentboundary->boundaryid, &bcIS);
            if (bcIS) {
                ISGetLocalSize(bcIS, &nbcfaces);
                ISGetIndices(bcIS, &bcfaces);
                for (bcface = 0; bcface < nbcfaces; ++bcface) {
                    if (bcfaces[bcface] < pstart || bcfaces[bcface] >= pend) continue;
                    ierr = DMPlexGetSupport(user->da_solution, bcfaces[bcface], &scells); CHKERRQ(ierr);
                    /* get neighbouring cell data */
                    offset = NULL;
                    ierr = DMPlexPointLocalRef(user->da_solution, scells[0], fdof, &offset); CHKERRQ(ierr);
                    F2IFUNC(slistL,&offset[AS_OFFSET]);
                    pcellL = &offset[PF_OFFSET];
                    chempotL = &offset[DP_OFFSET];
                    temperatureL = &offset[TM_OFFSET];
                    mobilitycvL = &mobilitycv_global[scells[0]*user->ncp];
                    compositionL = &composition_global[scells[0]*user->ncp];
            
                    /* get geometric data */
                    deltaL = user->cellgeom[scells[0]];
        
                    /* get common active phases */
                    if (!(currentboundary->chem_bctype == NONE_BC)) {
                        memset(fluxd,0,user->ndp*sizeof(PetscReal));
                        if (currentboundary->chem_bctype == NEUMANN_BC) {
                            for (c=0; c<user->ndp; c++) {
                                if (currentboundary->chem_bcbool[c]) {
                                    fluxd[c] =  currentboundary->chem_bcval[c] * deltaL;
                                }
                            }
                        } else {
                            for (c=0; c<user->ndp; c++) {
                                if (currentboundary->chem_bcbool[c]) {
                                    fluxd[c] = mobilitycvL[c]*(currentboundary->chem_bcval[c] - chempotL[c]);
                                }
                            }
                        }
                    }
                    if        (currentboundary->thermal_bctype == NEUMANN_BC  ) {
                        fluxt = tconductivity[scells[0]]* currentboundary->thermal_bcval * deltaL;
                    } else if (currentboundary->thermal_bctype == DIRICHLET_BC) {
                        fluxt = tconductivity[scells[0]]*(currentboundary->thermal_bcval - (*temperatureL));
                    }  

                    {
                        offset = NULL;
                        ierr = DMPlexPointLocalRef(user->da_solution, scells[0], lap,  &offset); CHKERRQ(ierr);
                        dlaplL = &offset[DP_OFFSET];
                        tlaplL = &offset[TM_OFFSET];
                        cfactor = 1.0/FastPow(deltaL,2);
                        if (!(currentboundary->chem_bctype == NONE_BC)) {
                            for (c=0; c<user->ndp; c++) dlaplL[c] += cfactor*fluxd[c];
                        }
                        if (!(currentboundary->thermal_bctype == NONE_BC)) *tlaplL += cfactor*fluxt;
                    }
                }    
                ISRestoreIndices(bcIS, &bcfaces);
                ISDestroy(&bcIS);
            }
        } 
    }   

    /* Loop over cells and compute local contribution to the RHS */
    for (localcell = 0; localcell < user->nlocalcells; ++localcell) {
        cell = user->localcells[localcell];

        /* get fields */
        offset = NULL;
        ierr = DMPlexPointLocalRef(user->da_solution, cell, fdof, &offset); CHKERRQ(ierr);
        F2IFUNC(slist,&offset[AS_OFFSET]);
        pcell = &offset[PF_OFFSET];
        chempot = &offset[DP_OFFSET];
        sitepot_ex = &offset[EX_OFFSET];
        temperature =  &offset[TM_OFFSET];
        sitefrac = &sitefrac_global[cell*PF_SIZE*SF_SIZE];
        composition = &composition_global[cell*user->ncp];
        offset = NULL;
        ierr = DMPlexPointLocalRef(user->da_solution, cell, lap,  &offset); CHKERRQ(ierr);
        plapl = &offset[PF_OFFSET];
        dlapl = &offset[DP_OFFSET];
        tlapl = &offset[TM_OFFSET];
        offset = NULL;
        ierr = DMPlexPointGlobalRef(user->da_solution,cell, rhs,  &offset); CHKERRQ(ierr);
        pfrhs = &offset[PF_OFFSET];
        dprhs = &offset[DP_OFFSET];
        cprhs = &offset[EX_OFFSET];
        tmrhs = &offset[TM_OFFSET];
        
        /* calculate phase interpolants */
        EvalInterpolant(interpolant,pcell,slist[0]);
                
        if (slist[0] > 1) {
            /* calculate interface quantities */ 
            for (gk=0; gk<slist[0]; gk++) {
                for (gj=gk+1; gj<slist[0]; gj++) {
                    interfacekj = gk*slist[0]+gj;
                    interface_energy  [interfacekj] = 1.0; 
                    interface_mobility[interfacekj] = 1.0;
                }
            }    
            if (user->gradient_calculation) {
                for (dim=0; dim<user->dim; dim++) {
                    offset = NULL;
                    ierr = DMPlexPointLocalRef(user->da_solution, cell, grad[dim], &offset);
                    grad_cell[dim] = &offset[PF_OFFSET];
                }    
                for (gk=0; gk<slist[0]; gk++) {
                    for (gj=gk+1; gj<slist[0]; gj++) {
                        for (dim=0, rhsval=0.0; dim<user->dim; dim++) {
                            unitvec[dim] = grad_cell[dim][gk] - grad_cell[dim][gj];
                            rhsval += unitvec[dim]*unitvec[dim];
                        }    
                        rhsval = sqrt(rhsval);
                        if (rhsval > 1e-16) for (dim=0; dim<user->dim; dim++) unitvec[dim] /= rhsval;
                        currentinterface = &user->interface[user->interfacelist[slist[gk+1]*user->npf+slist[gj+1]]];
                        interfacekj = gk*slist[0]+gj;
                        for (dir=0; dir<currentinterface->ienergy->n; dir++) {
                            for (dim=0, dot_product=0.0; dim<user->dim; dim++) 
                                dot_product += unitvec[dim]*currentinterface->ienergy->dir[user->dim*dir+dim];
                            interface_energy[interfacekj] += FastPow(dot_product,4)*currentinterface->ienergy->val[dir];
                        }
                        for (dir=0; dir<currentinterface->imobility->n; dir++) {
                            for (dim=0, dot_product=0.0; dim<user->dim; dim++) 
                                dot_product += unitvec[dim]*currentinterface->imobility->dir[user->dim*dir+dim];
                            interface_mobility[interfacekj] += FastPow(dot_product,4)*currentinterface->imobility->val[dir];
                        }
                    }
                }
            }
            for (gk=0; gk<slist[0]; gk++) {
                for (gj=gk+1; gj<slist[0]; gj++) {
                    currentinterface = &user->interface[user->interfacelist[slist[gk+1]*user->npf+slist[gj+1]]];
                    interfacekj = gk*slist[0]+gj;
                    interface_energy[interfacekj] *= currentinterface->ienergy->e->m0;
                    if (currentinterface->ienergy->e->unary->nTser)
                    interface_energy[interfacekj] *= exp(  SumTSeries((*temperature), currentinterface->ienergy->e->unary[0])
                                                         / R_GAS_CONST/(*temperature));
                    interface_mobility[interfacekj] *= currentinterface->imobility->m->m0;
                    if (currentinterface->imobility->m->unary->nTser)
                    interface_mobility[interfacekj] *= exp(  SumTSeries((*temperature), currentinterface->imobility->m->unary[0])
                                                           / R_GAS_CONST/(*temperature));
                    interface_width[interfacekj] = currentinterface->width;
                }
            }    

            /* phase capillary driving force */ 
            memset(caplsource,0,slist[0]*sizeof(PetscReal));
            for (gk=0; gk<slist[0]; gk++) {
                for (gj=gk+1; gj<slist[0]; gj++) {
                    interfacekj = gk*slist[0]+gj;
                    caplsource[gk] -= interface_energy[interfacekj]*(pcell[gj]/interface_width[interfacekj] + plapl[gj]*interface_width[interfacekj]/PETSC_PI/PETSC_PI);
                    caplsource[gj] -= interface_energy[interfacekj]*(pcell[gk]/interface_width[interfacekj] + plapl[gk]*interface_width[interfacekj]/PETSC_PI/PETSC_PI);
                    for (gi=gj+1; gi<slist[0]; gi++) {
                        triplejunctionenergy = 0.0;
                        triplejunctionenergy = triplejunctionenergy > interface_energy[gk*slist[0]+gj]/interface_width[gk*slist[0]+gj] ?
                                               triplejunctionenergy : interface_energy[gk*slist[0]+gj]/interface_width[gk*slist[0]+gj];
                        triplejunctionenergy = triplejunctionenergy > interface_energy[gk*slist[0]+gi]/interface_width[gk*slist[0]+gi] ?
                                               triplejunctionenergy : interface_energy[gk*slist[0]+gi]/interface_width[gk*slist[0]+gi];
                        triplejunctionenergy = triplejunctionenergy > interface_energy[gj*slist[0]+gi]/interface_width[gj*slist[0]+gi] ?
                                               triplejunctionenergy : interface_energy[gj*slist[0]+gi]/interface_width[gj*slist[0]+gi];
                        triplejunctionenergy *= user->solparams.junctionpenalty;
                        caplsource[gk] -= triplejunctionenergy*pcell[gj]*pcell[gi];
                        caplsource[gj] -= triplejunctionenergy*pcell[gk]*pcell[gi];
                        caplsource[gi] -= triplejunctionenergy*pcell[gk]*pcell[gj];
                    }
                }
                caplsource[gk] *= 8.0;
                caplsource[gk] += plapl[gk];
            }   

            /* phase chemical driving force */ 
            Chemenergy(chemsource,sitefrac,chempot,(*temperature),slist,user);
            memset(work_vec_DP,0,DP_SIZE*sizeof(PetscReal));
            for (g=0; g<slist[0]; g++) {
                currentmaterial = &user->material[user->phasematerialmapping[slist[g+1]]];
                for (s=0; s<currentmaterial->nsites; s++) {
                    for (c=0; c<user->ndp; c++) {
                        work_vec_DP[c] += interpolant[g]
                                        * currentmaterial->stochiometry[s]
                                        * sitefrac[g*SF_SIZE+s*user->ncp+c];
                    }
                }
            }
            for (gk=0; gk<slist[0]; gk++) {
                for (gj=gk+1; gj<slist[0]; gj++) {
                    currentinterface = &user->interface[user->interfacelist[slist[gk+1]*user->npf+slist[gj+1]]];
                    if (currentinterface->potential) {
                        rhsval = 0.0;
                        for (c=0; c<user->ndp; c++) rhsval += currentinterface->potential[c]*work_vec_DP[c];
                        chemsource[gk] += rhsval*interpolant[gj]; chemsource[gj] += rhsval*interpolant[gk];    
                    }    
                }
            }
            MatMulInterpolantDerivative(chemsource,pcell,slist[0]);

            /* build unconstrained RHS to calculate active set */ 
            nactivephases = 0.0;
            memset(rhs_unconstrained,0,slist[0]*sizeof(PetscReal));
            for (gk=0; gk<slist[0]; gk++) {
                for (gj=gk+1; gj<slist[0]; gj++) {
                    interfacekj = gk*slist[0]+gj;
                    rhsval = interface_mobility[interfacekj]*(  (caplsource[gk] - caplsource[gj])
                                                              - (chemsource[gk] - chemsource[gj])
                                                              * 8.0*sqrt(interpolant[gk]*interpolant[gj])/PETSC_PI);
                    rhs_unconstrained[gk] += rhsval; rhs_unconstrained[gj] -= rhsval;
                }
                active[gk] =    (pcell[gk] >       TOL && pcell            [gk] < 1.0 - TOL)
                             || (pcell[gk] <       TOL && rhs_unconstrained[gk] > 0.0      )
                             || (pcell[gk] > 1.0 - TOL && rhs_unconstrained[gk] < 0.0      );
                if (active[gk]) nactivephases += 1.0;
            }

            /* build constrained RHS from active set*/ 
            for (gk=0; gk<slist[0]; gk++) {
                if (active[gk]) {
                    for (gj=gk+1; gj<slist[0]; gj++) {
                        if (active[gj]) {
                            interfacekj = gk*slist[0]+gj;
                            rhsval = interface_mobility[interfacekj]*(  (caplsource[gk] - caplsource[gj])
                                                                      - (chemsource[gk] - chemsource[gj])
                                                                      * 8.0*sqrt(interpolant[gk]*interpolant[gj])/PETSC_PI);
                            pfrhs[gk] += rhsval; pfrhs[gj] -= rhsval;
                        }
                    }
                    pfrhs[gk] /= nactivephases;
                }        
            }
        }    
        
        /* temperature rate */
        specific_heat = 0.0; *tmrhs = *tlapl;
        for (g =0; g<slist[0];  g++) {
            currentmaterial = &user->material[user->phasematerialmapping[slist[g+1]]];
            if (currentmaterial->thermal.specific_heat.nTser) 
              specific_heat += interpolant[g]*SumTSeries(*temperature,currentmaterial->thermal.specific_heat);
            if (currentmaterial->thermal.latent_heat.nTser) 
              *tmrhs += pfrhs[g]*SumTSeries(*temperature,currentmaterial->thermal.latent_heat);
        }
        if (specific_heat) *tmrhs /= specific_heat;
        if (user->solparams.temperature_rate) *tmrhs += user->solparams.temperature_rate[user->solparams.currentloadcase];    
            
        if (user->ndp) {
            /* calculate explicit potential rate */
            for (g=0; g<slist[0]; g++) {
                currentmaterial = &user->material[user->phasematerialmapping[slist[g+1]]];
                SitepotentialExplicit(work_vec_SP,&sitefrac[g*SF_SIZE],(*temperature),slist[g+1],user);
                for (s=0; s<currentmaterial->nsites; s++) {
                    for (c=0; c<user->ndp; c++) {
                        cprhs[g*SP_SIZE+s*user->ndp+c] = currentmaterial->chempot_ex_kineticcoeff
                                                       * (work_vec_SP[s*user->ndp+c] - sitepot_ex[g*SP_SIZE+s*user->ndp+c]);
                    }
                }
            }

            /* calculate composition rate */
            memset(chempot_interface,0,user->ndp*sizeof(PetscReal));
            memset(work_vec_CP,0,PF_SIZE*DP_SIZE*sizeof(PetscReal));
            for (gk=0; gk<slist[0]; gk++) {
                for (gj=gk+1; gj<slist[0]; gj++) {
                    interfacekj = user->interfacelist[slist[gk+1]*user->npf+slist[gj+1]];
                    currentinterface = &user->interface[interfacekj];
                    if (currentinterface->potential) {
                        for (c=0; c<user->ndp; c++) {
                            chempot_interface[c] += interpolant[gk]*interpolant[gj]*currentinterface->potential[c];
                            work_vec_CP[gk*DP_SIZE+c] += currentinterface->potential[c]*interpolant[gj];
                            work_vec_CP[gj*DP_SIZE+c] += currentinterface->potential[c]*interpolant[gk];
                        }
                    }
                }
            }

            memset(dcdm,0,DP_SIZE*DP_SIZE*sizeof(PetscReal));
            for (g =0; g<slist[0];  g++) {
                currentmaterial = &user->material[user->phasematerialmapping[slist[g+1]]];
                for (s=0; s<currentmaterial->nsites; s++) {
                    for (c=0; c<user->ndp; c++) {
                        sitepot_im[s*user->ndp+c] = currentmaterial->stochiometry[s]
                                                  * (chempot[c] - chempot_interface[c]) 
                                                  - sitepot_ex[g*SP_SIZE+s*user->ndp+c];
                    }
                }
                SitefracTangent(&work_vec_CTPot[g*SP_SIZE*DP_SIZE], 
                                &work_vec_CTEx [g*SP_SIZE*DP_SIZE],
                                &work_vec_CTT  [g*SP_SIZE        ],
                                &sitefrac[g*SF_SIZE],sitepot_im,(*temperature),slist[g+1],user);
                for (s=0; s<currentmaterial->nsites; s++) {
                    for (c=0; c<user->ndp*user->ndp; c++) {
                        dcdm[c] += interpolant[g]
                                 * currentmaterial->stochiometry[s]
                                 * currentmaterial->stochiometry[s]
                                 * work_vec_CTPot[g*SP_SIZE*DP_SIZE+s*DP_SIZE*DP_SIZE+c];
                    }
                }        
            }
            
            for (g=0; g<slist[0]; g++) {
                memcpy(work_vec_DP,&work_vec_CP[g*DP_SIZE],DP_SIZE*sizeof(PetscReal));
                MatVecMult_CIPHER(&work_vec_CP[g*DP_SIZE],dcdm,work_vec_DP,DP_SIZE); 
            }    
            memset(work_vec_DP,0,DP_SIZE*sizeof(PetscReal));
            if (slist[0] > 1) {
                for (c=0; c<user->ndp; c++) {
                    for (g =0; g<slist[0];  g++) {
                        currentmaterial = &user->material[user->phasematerialmapping[slist[g+1]]];
                        work_vec_PF[g] = 0.0;
                        for (s=0; s<currentmaterial->nsites; s++) {
                            work_vec_PF[g] += currentmaterial->stochiometry[s]*sitefrac[g*SF_SIZE+s*user->ncp+c];
                        }
                        work_vec_PF[g] -= work_vec_CP[g*DP_SIZE+c];
                    }
                    MatMulInterpolantDerivative(work_vec_PF,pcell,slist[0]);
                    for (g =0; g<slist[0];  g++) {
                        work_vec_DP[c] -= pfrhs[g]*work_vec_PF[g];
                    }    
                }
            }
            Chemsource(compsource, composition, interpolant, slist, user);
            for (c=0; c<user->ndp; c++) work_vec_DP[c] += (dlapl[c] + compsource[c]);
            
            for (g =0; g<slist[0];  g++) {
                currentmaterial = &user->material[user->phasematerialmapping[slist[g+1]]];
                for (s=0; s<currentmaterial->nsites; s++) {
                    MatVecMult_CIPHER(&work_vec_CP[g*DP_SIZE],
                                      &work_vec_CTEx[g*SP_SIZE*DP_SIZE+s*DP_SIZE*DP_SIZE],
                                      &cprhs[g*SP_SIZE+s*DP_SIZE],DP_SIZE);
                    for (c=0; c<user->ndp; c++) {
                        work_vec_DP[c] -= interpolant[g]
                                        * currentmaterial->stochiometry[s]
                                        * (work_vec_CP[g*DP_SIZE+c] + work_vec_CTT[g*SP_SIZE+s*DP_SIZE+c]*(*tmrhs));
                    }
                }
            }
            Invertmatrix(dmdc,dcdm,DP_SIZE);  
            MatVecMult_CIPHER(dprhs,dmdc,work_vec_DP,DP_SIZE);  
        }
    }

    /* Restore precalculated gradient quantities */
    if (user->gradient_calculation) {
        for (dim=0; dim<user->dim; dim++) {
            ierr = VecRestoreArray(gradient_local[dim],&grad[dim]); CHKERRQ(ierr);
            ierr = DMRestoreLocalVector(user->da_solution,&gradient_local[dim]); CHKERRQ(ierr);
        }    
    }

    /* Restore FVM residuals */
    ierr = VecRestoreArray(localX, &fdof);
    ierr = DMRestoreLocalVector(user->da_solution,&localX); CHKERRQ(ierr);
    ierr = VecRestoreArray(F, &rhs);
    ierr = VecRestoreArray(laplacian,&lap); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(user->da_solution,&laplacian); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}    

/*
 PostStep - Update active phases and write output
 */
PetscErrorCode PostStep(TS ts)
{
    PetscErrorCode ierr;
    AppCtx         *user;
    Vec            solution, globalX, localX;
    PetscScalar    *fdof, *fadof, *offset, *offsetg;
    PetscScalar    *pcell, *dcell, *ccell;
    uint16_t       slist[AS_SIZE], glist[AS_SIZE], nlist[AS_SIZE];
    PetscInt       conesize, nsupp;
    const PetscInt *cone, *scells, *sitecells;
    PetscReal      phiUP[PF_SIZE], phiSS[PF_SIZE], sitepot_exSS[EX_SIZE], sitepot_im[SP_SIZE], chempot_interface[DP_SIZE];
    PetscInt       localcell, cell, face, supp;
    PetscInt       g, gk, gj, c, s;
    MATERIAL       *currentmaterial, *sitematerial;
    uint16_t       superset[AS_SIZE], setunion[AS_SIZE], setintersection[AS_SIZE], injectionL[AS_SIZE], injectionR[AS_SIZE];
    PetscReal      currenttime, currenttimestep, *temperature;
    DMLabel        slabel;
    IS             siteIS;
    PetscReal      cellvolume, interpolant[PF_SIZE];
    PetscInt       site, site_phase, nsitecells, interfacekj;
    NUCLEUS        *currentnucleus;
    INTERFACE      *currentinterface;
    
    PetscFunctionBeginUser;    
    ierr = TSGetSolution(ts, &solution);CHKERRQ(ierr);    
    ierr = TSGetApplicationContext(ts,&user);CHKERRQ(ierr);
    ierr = TSGetTime(ts,&currenttime);CHKERRQ(ierr);
    ierr = TSGetTimeStep(ts,&currenttimestep);CHKERRQ(ierr);
    
    /* Update nucleation events */
    if (user->nsites) {
        PetscReal gv_leaf[user->nsites], gv_root[user->nsites_local];
        PetscReal volume_leaf[user->nsites], volume_root[user->nsites_local];
        PetscReal temperature_leaf[user->nsites], temperature_root[user->nsites_local];
        PetscReal diffusivity_leaf[user->ncp][user->nsites], diffusivity_root[user->ncp][user->nsites_local];
        char deactive_leaf[user->nsites], deactive_root[user->nsites_local];
        PetscReal barrier, diffusivity[user->ncp];
        PetscInt  pstart, pend, ngdof;
        PetscSection gsection;

        ierr = DMPlexGetHeightStratum(user->da_solution, 0, &pstart, &pend); CHKERRQ(ierr);
        ierr = DMGetGlobalSection(user->da_solution,&gsection);CHKERRQ(ierr);
        ierr = DMGetLabel(user->da_solution, "site", &slabel);
        ierr = VecGetArray(solution,&fdof);
        memset(gv_leaf,0,user->nsites*sizeof(PetscReal));
        memset(volume_leaf,0,user->nsites*sizeof(PetscReal));
        memset(temperature_leaf,0,user->nsites*sizeof(PetscReal));
        for (c=0; c<user->ncp; c++) memset(diffusivity_leaf[c],0,user->nsites*sizeof(PetscReal));
        for (site = 0; site<user->nsites; site++) {
            if (user->siteactivity_global[site]) {
                DMLabelGetStratumIS(slabel, site+1, &siteIS);
                nsitecells = 0;
                if (siteIS) {
                    ISGetLocalSize(siteIS, &nsitecells);
                    ISGetIndices(siteIS, &sitecells);
                }    
                site_phase = user->sitephasemapping[site];
                currentnucleus = &user->nucleus[user->sitenucleusmapping[site]];
                sitematerial = &user->material[user->phasematerialmapping[site_phase]];
                for (cell = 0; cell < nsitecells; ++cell) {
                    if (sitecells[cell] < pstart || sitecells[cell] >= pend) continue;
                    ierr = PetscSectionGetDof(gsection, sitecells[cell], &ngdof);
                    if (ngdof > 0) {
                        /* get cell state */
                        offset = NULL;
                        ierr = DMPlexPointGlobalRef(user->da_solution, sitecells[cell], fdof, &offset);
                        F2IFUNC(slist,&offset[AS_OFFSET]);
                        pcell = &offset[PF_OFFSET];
                        dcell = &offset[DP_OFFSET];
                        ccell = &offset[EX_OFFSET];
                        temperature = &offset[TM_OFFSET];

                        /* check that active phases are a subset of matrix phases */
                        SetIntersection(setintersection,injectionL,injectionR,slist,currentnucleus->matrixlist);
                        if (setintersection[0] != slist[0]) {user->siteactivity_global[site] = 0; break;}
                
                        /* calculate phase interpolants */
                        EvalInterpolant(interpolant,pcell,slist[0]);
                
                        /* calculate nucleation barrier */
                        NucleationBarrier(&barrier,diffusivity,dcell,ccell,(*temperature),interpolant,site,slist,user);
                        cellvolume = FastPow(user->cellgeom[sitecells[cell]],user->dim); 
                        volume_leaf[site] += cellvolume;
                        temperature_leaf[site] += cellvolume * (*temperature);
                        gv_leaf[site] += cellvolume * barrier;
                        for (c=0; c<user->ncp; c++) diffusivity_leaf[c][site] += cellvolume * diffusivity[c];
                    }
                }    
                if (siteIS) {
                    ISRestoreIndices(siteIS, &sitecells);
                    ISDestroy(&siteIS);
                }
            }
        }       
        PetscSFReduceBegin(user->nucleation_sf,MPI_CHAR,user->siteactivity_global,user->siteactivity_local,MPI_LAND);
        PetscSFReduceEnd(user->nucleation_sf,MPI_CHAR,user->siteactivity_global,user->siteactivity_local,MPI_LAND);
        memset(user->siteactivity_global,0,user->nsites*sizeof(char));
        PetscSFBcastBegin(user->nucleation_sf,MPI_CHAR,user->siteactivity_local,user->siteactivity_global,MPI_REPLACE);
        PetscSFBcastEnd(user->nucleation_sf,MPI_CHAR,user->siteactivity_local,user->siteactivity_global,MPI_REPLACE);
    
        memset(gv_root,0,user->nsites_local*sizeof(PetscReal));
        PetscSFReduceBegin(user->nucleation_sf,MPIU_SCALAR,gv_leaf,gv_root,MPI_SUM);
        PetscSFReduceEnd(user->nucleation_sf,MPIU_SCALAR,gv_leaf,gv_root,MPI_SUM);
        for (c=0; c<user->ncp; c++)  {
            memset(diffusivity_root[c],0,user->nsites_local*sizeof(PetscReal));
            PetscSFReduceBegin(user->nucleation_sf,MPIU_SCALAR,diffusivity_leaf[c],diffusivity_root[c],MPI_SUM);
            PetscSFReduceEnd(user->nucleation_sf,MPIU_SCALAR,diffusivity_leaf[c],diffusivity_root[c],MPI_SUM);
        }
        memset(volume_root,0,user->nsites_local*sizeof(PetscReal));
        PetscSFReduceBegin(user->nucleation_sf,MPIU_SCALAR,volume_leaf,volume_root,MPI_SUM);
        PetscSFReduceEnd(user->nucleation_sf,MPIU_SCALAR,volume_leaf,volume_root,MPI_SUM);
        memset(temperature_root,0,user->nsites_local*sizeof(PetscReal));
        PetscSFReduceBegin(user->nucleation_sf,MPIU_SCALAR,temperature_leaf,temperature_root,MPI_SUM);
        PetscSFReduceEnd(user->nucleation_sf,MPIU_SCALAR,temperature_leaf,temperature_root,MPI_SUM);
        memset(deactive_root,0,user->nsites_local*sizeof(char));
        for (site=0; site<user->nsites_local; site++) {
            if (user->siteactivity_local[site]) {
                currentnucleus = &user->nucleus[user->sitenucleusmapping[user->siteoffset[user->worldrank]+site]];
                for (c=0, *diffusivity = LARGE; c<user->ncp; c++) {
                    if (currentnucleus->activesolutes[c]) *diffusivity = *diffusivity < diffusivity_root[c][site] 
                                                                       ? *diffusivity : diffusivity_root[c][site];
                }
                *diffusivity /= volume_root[site];
                gv_root[site] /= volume_root[site];
                temperature_root[site] /= volume_root[site];
                deactive_root[site] = NucleationEvent(currenttime,currenttimestep,
                                                      temperature_root[site],volume_root[site],gv_root[site],*diffusivity,
                                                      user->siteoffset[user->worldrank]+site,user);
            }    
        }  
        memset(deactive_leaf,0,user->nsites*sizeof(char));
        PetscSFBcastBegin(user->nucleation_sf,MPI_CHAR,deactive_root,deactive_leaf,MPI_REPLACE);
        PetscSFBcastEnd(user->nucleation_sf,MPI_CHAR,deactive_root,deactive_leaf,MPI_REPLACE);

        char reset = 0;
        for (site=0; site<user->nsites; site++) {
            if (deactive_leaf[site]) {
                reset = 1;
                user->siteactivity_global[site] = 0;
                site_phase = user->sitephasemapping[site];
                DMLabelGetStratumIS(slabel, site+1, &siteIS);
                nsitecells = 0;
                if (siteIS) {
                    ISGetLocalSize(siteIS, &nsitecells);
                    ISGetIndices(siteIS, &sitecells);
                }    
                sitematerial = &user->material[user->phasematerialmapping[site_phase]];
                for (cell = 0; cell < nsitecells; ++cell) {
                    if (sitecells[cell] < pstart || sitecells[cell] >= pend) continue;
                    ierr = PetscSectionGetDof(gsection, sitecells[cell], &ngdof);
                    if (ngdof > 0) {
                        /* get cell state */
                        offset = NULL;
                        ierr = DMPlexPointGlobalRef(user->da_solution, sitecells[cell], fdof, &offset); CHKERRQ(ierr);
                        F2IFUNC(slist,&offset[AS_OFFSET]);
                        pcell = &offset[PF_OFFSET];
                        dcell = &offset[DP_OFFSET];
                        ccell = &offset[EX_OFFSET];
                        temperature = &offset[TM_OFFSET];

                        /* update cell state */
                        slist[0] = 1; slist[1] = site_phase;
                        I2FFUNC(&offset[AS_OFFSET],slist);
                        pcell[0] = 1.0;
                        SitepotentialExplicit(ccell,sitematerial->c0,(*temperature),site_phase,user);
                    }
                }
                if (siteIS) {
                    ISRestoreIndices(siteIS, &sitecells);
                    ISDestroy(&siteIS);
                }
            }    
        }        
        ierr = VecRestoreArray(solution,&fdof);
        MPI_Allreduce(MPI_IN_PLACE,&reset,1,MPI_CHAR,MPI_LOR,PETSC_COMM_WORLD);
        if (reset) ierr = TSSetTimeStep(ts,user->solparams.mintimestep);
        PetscSFReduceBegin(user->nucleation_sf,MPI_CHAR,user->siteactivity_global,user->siteactivity_local,MPI_LAND);
        PetscSFReduceEnd(user->nucleation_sf,MPI_CHAR,user->siteactivity_global,user->siteactivity_local,MPI_LAND);
        memset(user->siteactivity_global,0,user->nsites*sizeof(char));
        PetscSFBcastBegin(user->nucleation_sf,MPI_CHAR,user->siteactivity_local,user->siteactivity_global,MPI_REPLACE);
        PetscSFBcastEnd(user->nucleation_sf,MPI_CHAR,user->siteactivity_local,user->siteactivity_global,MPI_REPLACE);
    }

    /* Determine active phase set, update composition */
    { 
     ierr = DMGetGlobalVector(user->da_solution,&globalX); CHKERRQ(ierr);
     ierr = VecCopy(solution,globalX);
     ierr = VecGetArray(globalX,&fadof);
     ierr = VecGetArray(solution,&fdof);
     for (localcell = 0; localcell < user->nlocalcells; ++localcell) {
         cell = user->localcells[localcell];

         /* get cell state */
         offset = NULL;
         ierr = DMPlexPointGlobalRef(user->da_solution,cell, fdof, &offset); CHKERRQ(ierr);
         F2IFUNC(slist,&offset[AS_OFFSET]);
         pcell = &offset[PF_OFFSET];
         offsetg = NULL;
         ierr = DMPlexPointGlobalRef(user->da_solution,cell, fadof, &offsetg); CHKERRQ(ierr);
         F2IFUNC(glist,&offsetg[AS_OFFSET]);

         if (slist[0] > 1) {
             /* project phase fields back to Gibbs simplex */
             memcpy(phiUP,pcell,slist[0]*sizeof(PetscReal));
             SimplexProjection(pcell,phiUP,slist[0]);
        
             /* update active set */
             glist[0] = 0;
             for (g=0; g<slist[0];  g++) {
                 if (pcell[g] > TOL) {
                     glist[++(glist[0])] = slist[g+1];
                 }
             }
             I2FFUNC(&offsetg[AS_OFFSET],glist);
    
         } else {
             /* project phase fields back to Gibbs simplex */
             pcell[0] = 1.0;
    
             /* update active set */
             glist[0] = 1; glist[1] = slist[1];
             I2FFUNC(&offsetg[AS_OFFSET],glist);
         }
     }
     ierr = VecRestoreArray(solution,&fdof);
     ierr = VecRestoreArray(globalX,&fadof);
     ierr = DMGetLocalVector(user->da_solution,&localX); CHKERRQ(ierr);
     ierr = DMGlobalToLocalBegin(user->da_solution,globalX,INSERT_VALUES,localX); CHKERRQ(ierr);
     ierr = DMGlobalToLocalEnd(user->da_solution,globalX,INSERT_VALUES,localX); CHKERRQ(ierr);
     ierr = DMRestoreGlobalVector(user->da_solution,&globalX); CHKERRQ(ierr);

     /* Determine active phase super set, reorder dofs in solution, initialize new phase frac & comp */
     ierr = VecGetArray(solution,&fdof); CHKERRQ(ierr);
     ierr = VecGetArray(localX,&fadof); CHKERRQ(ierr);
     for (localcell = 0; localcell < user->nlocalcells; ++localcell) {
         cell = user->localcells[localcell];

         /* get cell state */
         offsetg = NULL;
         ierr = DMPlexPointGlobalRef(user->da_solution, cell, fdof, &offsetg); CHKERRQ(ierr);
         F2IFUNC(slist,&offsetg[AS_OFFSET]);
         pcell = &offsetg[PF_OFFSET];
         dcell = &offsetg[DP_OFFSET];
         ccell = &offsetg[EX_OFFSET];
         temperature = &offsetg[TM_OFFSET];
         offset = NULL;
         ierr = DMPlexPointLocalRef(user->da_solution, cell, fadof, &offset); CHKERRQ(ierr);
         F2IFUNC(superset,&offset[AS_OFFSET]);

         /* add union of neighbouring cells */
         ierr = DMPlexGetConeSize(user->da_solution,cell,&conesize);
         ierr = DMPlexGetCone    (user->da_solution,cell,&cone    );
         for (face=0; face<conesize; face++) {
             ierr = DMPlexGetSupportSize(user->da_solution, cone[face], &nsupp);
             ierr = DMPlexGetSupport(user->da_solution, cone[face], &scells);
             for (supp=0; supp<nsupp; supp++) {
                 if (scells[supp] != cell && scells[supp] < user->ninteriorcells) {
                     offset = NULL;
                     ierr = DMPlexPointLocalRef(user->da_solution, scells[supp], fadof, &offset); CHKERRQ(ierr);
                     F2IFUNC(nlist,&offset[AS_OFFSET]);
                     SetUnion(setunion,injectionL,injectionR,superset,nlist);
                     memcpy(superset,setunion,(setunion[0]+1)*sizeof(uint16_t));
                 }
             }
         }

         /* initialize new phase states */
         for (g=0; g<superset[0];  g++) {
             phiSS[g] = 0.0;
             currentmaterial = &user->material[user->phasematerialmapping[superset[g+1]]];
             SitepotentialExplicit(&sitepot_exSS[g*SP_SIZE],currentmaterial->c0,(*temperature),superset[g+1],user);
         }
        
         /* reorder dofs to new active phase superset */
         SetIntersection(setintersection,injectionL,injectionR,slist,superset);
         for (g=0; g<setintersection[0];  g++) {
             phiSS[injectionR[g]] = pcell[injectionL[g]];
             memcpy(&sitepot_exSS[injectionR[g]*SP_SIZE],&ccell[injectionL[g]*SP_SIZE],SP_SIZE*sizeof(PetscReal));
         }
         memcpy(pcell,phiSS       ,superset[0]        *sizeof(PetscReal));
         memcpy(ccell,sitepot_exSS,superset[0]*SP_SIZE*sizeof(PetscReal));
         I2FFUNC(&offsetg[AS_OFFSET],superset);
         assert(superset[0] < MAXAP);
     }    
     ierr = VecRestoreArray(localX,&fadof); CHKERRQ(ierr);
     ierr = VecRestoreArray(solution,&fdof); CHKERRQ(ierr);
     ierr = DMRestoreLocalVector(user->da_solution,&localX); CHKERRQ(ierr);
     ierr = TSSetSolution(ts,solution); CHKERRQ(ierr);
    } 
    
    /* write output */
    ierr = TSGetStepNumber(ts, &user->solparams.step);CHKERRQ(ierr);
    if (user->solparams.step%user->solparams.outputfreq == 0) {
       Vec               Xout;
       PetscScalar       max, *xout, avgcomp[user->ncp], sitefrac[PF_SIZE*SF_SIZE], interpolant[PF_SIZE];
       char              name[256];
       PetscViewer       viewer;
       PetscInt          o;
       
       ierr = DMGetGlobalVector(user->da_output,&Xout); CHKERRQ(ierr);
       ierr = PetscObjectSetName((PetscObject) Xout, user->solparams.outfile);CHKERRQ(ierr);
       ierr = VecZeroEntries(Xout);
       ierr = VecGetArray(solution,&fdof); CHKERRQ(ierr);
       ierr = VecGetArray(Xout,&xout); CHKERRQ(ierr);
       for (localcell = 0; localcell < user->nlocalcells; ++localcell) {
           cell = user->localcells[localcell];
           /* get cell state */
           offset = NULL;
           ierr = DMPlexPointGlobalRef(user->da_solution, cell, fdof, &offset); CHKERRQ(ierr);
           F2IFUNC(slist,&offset[AS_OFFSET]);
           pcell = &offset[PF_OFFSET];
           dcell = &offset[DP_OFFSET];
           ccell = &offset[EX_OFFSET];
           temperature = &offset[TM_OFFSET];
           offset = NULL;
           ierr = DMPlexPointGlobalRef(user->da_output,cell, xout, &offset); CHKERRQ(ierr);
       
           EvalInterpolant(interpolant,pcell,slist[0]);

           for (o=0; o<user->noutputs; o++, offset++) {
               strcpy(name,user->outputname[o]);
               if (!strcmp(name,"phaseid")) {
                   max = -LARGE;
                   for (g=0; g<slist[0]; g++) {
                       if (interpolant[g] > max) {
                          max = interpolant[g];
                          *offset = (PetscScalar) slist[g+1];
                       }
                   }
               } else if (!strcmp(name,"matid")) {
                   max = -LARGE;
                   for (g=0; g<slist[0]; g++) {
                       if (interpolant[g] > max) {
                          max = interpolant[g];
                          *offset = (PetscScalar) user->phasematerialmapping[slist[g+1]];
                       }
                   }
               } else if (!strcmp(name,"interfaceid")) {
                   *offset = -1.0;
                   for (gk=0; gk<slist[0]; gk++) {
                       for (gj=gk+1; gj<slist[0]; gj++) {
                           *offset = (PetscScalar) user->interfacelist[slist[gk+1]*user->npf+slist[gj+1]];
                       }
                   }        
               } else if (strstr(name,"_chempot")) {
                   char *tok, *savetok;
                   tok = strtok_r(name, "_", &savetok);
                   for (c=0; c<user->ndp; c++) {
                       if (!strcmp(user->componentname[c],tok)) *offset = dcell[c];
                   }
               } else if (strstr(name,"_c")) {
                   memset(chempot_interface,0,user->ndp*sizeof(PetscReal));
                   for (gk=0; gk<slist[0]; gk++) {
                       for (gj=gk+1; gj<slist[0]; gj++) {
                           interfacekj = user->interfacelist[slist[gk+1]*user->npf+slist[gj+1]];
                           currentinterface = &user->interface[interfacekj];
                           if (currentinterface->potential) {
                               for (c=0; c<user->ndp; c++) chempot_interface[c] += interpolant[gk]*interpolant[gj]
                                                                                 * currentinterface->potential[c];
                           }    
                       }
                   }    
                   memset(avgcomp,0,user->ncp*sizeof(PetscReal));
                   for (g =0; g<slist[0];  g++) {
                       currentmaterial = &user->material[user->phasematerialmapping[slist[g+1]]];
                       for (s=0; s<currentmaterial->nsites; s++) {
                           for (c=0; c<user->ndp; c++) {
                               sitepot_im[s*user->ndp+c] = currentmaterial->stochiometry[s]*(dcell[c] - chempot_interface[c]) 
                                                         - ccell[g*SP_SIZE+s*user->ndp+c];
                           }
                       }
                       Sitefrac(&sitefrac[g*SF_SIZE],sitepot_im,(*temperature),slist[g+1],user);
                       for (s=0; s<currentmaterial->nsites; s++) {
                           for (c=0; c<user->ncp; c++) {
                               avgcomp[c] += interpolant[g]
                                           * currentmaterial->stochiometry[s]
                                           * sitefrac[g*SF_SIZE+s*user->ncp+c];
                           }
                       }
                   }
                   char *tok, *savetok;
                   tok = strtok_r(name, "_", &savetok);
                   for (c=0; c<user->ncp; c++) {
                       if (!strcmp(user->componentname[c],tok)) *offset = avgcomp[c];
                   }
               } else if (strstr(name,"_phi")) {
                   char *tok, *savetok;
                   tok = strtok_r(name, "_", &savetok);
                   *offset = 0.0;
                   for (g=0; g<slist[0]; g++) {
                       if (slist[g+1] == atoi(tok)) *offset = pcell[g];
                   }
               } else if (!strcmp(name,"temperature")) {
                   *offset = (*temperature);
               }
           }
       }
       ierr = VecRestoreArray(solution,&fdof); CHKERRQ(ierr);
       ierr = VecRestoreArray(Xout,&xout); CHKERRQ(ierr);

       sprintf(name, "%s_%d.vtu",user->solparams.outfile,user->solparams.step);
       PetscPrintf(PETSC_COMM_WORLD,"writing output at time %e to %s\n",currenttime,name);
       ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD, name, FILE_MODE_WRITE, &viewer);CHKERRQ(ierr);
       ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_VTK_VTU);CHKERRQ(ierr);
       ierr = VecView(Xout,viewer);
       ierr = PetscViewerDestroy(&viewer);
       ierr = DMRestoreGlobalVector(user->da_output,&Xout); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

static PetscErrorCode InitializeTS(DM dm, AppCtx *user, TS *ts)
{
  TSAdapt adapt;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = TSCreate(PetscObjectComm((PetscObject)dm), ts);CHKERRQ(ierr);
  ierr = TSSetType(*ts, TSRK);CHKERRQ(ierr);
  ierr = TSSetDM(*ts, dm);CHKERRQ(ierr);
  ierr = TSSetApplicationContext(*ts, user);CHKERRQ(ierr);
  ierr = TSSetRHSFunction(*ts, NULL, RHSFunction, user);CHKERRQ(ierr);
  ierr = TSSetPostStep(*ts,PostStep); CHKERRQ(ierr);
  ierr = TSGetAdapt(*ts,&adapt); CHKERRQ(ierr);
  ierr = TSAdaptSetType(adapt,TSADAPTDSP); CHKERRQ(ierr);
  ierr = TSAdaptSetStepLimits(adapt,user->solparams.mintimestep,user->solparams.maxtimestep); CHKERRQ(ierr);
  ierr = TSSetTolerances(*ts,user->solparams.abstol,NULL,user->solparams.reltol,NULL); CHKERRQ(ierr);
  ierr = TSSetFromOptions(*ts);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int main(int argc,char **args)
{
  /* user-defined work context */
  AppCtx         ctx;         
  /* solution vector */
  Vec            solution;
  /* time stepping context */
  TS             ts = NULL;
  /* numerical parameters */
  PetscErrorCode ierr;
      
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Initialize program
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = PetscInitialize(&argc,&args,(char*)0,help); CHKERRQ(ierr);
  MPI_Comm_rank(PETSC_COMM_WORLD,&ctx.worldrank);
  MPI_Comm_size(PETSC_COMM_WORLD,&ctx.worldsize);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Initialize problem parameters
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SetUpConfig(&ctx); CHKERRQ(ierr);
  utility_init(&ctx);
  material_init(&ctx);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create distributed mesh (DMPLEX) to manage 
   parallel vectors for the multi-phase field PDE
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = PetscOptionsInsertString(NULL,ctx.solparams.petscoptions);

  /* Create and distribute mesh */
  ierr = DMCreate(PETSC_COMM_WORLD,&ctx.da_solforest);CHKERRQ(ierr);
  ierr = DMSetType(ctx.da_solforest,ctx.dim == 2 ? DMP4EST : DMP8EST);CHKERRQ(ierr);
  ierr = DMForestSetTopology(ctx.da_solforest,"brick");CHKERRQ(ierr);
  ierr = DMForestSetInitialRefinement(ctx.da_solforest,ctx.amrparams.initrefine);CHKERRQ(ierr);
  ierr = DMForestSetMaximumRefinement(ctx.da_solforest,ctx.amrparams.maxnrefine);CHKERRQ(ierr);
  ierr = DMForestSetMinimumRefinement(ctx.da_solforest,ctx.amrparams.minnrefine);CHKERRQ(ierr);
  ierr = DMForestSetPartitionOverlap(ctx.da_solforest,1);CHKERRQ(ierr);
  ierr = DMSetFromOptions(ctx.da_solforest);CHKERRQ(ierr);
  ierr = DMSetUp(ctx.da_solforest);CHKERRQ(ierr);
  ierr = DMConvert(ctx.da_solforest,DMPLEX,&ctx.da_solution);CHKERRQ(ierr);
  ierr = DMLocalizeCoordinates(ctx.da_solution);CHKERRQ(ierr);
  
  /* Create discretisation for solution */
  {
    PetscFE solution_fe;
    ierr = PetscFECreateDefault(PETSC_COMM_WORLD,ctx.dim,
                                AS_SIZE+PF_SIZE+DP_SIZE+EX_SIZE+TM_SIZE,
                                PETSC_FALSE,NULL,PETSC_DEFAULT,&solution_fe);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) solution_fe, "solution");CHKERRQ(ierr);
    ierr = DMSetField(ctx.da_solution,0, NULL, (PetscObject) solution_fe);CHKERRQ(ierr);
    ierr = DMCreateDS(ctx.da_solution);CHKERRQ(ierr);
    ierr = PetscFEDestroy(&solution_fe);CHKERRQ(ierr);
    ierr = DMCopyDisc(ctx.da_solution,ctx.da_solforest);CHKERRQ(ierr);
  }
  
  /* Create discretisation for output */
  {
    PetscFE output_fe;
    ierr = DMClone(ctx.da_solution,&ctx.da_output);CHKERRQ(ierr);
    ierr = PetscFECreateDefault(PETSC_COMM_WORLD,ctx.dim,
                                ctx.noutputs,
                                PETSC_FALSE,NULL,PETSC_DEFAULT,&output_fe);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) output_fe, " output");CHKERRQ(ierr);
    ierr = DMSetField(ctx.da_output, 0, NULL, (PetscObject) output_fe);CHKERRQ(ierr);
    ierr = DMCreateDS(ctx.da_output);CHKERRQ(ierr);
    ierr = PetscFEDestroy(&output_fe);CHKERRQ(ierr);
  }

  /* Pre-calculate local cells and geometry */
  {
    DMLabel bclabel;
    PetscReal cvolume, fcentroid[ctx.dim];
    PetscSection lsection, gsection;
    PetscInt point, pstart, pend, ngdof, nldof, nsupp, nchild, dim, boundaryid;
    
    ierr = DMGetSection(ctx.da_solution,&lsection);CHKERRQ(ierr);
    ierr = DMGetGlobalSection(ctx.da_solution,&gsection);CHKERRQ(ierr);
    ierr = DMPlexGetHeightStratum(ctx.da_solution, 0, &pstart, &pend); CHKERRQ(ierr);
    ctx.localcells = malloc((pend-pstart)*sizeof(PetscInt)); 
    ctx.nlocalcells = 0;
    ctx.ninteriorcells = 0;
    for (point = pstart; point < pend; ++point) {
        ierr = PetscSectionGetDof(lsection, point, &nldof);
        ierr = PetscSectionGetDof(gsection, point, &ngdof);
        if (ngdof > 0) ctx.localcells[(ctx.nlocalcells)++] = point;
        if (nldof > 0) ++(ctx.ninteriorcells);
    }
    DMCreateLabel(ctx.da_solution, "boundary"); 
    DMGetLabel(ctx.da_solution, "boundary", &bclabel);
    ierr = DMPlexGetHeightStratum(ctx.da_solution, 1, &pstart, &pend); CHKERRQ(ierr);
    ctx.localfaces = malloc((pend-pstart)*sizeof(PetscInt)); ctx.nlocalfaces = 0;
    for (point = pstart; point < pend; ++point) {
        ierr = DMPlexGetSupportSize(ctx.da_solution, point, &nsupp);
        if (nsupp == 1) {
            DMPlexComputeCellGeometryFVM(ctx.da_solution, point, NULL, fcentroid, NULL);
            boundaryid = 0;
            for (dim=0; dim<ctx.dim; ++dim) {
                boundaryid++;
                if (fabs(fcentroid[dim]                ) < 1e-18) DMLabelSetValue(bclabel,point,boundaryid);
                boundaryid++;
                if (fabs(fcentroid[dim] - ctx.size[dim]) < 1e-18) DMLabelSetValue(bclabel,point,boundaryid);
            }
        }
        ierr = DMPlexGetTreeChildren(ctx.da_solution, point, &nchild, NULL);
        if (nsupp != 2 || nchild > 0) continue;
        ctx.localfaces[(ctx.nlocalfaces)++] = point;
    }
    ctx.cellgeom = malloc(ctx.ninteriorcells*sizeof(PetscReal));
    for (PetscInt cell = 0; cell < ctx.ninteriorcells; ++cell) {
        ierr = DMPlexComputeCellGeometryFVM(ctx.da_solution, cell, &cvolume, NULL, NULL);
        ctx.cellgeom[cell] = pow(cvolume,1.0/ctx.dim);
    }
  }
  
  /* Precalculate gradient matrix */
  if (ctx.gradient_calculation){
      PetscInt       conesize, nsupp;
      const PetscInt *cone, *scells;
      PetscInt       localcell, cell, face, supp, dim, dim_i, dim_j;
      PetscReal      val, centroid[ctx.dim], *gradient_matrix;
      PetscReal      a_mat[24*ctx.dim], ata_mat[ctx.dim*ctx.dim], atainv_mat[ctx.dim*ctx.dim];
      
      ctx.gradient_matrix = malloc(24*ctx.dim*ctx.nlocalcells*sizeof(PetscReal));
      ctx.gradient_nleastsq = malloc(ctx.nlocalcells*sizeof(PetscInt));
      for (localcell = 0; localcell < ctx.nlocalcells; ++localcell) {
          cell = ctx.localcells[localcell];
          ierr = DMPlexComputeCellGeometryFVM(ctx.da_solution, cell, NULL, centroid, NULL);
          ierr = DMPlexGetConeSize(ctx.da_solution,cell,&conesize);
          ierr = DMPlexGetCone(ctx.da_solution,cell,&cone);
          ctx.gradient_nleastsq[localcell] = 0;
          for (face=0; face<conesize; face++) {
              ierr = DMPlexGetSupportSize(ctx.da_solution, cone[face], &nsupp);
              ierr = DMPlexGetSupport(ctx.da_solution, cone[face], &scells);
              for (supp=0; supp<nsupp; supp++) {
                  if (scells[supp] != cell && scells[supp] < ctx.ninteriorcells) {
                      ierr = DMPlexComputeCellGeometryFVM( ctx.da_solution,scells[supp],NULL,
                                                          &a_mat[ctx.gradient_nleastsq[localcell]*ctx.dim],NULL);
                      for (dim=0; dim<ctx.dim; ++dim) a_mat[ctx.gradient_nleastsq[localcell]*ctx.dim+dim] -= centroid[dim];
                      (ctx.gradient_nleastsq[localcell])++;
                  }
              }
          }
          for (dim_i=0; dim_i<ctx.dim; ++dim_i) {
              for (dim_j=dim_i; dim_j<ctx.dim; ++dim_j) {
                  val = 0.0;
                  for (dim=0; dim<ctx.gradient_nleastsq[localcell]; ++dim) {
                      val += a_mat[dim*ctx.dim+dim_i]*a_mat[dim*ctx.dim+dim_j];
                  }
                  ata_mat[dim_i*ctx.dim+dim_j] = val;
                  ata_mat[dim_j*ctx.dim+dim_i] = val;
              }
          }
          Invertmatrix(atainv_mat,ata_mat,ctx.dim); 
          gradient_matrix = &ctx.gradient_matrix[24*ctx.dim*localcell];
          memset(gradient_matrix,0,ctx.dim*ctx.gradient_nleastsq[localcell]*sizeof(PetscReal));
          for (dim_i=0; dim_i<ctx.dim; ++dim_i) {
              for (dim_j=0; dim_j<ctx.gradient_nleastsq[localcell]; ++dim_j) {
                  for (dim=0; dim<ctx.dim; ++dim) 
                      gradient_matrix[dim_i*ctx.gradient_nleastsq[localcell]+dim_j] += atainv_mat[dim_i*ctx.dim+dim]
                                                                                     *      a_mat[dim_j*ctx.dim+dim];
              }
          }    
      }    
  }

  /* Add phase label */
  {
    DM celldm;
    Vec cellgeom;
    const PetscScalar *cgeom;
    PetscFVCellGeom *cg;
    PetscInt cell, localcell, off, dim, ioff[3] = {0};

    ierr = DMPlexGetDataFVM(ctx.da_solution, NULL, &cellgeom, NULL, NULL);
    ierr = VecGetDM(cellgeom,&celldm);
    ierr = VecGetArrayRead(cellgeom,&cgeom);CHKERRQ(ierr);
    ierr = DMCreateLabel(ctx.da_solution, "phase");CHKERRQ(ierr);
    ierr = DMCreateLabel(ctx.da_solution, "site");CHKERRQ(ierr);
    for (localcell = 0; localcell < ctx.nlocalcells; ++localcell) {
        cell = ctx.localcells[localcell];
        ierr = DMPlexPointLocalRead(celldm,cell,cgeom,&cg);
        for (dim=0; dim<ctx.dim; ++dim) {
            ioff[dim] = (PetscInt) (cg->centroid[dim]*ctx.resolution[dim]/ctx.size[dim]);
        }
        off = ioff[ctx.dim-1];
        for (dim=ctx.dim-1; dim>0; --dim) {
            off = (off*ctx.resolution[dim-1] + ioff[dim-1]);
        }
        ierr = DMSetLabelValue(ctx.da_solution, "phase", cell, ctx.voxelphasemapping[off]);CHKERRQ(ierr);
        ierr = DMSetLabelValue(ctx.da_solution, "site", cell, ctx.voxelsitemapping[off]);CHKERRQ(ierr);
    }
    ierr = VecRestoreArrayRead(cellgeom,&cgeom);CHKERRQ(ierr);
    free(ctx.voxelphasemapping);
    free(ctx.voxelsitemapping);
    ierr = DMCopyLabels(ctx.da_solution,ctx.da_solforest,PETSC_COPY_VALUES,PETSC_FALSE,DM_COPY_LABELS_REPLACE);CHKERRQ(ierr);
  }

  /* Set up star forest */
  if (ctx.nsites) {
      PetscInt *roots, sitepresent[ctx.nsites], site, rootctr;
      PetscInt pmin, pmax, pstart, pend;
      PetscReal sitesperproc;
      PetscSFNode *leaves;
      IS siteIS;
      DMLabel slabel;

      ierr = PetscSFCreate(PETSC_COMM_WORLD,&ctx.nucleation_sf);
      ierr = PetscSFSetFromOptions(ctx.nucleation_sf);
      ierr = DMGetLabel(ctx.da_solution, "site", &slabel);
      ierr = DMPlexGetHeightStratum(ctx.da_solution, 0, &pstart, &pend); CHKERRQ(ierr);
      memset(sitepresent,0,ctx.nsites*sizeof(PetscInt));
      for (site = 0, rootctr = 0; site<ctx.nsites; site++) {
          ierr = DMLabelGetStratumIS(slabel, site+1, &siteIS);CHKERRQ(ierr);
          if (siteIS) {
              ierr = ISGetMinMax(siteIS,&pmin,&pmax);CHKERRQ(ierr);
              if (pmin < pend && pmax >= pstart) {rootctr++; sitepresent[site] = 1;}
          }
      }
      ierr = PetscMalloc1(rootctr,&roots);CHKERRQ(ierr);
      ierr = PetscMalloc1(rootctr,&leaves);CHKERRQ(ierr);
      sitesperproc = ((PetscReal) ctx.nsites)/((PetscReal) ctx.worldsize);
      for (site = 0, rootctr = 0; site<ctx.nsites; site++) {
          if (sitepresent[site]) {
              roots[rootctr] = site;
              leaves[rootctr].rank = (PetscInt) (((PetscReal) site)/sitesperproc);
              leaves[rootctr].index = site - ctx.siteoffset[leaves[rootctr].rank];
              rootctr++;
          }
      }
      ierr = PetscSFSetGraph(ctx.nucleation_sf,ctx.nsites_local,rootctr,roots,PETSC_OWN_POINTER,leaves,PETSC_OWN_POINTER);
      ierr = PetscSFSetUp(ctx.nucleation_sf);
      memset(ctx.siteactivity_global,0,ctx.nsites*sizeof(char));
      PetscSFBcastBegin(ctx.nucleation_sf,MPI_CHAR,ctx.siteactivity_local,ctx.siteactivity_global,MPI_REPLACE);
      PetscSFBcastEnd(ctx.nucleation_sf,MPI_CHAR,ctx.siteactivity_local,ctx.siteactivity_global,MPI_REPLACE);
  }
  
  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Set up problem
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMCreateGlobalVector(ctx.da_solution, &solution);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) solution, "solution");CHKERRQ(ierr);
  ierr = SetUpProblem(solution,&ctx);CHKERRQ(ierr);
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Set up and perform time integration 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  for (PetscInt initcoarseniter = 0; initcoarseniter < ctx.amrparams.initcoarsen; ++initcoarseniter){
      PetscPrintf(PETSC_COMM_WORLD,"...initial refinement/coarsening (%d)...\n",initcoarseniter);
      /* Adapt mesh */
      PetscInt cell;
      PetscScalar *fdof, *offset;
      DMLabel adaptlabel;
      DM postsolforest;
      Vec lsolution, postvec;
  
      ierr = DMLabelCreate(PETSC_COMM_SELF,"adapt",&adaptlabel);CHKERRQ(ierr);
      ierr = DMLabelSetDefaultValue(adaptlabel,DM_ADAPT_COARSEN);CHKERRQ(ierr);

      ierr = DMGetLocalVector(ctx.da_solution,&lsolution);CHKERRQ(ierr);
      ierr = DMGlobalToLocalBegin(ctx.da_solution,solution,INSERT_VALUES,lsolution); CHKERRQ(ierr);
      ierr = DMGlobalToLocalEnd(ctx.da_solution,solution,INSERT_VALUES,lsolution); CHKERRQ(ierr);
      ierr = VecGetArray(lsolution, &fdof); CHKERRQ(ierr);
      for (cell = 0; cell < ctx.ninteriorcells; ++cell) {
          offset = NULL;
          ierr = DMPlexPointLocalRef(ctx.da_solution, cell, fdof, &offset); CHKERRQ(ierr);
          if (round(offset[AS_OFFSET]) > 1) {
              ierr = DMLabelSetValue(adaptlabel, cell, DM_ADAPT_REFINE); CHKERRQ(ierr);
          }    
      }
      ierr = VecRestoreArray(lsolution, &fdof); CHKERRQ(ierr);
      ierr = DMRestoreLocalVector(ctx.da_solution,&lsolution);CHKERRQ(ierr);
    
      if (ctx.nsites) { 
          PetscInt site, nsitecells;
          PetscInt pstart, pend;
          PetscInt face, supp, conesize, nsupp, labelvalue;
          const PetscInt *cone, *scells, *sitecells;
          IS siteIS;
          DMLabel slabel;
       
          ierr = DMPlexGetHeightStratum(ctx.da_solution, 0, &pstart, &pend); CHKERRQ(ierr);
          ierr = DMGetLabel(ctx.da_solution, "site", &slabel);
          for (site = 0; site<ctx.nsites; site++) {
              if (ctx.siteactivity_global[site]) {
                  DMLabelGetStratumIS(slabel, site+1, &siteIS);
                  if (siteIS) {
                      ISGetLocalSize(siteIS, &nsitecells);
                      ISGetIndices(siteIS, &sitecells);
                      for (cell = 0; cell < nsitecells; ++cell) {
                          if (sitecells[cell] < pstart || sitecells[cell] >= pend) continue;
                          /* check label of neighbouring cells */
                          ierr = DMPlexGetConeSize(ctx.da_solution,sitecells[cell],&conesize);
                          ierr = DMPlexGetCone    (ctx.da_solution,sitecells[cell],&cone    );
                          for (face=0; face<conesize; face++) {
                              ierr = DMPlexGetSupportSize(ctx.da_solution, cone[face], &nsupp);
                              ierr = DMPlexGetSupport(ctx.da_solution, cone[face], &scells);
                              for (supp=0; supp < nsupp; supp++) {
                                  ierr = DMGetLabelValue(ctx.da_solution, "site", scells[supp], &labelvalue);CHKERRQ(ierr);
                                  if (labelvalue != site+1) {
                                      ierr = DMLabelSetValue(adaptlabel, sitecells[cell], DM_ADAPT_REFINE); CHKERRQ(ierr);
                                  }
                              }
                          }
                      }
                      ISRestoreIndices(siteIS, &sitecells);
                      ISDestroy(&siteIS);
                  }
              }
          }
      } 

      postsolforest = NULL;
      ierr = DMAdaptLabel(ctx.da_solforest,adaptlabel,&postsolforest);CHKERRQ(ierr);
      if (postsolforest != NULL) {
          ierr = DMCreateGlobalVector(postsolforest, &postvec);CHKERRQ(ierr);
          ierr = DMForestTransferVec(ctx.da_solforest, solution, postsolforest, postvec, PETSC_TRUE, 0.0);CHKERRQ(ierr);
          ierr = VecDestroy(&solution);CHKERRQ(ierr);
          ierr = DMDestroy(&ctx.da_solforest);CHKERRQ(ierr);
          ierr = DMDestroy(&ctx.da_solution);CHKERRQ(ierr);
          ctx.da_solforest = postsolforest;
          ierr = DMForestSetAdaptivityForest(ctx.da_solforest,NULL);CHKERRQ(ierr);
          ierr = DMConvert(ctx.da_solforest,DMPLEX,&ctx.da_solution);CHKERRQ(ierr);
          ierr = DMCopyDisc(ctx.da_solforest,ctx.da_solution);CHKERRQ(ierr);
          ierr = DMCreateGlobalVector(ctx.da_solution,&solution);CHKERRQ(ierr);
          ierr = VecCopy(postvec,solution);CHKERRQ(ierr);
          ierr = VecDestroy(&postvec);CHKERRQ(ierr);

          /* Pre-calculate local cells and geometry */
          {
            PetscReal cvolume;
            PetscSection lsection, gsection;
            PetscInt point, pstart, pend, ngdof, nldof, nsupp, nchild, dim;
            
            ierr = DMGetSection(ctx.da_solution,&lsection);CHKERRQ(ierr);
            ierr = DMGetGlobalSection(ctx.da_solution,&gsection);CHKERRQ(ierr);
            ierr = DMPlexGetHeightStratum(ctx.da_solution, 0, &pstart, &pend); CHKERRQ(ierr);
            free(ctx.localcells);
            ctx.localcells = malloc((pend-pstart)*sizeof(PetscInt)); 
            ctx.nlocalcells = 0;
            ctx.ninteriorcells = 0;
            for (point = pstart; point < pend; ++point) {
                ierr = PetscSectionGetDof(lsection, point, &nldof);
                ierr = PetscSectionGetDof(gsection, point, &ngdof);
                if (ngdof > 0) ctx.localcells[(ctx.nlocalcells)++] = point;
                if (nldof > 0) ++(ctx.ninteriorcells);
            }
            ierr = DMPlexGetHeightStratum(ctx.da_solution, 1, &pstart, &pend); CHKERRQ(ierr);
            free(ctx.localfaces);
            ctx.localfaces = malloc((pend-pstart)*sizeof(PetscInt)); ctx.nlocalfaces = 0;
            for (point = pstart; point < pend; ++point) {
                ierr = DMPlexGetSupportSize(ctx.da_solution, point, &nsupp);
                ierr = DMPlexGetTreeChildren(ctx.da_solution, point, &nchild, NULL);
                if (nsupp != 2 || nchild > 0) continue;
                ctx.localfaces[(ctx.nlocalfaces)++] = point;
            }
            free(ctx.cellgeom);
            ctx.cellgeom = malloc(ctx.ninteriorcells*sizeof(PetscReal));
            for (PetscInt cell = 0; cell < ctx.ninteriorcells; ++cell) {
                ierr = DMPlexComputeCellGeometryFVM(ctx.da_solution, cell, &cvolume, NULL, NULL);
                ctx.cellgeom[cell] = pow(cvolume,1.0/ctx.dim);
            }
         }
          
          /* Precalculate gradient matrix */
          if (ctx.gradient_calculation){
              PetscInt       conesize, nsupp;
              const PetscInt *cone, *scells;
              PetscInt       localcell, cell, face, supp, dim, dim_i, dim_j;
              PetscReal      val, centroid[ctx.dim], *gradient_matrix;
              PetscReal      a_mat[24*ctx.dim], ata_mat[ctx.dim*ctx.dim], atainv_mat[ctx.dim*ctx.dim];
              
              free(ctx.gradient_matrix);ctx.gradient_matrix = malloc(24*ctx.dim*ctx.nlocalcells*sizeof(PetscReal));
              free(ctx.gradient_nleastsq);ctx.gradient_nleastsq = malloc(ctx.nlocalcells*sizeof(PetscInt));
              for (localcell = 0; localcell < ctx.nlocalcells; ++localcell) {
                  cell = ctx.localcells[localcell];
                  ierr = DMPlexComputeCellGeometryFVM(ctx.da_solution, cell, NULL, centroid, NULL);
                  ierr = DMPlexGetConeSize(ctx.da_solution,cell,&conesize);
                  ierr = DMPlexGetCone(ctx.da_solution,cell,&cone);
                  ctx.gradient_nleastsq[localcell] = 0;
                  for (face=0; face<conesize; face++) {
                      ierr = DMPlexGetSupportSize(ctx.da_solution, cone[face], &nsupp);
                      ierr = DMPlexGetSupport(ctx.da_solution, cone[face], &scells);
                      for (supp=0; supp<nsupp; supp++) {
                          if (scells[supp] != cell && scells[supp] < ctx.ninteriorcells) {
                              ierr = DMPlexComputeCellGeometryFVM( ctx.da_solution,scells[supp],NULL,
                                                                  &a_mat[ctx.gradient_nleastsq[localcell]*ctx.dim],NULL);
                              for (dim=0; dim<ctx.dim; ++dim) a_mat[ctx.gradient_nleastsq[localcell]*ctx.dim+dim] -= centroid[dim];
                              (ctx.gradient_nleastsq[localcell])++;
                          }
                      }
                  }
                  for (dim_i=0; dim_i<ctx.dim; ++dim_i) {
                      for (dim_j=dim_i; dim_j<ctx.dim; ++dim_j) {
                          val = 0.0;
                          for (dim=0; dim<ctx.gradient_nleastsq[localcell]; ++dim) {
                              val += a_mat[dim*ctx.dim+dim_i]*a_mat[dim*ctx.dim+dim_j];
                          }
                          ata_mat[dim_i*ctx.dim+dim_j] = val;
                          ata_mat[dim_j*ctx.dim+dim_i] = val;
                      }
                  }
                  Invertmatrix(atainv_mat,ata_mat,ctx.dim); 
                  gradient_matrix = &ctx.gradient_matrix[24*ctx.dim*localcell];
                  memset(gradient_matrix,0,ctx.dim*ctx.gradient_nleastsq[localcell]*sizeof(PetscReal));
                  for (dim_i=0; dim_i<ctx.dim; ++dim_i) {
                      for (dim_j=0; dim_j<ctx.gradient_nleastsq[localcell]; ++dim_j) {
                          for (dim=0; dim<ctx.dim; ++dim) 
                              gradient_matrix[dim_i*ctx.gradient_nleastsq[localcell]+dim_j] += atainv_mat[dim_i*ctx.dim+dim]
                                                                                             *      a_mat[dim_j*ctx.dim+dim];
                      }
                  }    
              }    
          }

          /* Set up star forest */
          if (ctx.nsites) {
              PetscInt *roots, sitepresent[ctx.nsites], site, rootctr;
              PetscInt pmin, pmax, pstart, pend;
              PetscReal sitesperproc;
              PetscSFNode *leaves;
              IS siteIS;
              DMLabel slabel;

              ierr = PetscSFDestroy(&ctx.nucleation_sf);
              ierr = PetscSFCreate(PETSC_COMM_WORLD,&ctx.nucleation_sf);
              ierr = PetscSFSetFromOptions(ctx.nucleation_sf);
              ierr = DMGetLabel(ctx.da_solution, "site", &slabel);
              ierr = DMPlexGetHeightStratum(ctx.da_solution, 0, &pstart, &pend); CHKERRQ(ierr);
              memset(sitepresent,0,ctx.nsites*sizeof(PetscInt));
              for (site = 0, rootctr = 0; site<ctx.nsites; site++) {
                  DMLabelGetStratumIS(slabel, site+1, &siteIS);
                  if (siteIS) {
                      ierr = ISGetMinMax(siteIS,&pmin,&pmax);CHKERRQ(ierr);
                      if (pmin < pend && pmax >= pstart) {rootctr++; sitepresent[site] = 1;}
                  }
              }
              ierr = PetscMalloc1(rootctr,&roots);CHKERRQ(ierr);
              ierr = PetscMalloc1(rootctr,&leaves);CHKERRQ(ierr);
              sitesperproc = ((PetscReal) ctx.nsites)/((PetscReal) ctx.worldsize);
              for (site = 0, rootctr = 0; site<ctx.nsites; site++) {
                  if (sitepresent[site]) {
                      roots[rootctr] = site;
                      leaves[rootctr].rank = (PetscInt) (((PetscReal) site)/sitesperproc);
                      leaves[rootctr].index = site - ctx.siteoffset[leaves[rootctr].rank];
                      rootctr++;
                  }
              }
              ierr = PetscSFSetGraph(ctx.nucleation_sf,ctx.nsites_local,rootctr,roots,PETSC_OWN_POINTER,leaves,PETSC_OWN_POINTER);
              ierr = PetscSFSetUp(ctx.nucleation_sf);
              memset(ctx.siteactivity_global,0,ctx.nsites*sizeof(char));
              PetscSFBcastBegin(ctx.nucleation_sf,MPI_CHAR,ctx.siteactivity_local,ctx.siteactivity_global,MPI_REPLACE);
              PetscSFBcastEnd(ctx.nucleation_sf,MPI_CHAR,ctx.siteactivity_local,ctx.siteactivity_global,MPI_REPLACE);
          }

         {
             PetscFE output_fe;
             ierr = DMDestroy(&ctx.da_output);CHKERRQ(ierr);
             ierr = DMClone(ctx.da_solution,&ctx.da_output);CHKERRQ(ierr);
             ierr = PetscFECreateDefault(PETSC_COMM_WORLD,ctx.dim,
                                         ctx.noutputs,
                                         PETSC_FALSE,NULL,PETSC_DEFAULT,&output_fe);CHKERRQ(ierr);
             ierr = PetscObjectSetName((PetscObject) output_fe, " output");CHKERRQ(ierr);
             ierr = DMSetField(ctx.da_output, 0, NULL, (PetscObject) output_fe);CHKERRQ(ierr);
             ierr = DMCreateDS(ctx.da_output);CHKERRQ(ierr);
             ierr = PetscFEDestroy(&output_fe);CHKERRQ(ierr);
          }
      }
  }    
  ierr = InitializeTS(ctx.da_solution, &ctx, &ts);CHKERRQ(ierr);

  PetscReal currenttime = 0.0, loadcasetime = 0.0, loadcaseendtime = 0.0;
  PetscInt  nsteps = 0;
  for (ctx.solparams.currentloadcase = 0;
       ctx.solparams.currentloadcase < ctx.solparams.nloadcases;
       loadcasetime += ctx.solparams.time[ctx.solparams.currentloadcase++]) {
      ctx.solparams.timestep = ctx.solparams.mintimestep;
      loadcaseendtime = loadcasetime+ctx.solparams.time[ctx.solparams.currentloadcase];
      for (;currenttime < loadcasetime+ctx.solparams.time[ctx.solparams.currentloadcase];) {
          ierr = TSSetStepNumber(ts,nsteps);CHKERRQ(ierr);
          ierr = TSSetTime(ts,currenttime);CHKERRQ(ierr);
          ierr = TSSetTimeStep(ts,ctx.solparams.timestep);CHKERRQ(ierr);
          ierr = TSSetMaxSteps(ts,nsteps+ctx.amrparams.amrinterval-nsteps%ctx.amrparams.amrinterval);CHKERRQ(ierr);
          ierr = TSSetMaxTime(ts,loadcaseendtime);CHKERRQ(ierr);
          ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);CHKERRQ(ierr);
          ierr = TSSolve(ts,solution);CHKERRQ(ierr);
          ierr = TSGetSolveTime(ts,&currenttime);CHKERRQ(ierr);
          ierr = TSGetStepNumber(ts,&nsteps);CHKERRQ(ierr);
          ierr = TSGetTimeStep(ts,&ctx.solparams.timestep);CHKERRQ(ierr);
          if (!(nsteps%ctx.amrparams.amrinterval)) {
              PetscPrintf(PETSC_COMM_WORLD,"...remeshing...\n");
              /* Adapt mesh */
              PetscInt cell;
              PetscScalar *fdof, *offset;
              DMLabel adaptlabel;
              DM postsolforest;
              Vec lsolution, postvec;
      
              ierr = DMLabelCreate(PETSC_COMM_SELF,"adapt",&adaptlabel);CHKERRQ(ierr);
              ierr = DMLabelSetDefaultValue(adaptlabel,DM_ADAPT_COARSEN);CHKERRQ(ierr);
              ierr = DMGetLocalVector(ctx.da_solution,&lsolution);CHKERRQ(ierr);
              ierr = DMGlobalToLocalBegin(ctx.da_solution,solution,INSERT_VALUES,lsolution); CHKERRQ(ierr);
              ierr = DMGlobalToLocalEnd(ctx.da_solution,solution,INSERT_VALUES,lsolution); CHKERRQ(ierr);
              ierr = VecGetArray(lsolution, &fdof); CHKERRQ(ierr);
              for (cell = 0; cell < ctx.ninteriorcells; ++cell) {
                  offset = NULL;
                  ierr = DMPlexPointLocalRef(ctx.da_solution, cell, fdof, &offset); CHKERRQ(ierr);
                  if (round(offset[AS_OFFSET]) > 1) {
                      ierr = DMLabelSetValue(adaptlabel, cell, DM_ADAPT_REFINE); CHKERRQ(ierr);
                  }    
              }
              ierr = VecRestoreArray(lsolution, &fdof); CHKERRQ(ierr);
              ierr = DMRestoreLocalVector(ctx.da_solution,&lsolution);CHKERRQ(ierr);

              if (ctx.nsites) { 
                  PetscInt site, nsitecells;
                  PetscInt pstart, pend;
                  PetscInt face, supp, conesize, nsupp, labelvalue;
                  const PetscInt *cone, *scells, *sitecells;
                  IS siteIS;
                  DMLabel slabel;
       
                  ierr = DMPlexGetHeightStratum(ctx.da_solution, 0, &pstart, &pend); CHKERRQ(ierr);
                  ierr = DMGetLabel(ctx.da_solution, "site", &slabel);
                  for (site = 0; site<ctx.nsites; site++) {
                      if (ctx.siteactivity_global[site]) {
                          DMLabelGetStratumIS(slabel, site+1, &siteIS);
                          if (siteIS) {
                              ISGetLocalSize(siteIS, &nsitecells);
                              ISGetIndices(siteIS, &sitecells);
                              for (cell = 0; cell < nsitecells; ++cell) {
                                  if (sitecells[cell] < pstart || sitecells[cell] >= pend) continue;
                                  /* check label of neighbouring cells */
                                  ierr = DMPlexGetConeSize(ctx.da_solution,sitecells[cell],&conesize);
                                  ierr = DMPlexGetCone    (ctx.da_solution,sitecells[cell],&cone    );
                                  for (face=0; face<conesize; face++) {
                                      ierr = DMPlexGetSupportSize(ctx.da_solution, cone[face], &nsupp);
                                      ierr = DMPlexGetSupport(ctx.da_solution, cone[face], &scells);
                                      for (supp=0; supp < nsupp; supp++) {
                                          ierr = DMGetLabelValue(ctx.da_solution, "site", scells[supp], &labelvalue);CHKERRQ(ierr);
                                          if (labelvalue != site+1) {
                                              ierr = DMLabelSetValue(adaptlabel, sitecells[cell], DM_ADAPT_REFINE); CHKERRQ(ierr);
                                          }
                                      }
                                  }
                              }
                              ISRestoreIndices(siteIS, &sitecells);
                              ISDestroy(&siteIS);
                          }
                      }
                  }
              } 

              postsolforest = NULL;
              ierr = DMAdaptLabel(ctx.da_solforest,adaptlabel,&postsolforest);CHKERRQ(ierr);
              if (postsolforest != NULL) {
                  ierr = DMCreateGlobalVector(postsolforest, &postvec);CHKERRQ(ierr);
                  ierr = DMForestTransferVec(ctx.da_solforest, solution, postsolforest, postvec, PETSC_TRUE, 0.0);CHKERRQ(ierr);
                  ierr = VecDestroy(&solution);CHKERRQ(ierr);
                  ierr = DMDestroy(&ctx.da_solforest);CHKERRQ(ierr);
                  ierr = DMDestroy(&ctx.da_solution);CHKERRQ(ierr);
                  ctx.da_solforest = postsolforest;
                  ierr = DMForestSetAdaptivityForest(ctx.da_solforest,NULL);CHKERRQ(ierr);
                  ierr = DMConvert(ctx.da_solforest,DMPLEX,&ctx.da_solution);CHKERRQ(ierr);
                  ierr = DMCopyDisc(ctx.da_solforest,ctx.da_solution);CHKERRQ(ierr);
                  ierr = DMCreateGlobalVector(ctx.da_solution,&solution);CHKERRQ(ierr);
                  ierr = VecCopy(postvec,solution);CHKERRQ(ierr);
                  ierr = VecDestroy(&postvec);CHKERRQ(ierr);

                  /* Pre-calculate local cells and geometry */
                  {
                    PetscReal cvolume;
                    PetscSection lsection, gsection;
                    PetscInt point, pstart, pend, ngdof, nldof, nsupp, nchild, dim;
            
                    ierr = DMGetSection(ctx.da_solution,&lsection);CHKERRQ(ierr);
                    ierr = DMGetGlobalSection(ctx.da_solution,&gsection);CHKERRQ(ierr);
                    ierr = DMPlexGetHeightStratum(ctx.da_solution, 0, &pstart, &pend); CHKERRQ(ierr);
                    free(ctx.localcells);
                    ctx.localcells = malloc((pend-pstart)*sizeof(PetscInt)); 
                    ctx.nlocalcells = 0;
                    ctx.ninteriorcells = 0;
                    for (point = pstart; point < pend; ++point) {
                        ierr = PetscSectionGetDof(lsection, point, &nldof);
                        ierr = PetscSectionGetDof(gsection, point, &ngdof);
                        if (ngdof > 0) ctx.localcells[(ctx.nlocalcells)++] = point;
                        if (nldof > 0) ++(ctx.ninteriorcells);
                    }
                    ierr = DMPlexGetHeightStratum(ctx.da_solution, 1, &pstart, &pend); CHKERRQ(ierr);
                    free(ctx.localfaces);
                    ctx.localfaces = malloc((pend-pstart)*sizeof(PetscInt)); ctx.nlocalfaces = 0;
                    for (point = pstart; point < pend; ++point) {
                        ierr = DMPlexGetSupportSize(ctx.da_solution, point, &nsupp);
                        ierr = DMPlexGetTreeChildren(ctx.da_solution, point, &nchild, NULL);
                        if (nsupp != 2 || nchild > 0) continue;
                        ctx.localfaces[(ctx.nlocalfaces)++] = point;
                    }
                    free(ctx.cellgeom);
                    ctx.cellgeom = malloc(ctx.ninteriorcells*sizeof(PetscReal));
                    for (PetscInt cell = 0; cell < ctx.ninteriorcells; ++cell) {
                        ierr = DMPlexComputeCellGeometryFVM(ctx.da_solution, cell, &cvolume, NULL, NULL);
                        ctx.cellgeom[cell] = pow(cvolume,1.0/ctx.dim);
                    }
                  }
          
                  /* Precalculate gradient matrix */
                  if (ctx.gradient_calculation){
                      PetscInt       conesize, nsupp;
                      const PetscInt *cone, *scells;
                      PetscInt       localcell, cell, face, supp, dim, dim_i, dim_j;
                      PetscReal      val, centroid[ctx.dim], *gradient_matrix;
                      PetscReal      a_mat[24*ctx.dim], ata_mat[ctx.dim*ctx.dim], atainv_mat[ctx.dim*ctx.dim];
                      
                      free(ctx.gradient_matrix);ctx.gradient_matrix = malloc(24*ctx.dim*ctx.nlocalcells*sizeof(PetscReal));
                      free(ctx.gradient_nleastsq);ctx.gradient_nleastsq = malloc(ctx.nlocalcells*sizeof(PetscInt));
                      for (localcell = 0; localcell < ctx.nlocalcells; ++localcell) {
                          cell = ctx.localcells[localcell];
                          ierr = DMPlexComputeCellGeometryFVM(ctx.da_solution, cell, NULL, centroid, NULL);
                          ierr = DMPlexGetConeSize(ctx.da_solution,cell,&conesize);
                          ierr = DMPlexGetCone(ctx.da_solution,cell,&cone);
                          ctx.gradient_nleastsq[localcell] = 0;
                          for (face=0; face<conesize; face++) {
                              ierr = DMPlexGetSupportSize(ctx.da_solution, cone[face], &nsupp);
                              ierr = DMPlexGetSupport(ctx.da_solution, cone[face], &scells);
                              for (supp=0; supp<nsupp; supp++) {
                                  if (scells[supp] != cell && scells[supp] < ctx.ninteriorcells) {
                                      ierr = DMPlexComputeCellGeometryFVM( ctx.da_solution,scells[supp],NULL,
                                                                          &a_mat[ctx.gradient_nleastsq[localcell]*ctx.dim],NULL);
                                      for (dim=0; dim<ctx.dim; ++dim) a_mat[ctx.gradient_nleastsq[localcell]*ctx.dim+dim] -= centroid[dim];
                                      (ctx.gradient_nleastsq[localcell])++;
                                  }
                              }
                          }
                          for (dim_i=0; dim_i<ctx.dim; ++dim_i) {
                              for (dim_j=dim_i; dim_j<ctx.dim; ++dim_j) {
                                  val = 0.0;
                                  for (dim=0; dim<ctx.gradient_nleastsq[localcell]; ++dim) {
                                      val += a_mat[dim*ctx.dim+dim_i]*a_mat[dim*ctx.dim+dim_j];
                                  }
                                  ata_mat[dim_i*ctx.dim+dim_j] = val;
                                  ata_mat[dim_j*ctx.dim+dim_i] = val;
                              }
                          }
                          Invertmatrix(atainv_mat,ata_mat,ctx.dim); 
                          gradient_matrix = &ctx.gradient_matrix[24*ctx.dim*localcell];
                          memset(gradient_matrix,0,ctx.dim*ctx.gradient_nleastsq[localcell]*sizeof(PetscReal));
                          for (dim_i=0; dim_i<ctx.dim; ++dim_i) {
                              for (dim_j=0; dim_j<ctx.gradient_nleastsq[localcell]; ++dim_j) {
                                  for (dim=0; dim<ctx.dim; ++dim) 
                                      gradient_matrix[dim_i*ctx.gradient_nleastsq[localcell]+dim_j] += atainv_mat[dim_i*ctx.dim+dim]
                                                                                                     *      a_mat[dim_j*ctx.dim+dim];
                              }
                          }    
                      }    
                  }

                  /* Set up star forest */
                  if (ctx.nsites){
                      PetscInt *roots, sitepresent[ctx.nsites], site, rootctr;
                      PetscInt pmin, pmax, pstart, pend;
                      PetscReal sitesperproc;
                      PetscSFNode *leaves;
                      IS siteIS;
                      DMLabel slabel;

                      ierr = PetscSFDestroy(&ctx.nucleation_sf);
                      ierr = PetscSFCreate(PETSC_COMM_WORLD,&ctx.nucleation_sf);
                      ierr = PetscSFSetFromOptions(ctx.nucleation_sf);
                      ierr = DMGetLabel(ctx.da_solution, "site", &slabel);
                      ierr = DMPlexGetHeightStratum(ctx.da_solution, 0, &pstart, &pend); CHKERRQ(ierr);
                      memset(sitepresent,0,ctx.nsites*sizeof(PetscInt));
                      for (site = 0, rootctr = 0; site<ctx.nsites; site++) {
                          DMLabelGetStratumIS(slabel, site+1, &siteIS);
                          if (siteIS) {
                              ierr = ISGetMinMax(siteIS,&pmin,&pmax);CHKERRQ(ierr);
                              if (pmin < pend && pmax >= pstart) {rootctr++; sitepresent[site] = 1;}
                          }
                      }
                      ierr = PetscMalloc1(rootctr,&roots);CHKERRQ(ierr);
                      ierr = PetscMalloc1(rootctr,&leaves);CHKERRQ(ierr);
                      sitesperproc = ((PetscReal) ctx.nsites)/((PetscReal) ctx.worldsize);
                      for (site = 0, rootctr = 0; site<ctx.nsites; site++) {
                          if (sitepresent[site]) {
                              roots[rootctr] = site;
                              leaves[rootctr].rank = (PetscInt) (((PetscReal) site)/sitesperproc);
                              leaves[rootctr].index = site - ctx.siteoffset[leaves[rootctr].rank];
                              rootctr++;
                          }
                      }
                      ierr = PetscSFSetGraph(ctx.nucleation_sf,ctx.nsites_local,rootctr,roots,PETSC_OWN_POINTER,leaves,PETSC_OWN_POINTER);
                      ierr = PetscSFSetUp(ctx.nucleation_sf);
                      memset(ctx.siteactivity_global,0,ctx.nsites*sizeof(char));
                      PetscSFBcastBegin(ctx.nucleation_sf,MPI_CHAR,ctx.siteactivity_local,ctx.siteactivity_global,MPI_REPLACE);
                      PetscSFBcastEnd(ctx.nucleation_sf,MPI_CHAR,ctx.siteactivity_local,ctx.siteactivity_global,MPI_REPLACE);
                  }
              
                  {
                    PetscFE output_fe;
                    ierr = DMDestroy(&ctx.da_output);CHKERRQ(ierr);
                    ierr = DMClone(ctx.da_solution,&ctx.da_output);CHKERRQ(ierr);
                    ierr = PetscFECreateDefault(PETSC_COMM_WORLD,ctx.dim,
                                                ctx.noutputs,
                                                PETSC_FALSE,NULL,PETSC_DEFAULT,&output_fe);CHKERRQ(ierr);
                    ierr = PetscObjectSetName((PetscObject) output_fe, " output");CHKERRQ(ierr);
                    ierr = DMSetField(ctx.da_output, 0, NULL, (PetscObject) output_fe);CHKERRQ(ierr);
                    ierr = DMCreateDS(ctx.da_output);CHKERRQ(ierr);
                    ierr = PetscFEDestroy(&output_fe);CHKERRQ(ierr);
                  }
                  ierr = TSDestroy(&ts);CHKERRQ(ierr);
                  ierr = InitializeTS(ctx.da_solution, &ctx, &ts);CHKERRQ(ierr);
              }
          }    
      }
  }

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Free work space.  All PETSc objects should be destroyed when they
   are no longer needed.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = VecDestroy(&solution);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx.da_solution);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}
