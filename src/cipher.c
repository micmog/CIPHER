
static char help[] = "Solves multi phase field equations \n";

#include <stdio.h>
#include <stdint.h>
#include <math.h>
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
    PetscScalar       *fdof, *rhs, *lap, *mat, *offset;
    PetscScalar       *pcell, *pcellL, *pcellR, *plapl, *plaplL, *plaplR, *pfrhs, fluxp[PF_SIZE];
    PetscScalar       *dcell, *dcellL, *dcellR, *dlapl, *dlaplL, *dlaplR, *dprhs, fluxd[DP_SIZE];
    uint16_t          slist[AS_SIZE], slistL[AS_SIZE], slistR[AS_SIZE];
    PetscScalar       *mcell;
    const PetscInt    *scells;
    PetscReal         ffactor, cfactor, deltaL, deltaR;
    PetscInt          localcell, cell, localface, face; 
    PetscReal         interpolant[PF_SIZE], caplflux[PF_SIZE], caplsource[PF_SIZE];
    PetscReal         chemsource[PF_SIZE], rhs_unconstrained[PF_SIZE], cvgdot_phase[PF_SIZE];
    PetscBool         active[PF_SIZE];
    PetscReal         chemicalpotential[DP_SIZE], composition[CP_SIZE];
    PetscReal         mobilitycv[DP_SIZE], cavgdot[DP_SIZE];
    PetscReal         dcdm[DP_SIZE*DP_SIZE], dmdc[DP_SIZE*DP_SIZE];
    uint16_t          setintersect[AS_SIZE], injectionL[AS_SIZE], injectionR[AS_SIZE];
    PetscReal         nactivephases, intval, rhsval, triplejunctionenergy;
    PetscInt          g, gi, gj, gk, c, ci, cj, interfacekj, interfaceji, interfaceki;
    INTERFACE         *currentinterface;
    
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
    ierr = VecGetArray(user->matstate,&mat); CHKERRQ(ierr);

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
        dcellL = &offset[DP_OFFSET];
        offset = NULL;
        ierr = DMPlexPointLocalRef(user->da_solution, scells[1], fdof, &offset); CHKERRQ(ierr);
        F2IFUNC(slistR,&offset[AS_OFFSET]);
        pcellR = &offset[PF_OFFSET];
        dcellR = &offset[DP_OFFSET];
        if (slistL[0] <= 1 && slistR[0] <= 1 && user->ncp <= 1) continue;

        /* get geometric data */
        deltaL = user->cellgeom[scells[0]]; deltaR = user->cellgeom[scells[0]];
        ffactor = 2.0*(deltaL < deltaR ? FastPow(deltaL,user->dim-1) : FastPow(deltaR,user->dim-1))/(deltaL + deltaR);
        
        /* get common active phases */
        if (slistL[0] > 1 || slistR[0] > 1) {
            SetIntersection(setintersect,injectionL,injectionR,slistL,slistR);
            for (g=0; g<setintersect[0]; g++) 
                fluxp[g] = (pcellR[injectionR[g]] - pcellL[injectionL[g]]);
        }
        for (c=0; c<user->ndp; c++) fluxd[c] = (dcellR[c] - dcellL[c]);
        
        {
            offset = NULL;
            ierr = DMPlexPointLocalRef(user->da_solution, scells[0], lap,  &offset); CHKERRQ(ierr);
            plaplL = &offset[PF_OFFSET];
            dlaplL = &offset[DP_OFFSET];
            cfactor = ffactor/FastPow(deltaL,user->dim);
            if (slistL[0] > 1) {
                for (g=0; g<setintersect[0]; g++)
                    plaplL[injectionL[g]] += cfactor*fluxp[g];
            }
            for (c=0; c<user->ndp; c++) dlaplL[c] += cfactor*fluxd[c];
        }
        {
            offset = NULL;
            ierr = DMPlexPointLocalRef(user->da_solution, scells[1], lap,  &offset); CHKERRQ(ierr);
            plaplR = &offset[PF_OFFSET];
            dlaplR = &offset[DP_OFFSET];
            cfactor = ffactor/FastPow(deltaR,user->dim);
            if (slistR[0] > 1) {
                for (g=0; g<setintersect[0]; g++)
                    plaplR[injectionR[g]] -= cfactor*fluxp[g];
            }
            for (c=0; c<user->ndp; c++) dlaplR[c] -= cfactor*fluxd[c];
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
        dcell = &offset[DP_OFFSET];
        offset = NULL;
        ierr = DMPlexPointLocalRef(user->da_solution, cell, lap,  &offset); CHKERRQ(ierr);
        plapl = &offset[PF_OFFSET];
        dlapl = &offset[DP_OFFSET];
        offset = NULL;
        ierr = DMPlexPointGlobalRef(user->da_solution,cell, rhs,  &offset); CHKERRQ(ierr);
        pfrhs = &offset[PF_OFFSET];
        dprhs = &offset[DP_OFFSET];
        mcell = NULL;
        ierr = DMPlexPointGlobalRef(user->da_matstate,cell, mat,  &mcell); CHKERRQ(ierr);
        
        /* calculate phase interpolants */
        EvalInterpolant(interpolant,pcell,slist[0]);
                
        /* update material state */
        memcpy(chemicalpotential,dcell,user->ndp*sizeof(PetscScalar));
        
        /* get chemical part of diffusion potential */
        for (gk=0; gk<slist[0]; gk++) {
            for (gj=gk+1; gj<slist[0]; gj++) {
                interfacekj = user->interfacelist[slist[gk+1]*user->npf + slist[gj+1]];
                currentinterface = &user->interface[interfacekj];
                intval = sqrt(interpolant[gk]*interpolant[gj]);
                for (c=0; c<user->ndp; c++) {
                    chemicalpotential[c] -= intval*currentinterface->potential[c];
                }
            }
        }    
        
        /* update composition for given chemical potential */
        memcpy(composition,mcell,CP_SIZE*sizeof(PetscScalar));
        user->rejectstage = Composition(composition,chemicalpotential,slist,user);
    
        if (slist[0] > 1) {
            /* phase capillary driving force */ 
            memset(caplflux  ,0,slist[0]*sizeof(PetscReal));
            memset(caplsource,0,slist[0]*sizeof(PetscReal));
            for (gk=0; gk<slist[0]; gk++) {
                for (gj=gk+1; gj<slist[0]; gj++) {
                    interfacekj = user->interfacelist[slist[gk+1]*user->npf + slist[gj+1]];
                    currentinterface = &user->interface[interfacekj];
                    caplflux  [gk] -= currentinterface->energy*plapl[gj];
                    caplflux  [gj] -= currentinterface->energy*plapl[gk];
                    caplsource[gk] -= currentinterface->energy*pcell[gj];
                    caplsource[gj] -= currentinterface->energy*pcell[gk];
                    for (gi=gj+1; gi<slist[0]; gi++) {
                        triplejunctionenergy = currentinterface->energy;
                        interfaceji = user->interfacelist[slist[gj+1]*user->npf + slist[gi+1]];
                        currentinterface = &user->interface[interfaceji];
                        triplejunctionenergy = triplejunctionenergy > currentinterface->energy ?
                                               triplejunctionenergy : currentinterface->energy;
                        interfaceki = user->interfacelist[slist[gk+1]*user->npf + slist[gi+1]];
                        currentinterface = &user->interface[interfaceki];
                        triplejunctionenergy = triplejunctionenergy > currentinterface->energy ?
                                               triplejunctionenergy : currentinterface->energy;
                        caplsource[gk] -= triplejunctionenergy*pcell[gj]*pcell[gi];
                        caplsource[gj] -= triplejunctionenergy*pcell[gk]*pcell[gi];
                        caplsource[gi] -= triplejunctionenergy*pcell[gk]*pcell[gj];
                    }
                }
                caplflux  [gk] *= 8.0*user->params.interfacewidth/PETSC_PI/PETSC_PI;
                caplsource[gk] *= 8.0/user->params.interfacewidth;
            }   

            /* phase chemical driving force */ 
            Chemenergy(chemsource,composition,dcell,slist,user);
            MatMulInterpolantDerivative(chemsource,pcell,slist[0]);

            /* build unconstrained RHS to calculate active set */ 
            nactivephases = 0.0;
            memset(rhs_unconstrained,0,slist[0]*sizeof(PetscReal));
            for (gk=0; gk<slist[0]; gk++) {
                for (gj=gk+1; gj<slist[0]; gj++) {
                    interfacekj = user->interfacelist[slist[gk+1]*user->npf + slist[gj+1]];
                    currentinterface = &user->interface[interfacekj];
                    rhsval = currentinterface->mobility*(  (caplflux  [gk] - caplflux  [gj])
                                                         + (caplsource[gk] - caplsource[gj])
                                                         - (chemsource[gk] - chemsource[gj])
                                                         * sqrt(interpolant[gk]*interpolant[gj]));
                    rhs_unconstrained[gk] += rhsval; rhs_unconstrained[gj] -= rhsval;
                }
                active[gk] =    (pcell[gk] >       TOL && pcell         [gk] < 1.0 - TOL)
                             || (pcell[gk] <       TOL && rhs_unconstrained[gk] > 0.0      )
                             || (pcell[gk] > 1.0 - TOL && rhs_unconstrained[gk] < 0.0      );
                if (active[gk]) nactivephases += 1.0;
            }

            /* build constrained RHS from active set*/ 
            for (gk=0; gk<slist[0]; gk++) {
                if (active[gk]) {
                    for (gj=gk+1; gj<slist[0]; gj++) {
                        if (active[gj]) {
                            interfacekj = user->interfacelist[slist[gk+1]*user->npf + slist[gj+1]];
                            currentinterface = &user->interface[interfacekj];
                            rhsval = currentinterface->mobility*(  (caplflux  [gk] - caplflux  [gj])
                                                                 + (caplsource[gk] - caplsource[gj])
                                                                 - (chemsource[gk] - chemsource[gj])
                                                                 * sqrt(interpolant[gk]*interpolant[gj]));
                            pfrhs[gk] += rhsval; pfrhs[gj] -= rhsval;
                        }
                    }
                    pfrhs[gk] /= nactivephases;
                }        
            }
        }    
            
        if (user->ndp) {
            /* calculate solute mobility matrix (volume-fixed frame of reference) */
            CompositionMobility(mobilitycv,composition,interpolant,slist,user);
  
            /* calculate avg composition rate */
            for (c=0; c<user->ndp; c++) {
                cavgdot[c] = mobilitycv[c]*dlapl[c];
                for (g =0; g<slist[0];  g++) cvgdot_phase[g] = composition[g*user->ncp+c];
                MatMulInterpolantDerivative(cvgdot_phase,pcell,slist[0]);
                for (g =0; g<slist[0];  g++) cavgdot[c] -= pfrhs[g]*cvgdot_phase[g];
            }
  
            /* calculate composition tangent wrt chemical potential */
            CompositionTangent(dcdm,composition,interpolant,slist,user);
            Invertmatrix(dmdc,dcdm);
  
            /* chemical potential RHS */
            for (cj=0; cj<user->ndp; cj++) {
                for (ci=0; ci<user->ndp; ci++) {
                    dprhs[cj] += dmdc[cj*user->ndp+ci]*cavgdot[ci];
                }
            }
        } 
    }

    /* Restore FVM residuals */
    ierr = VecRestoreArray(localX, &fdof);
    ierr = DMRestoreLocalVector(user->da_solution,&localX); CHKERRQ(ierr);
    ierr = VecRestoreArray(F, &rhs);
    ierr = VecRestoreArray(laplacian,&lap); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(user->da_solution,&laplacian); CHKERRQ(ierr);
    ierr = VecRestoreArray(user->matstate,&mat); CHKERRQ(ierr);
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
    PetscScalar    *fdof, *fadof, *mat, *offset, *offsetg;
    PetscScalar    *pcell;
    PetscScalar    *dcell;
    uint16_t       slist[AS_SIZE], glist[AS_SIZE], nlist[AS_SIZE];
    PetscScalar    *mcell;
    PetscInt       conesize, nsupp;
    const PetscInt *cone, *scells;
    PetscReal      intval;
    PetscReal      phiUP[PF_SIZE], interpolant[PF_SIZE], phiSS[PF_SIZE], compSS[CP_SIZE];
    PetscReal      chemicalpotential[DP_SIZE];
    PetscInt       localcell, cell, face, supp;
    PetscInt       g, gj, gk, c, interfacekj;
    INTERFACE      *currentinterface;
    MATERIAL       *currentmaterial;
    uint16_t       superset[AS_SIZE], setunion[AS_SIZE], setintersection[AS_SIZE], injectionL[AS_SIZE], injectionR[AS_SIZE];
    
    PetscFunctionBeginUser;    
    ierr = TSGetSolution(ts, &solution);CHKERRQ(ierr);    
    ierr = TSGetApplicationContext(ts,&user);CHKERRQ(ierr);
    
    /* Determine active phase set, update composition */
    ierr = DMGetGlobalVector(user->da_solution,&globalX); CHKERRQ(ierr);
    ierr = VecCopy(solution,globalX);
    ierr = VecGetArray(globalX,&fadof);
    ierr = VecGetArray(solution,&fdof);
    ierr = VecGetArray(user->matstate,&mat);
    for (localcell = 0; localcell < user->nlocalcells; ++localcell) {
        cell = user->localcells[localcell];

        /* get cell state */
        offset = NULL;
        ierr = DMPlexPointGlobalRef(user->da_solution,cell, fdof, &offset); CHKERRQ(ierr);
        F2IFUNC(slist,&offset[AS_OFFSET]);
        pcell = &offset[PF_OFFSET];
        dcell = &offset[DP_OFFSET];
        offsetg = NULL;
        ierr = DMPlexPointGlobalRef(user->da_solution,cell, fadof, &offsetg); CHKERRQ(ierr);
        F2IFUNC(glist,&offsetg[AS_OFFSET]);
        mcell = NULL;
        ierr = DMPlexPointGlobalRef(user->da_matstate,cell, mat,  &mcell); CHKERRQ(ierr);

        if (slist[0] > 1) {
            /* project phase fields back to Gibbs simplex */
            memcpy(phiUP,pcell,slist[0]*sizeof(PetscReal));
            SimplexProjection(pcell,phiUP,slist[0]);
            EvalInterpolant(interpolant,pcell,slist[0]);
        
            /* update active set */
            glist[0] = 0;
            for (g=0; g<slist[0];  g++) {
                if (pcell[g] > TOL) {
                    glist[++(glist[0])] = slist[g+1];
                }
            }
            I2FFUNC(&offsetg[AS_OFFSET],glist);
    
            if (user->ndp) {
                /* update material state */
                memcpy(chemicalpotential,dcell,user->ndp*sizeof(PetscScalar));
                
                /* get chemical part of diffusion potential */
                for (gk=0; gk<slist[0]; gk++) {
                    for (gj=gk+1; gj<slist[0]; gj++) {
                        interfacekj = user->interfacelist[slist[gk+1]*user->npf + slist[gj+1]];
                        currentinterface = &user->interface[interfacekj];
                        intval = sqrt(interpolant[gk]*interpolant[gj]);
                        for (c=0; c<user->ndp; c++) {
                            chemicalpotential[c] -= intval*currentinterface->potential[c];
                        }
                    }
                }
                
                /* update composition for given chemical potential */
                Composition(mcell,chemicalpotential,slist,user);
            }    
        } else {
            /* project phase fields back to Gibbs simplex */
            pcell[0] = 1.0;
    
            /* update material state */
            if (user->ndp) Composition(mcell,dcell,slist,user);
        
            /* update active set */
            glist[0] = 1; glist[1] = slist[1];
            I2FFUNC(&offsetg[AS_OFFSET],glist);
        }
    }
    ierr = VecRestoreArray(user->matstate,&mat);
    ierr = VecRestoreArray(solution,&fdof);
    ierr = VecRestoreArray(globalX,&fadof);
    ierr = DMGetLocalVector(user->da_solution,&localX); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(user->da_solution,globalX,INSERT_VALUES,localX); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->da_solution,globalX,INSERT_VALUES,localX); CHKERRQ(ierr);
    ierr = DMRestoreGlobalVector(user->da_solution,&globalX); CHKERRQ(ierr);

    /* Determine active phase super set, reorder dofs in solution, initialize new phase frac & comp */
    ierr = VecGetArray(solution,&fdof); CHKERRQ(ierr);
    ierr = VecGetArray(localX,&fadof); CHKERRQ(ierr);
    ierr = VecGetArray(user->matstate,&mat);
    for (localcell = 0; localcell < user->nlocalcells; ++localcell) {
        cell = user->localcells[localcell];

        /* get cell state */
        offsetg = NULL;
        ierr = DMPlexPointGlobalRef(user->da_solution, cell, fdof, &offsetg); CHKERRQ(ierr);
        F2IFUNC(slist,&offsetg[AS_OFFSET]);
        pcell = &offsetg[PF_OFFSET];
        mcell = NULL;
        ierr = DMPlexPointGlobalRef(user->da_matstate,cell, mat, &mcell); CHKERRQ(ierr);
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
            memcpy(&compSS[g*user->ncp],currentmaterial->c0,user->ncp*sizeof(PetscReal));
        }
        
        /* reorder dofs to new active phase superset */
        SetIntersection(setintersection,injectionL,injectionR,slist,superset);
        for (g=0; g<setintersection[0];  g++) {
            phiSS[injectionR[g]] = pcell[injectionL[g]];
            memcpy(&compSS[injectionR[g]*user->ncp],&mcell[injectionL[g]*user->ncp],user->ncp*sizeof(PetscReal));
        }
        memcpy(pcell,phiSS,superset[0]*sizeof(PetscReal));
        memcpy(mcell,compSS,superset[0]*user->ncp*sizeof(PetscReal));
        I2FFUNC(&offsetg[AS_OFFSET],superset);
    }    
    ierr = VecRestoreArray(user->matstate,&mat);
    ierr = VecRestoreArray(localX,&fadof); CHKERRQ(ierr);
    ierr = VecRestoreArray(solution,&fdof); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(user->da_solution,&localX); CHKERRQ(ierr);
    ierr = TSSetSolution(ts,solution); CHKERRQ(ierr);
    
    /* write output */
    ierr = TSGetStepNumber(ts, &user->params.step);CHKERRQ(ierr);
    if (user->params.step%user->params.outputfreq == 0) {
       Vec               Xout;
       PetscScalar       *xout, *phase, *composition;
       char              name[256];
       PetscViewer       viewer;
       
       ierr = DMGetGlobalVector(user->da_output,&Xout); CHKERRQ(ierr);
       sprintf(name, "step %d",user->params.step);
       ierr = PetscObjectSetName((PetscObject) Xout, name);CHKERRQ(ierr);
       ierr = VecZeroEntries(Xout);
       ierr = VecGetArray(solution,&fdof); CHKERRQ(ierr);
       ierr = VecGetArray(user->matstate,&mat); CHKERRQ(ierr);
       ierr = VecGetArray(Xout,&xout); CHKERRQ(ierr);
       for (localcell = 0; localcell < user->nlocalcells; ++localcell) {
           cell = user->localcells[localcell];
           /* get cell state */
           offset = NULL;
           ierr = DMPlexPointGlobalRef(user->da_solution, cell, fdof, &offset); CHKERRQ(ierr);
           F2IFUNC(slist,&offset[AS_OFFSET]);
           pcell = &offset[PF_OFFSET];
           mcell = NULL;
           ierr = DMPlexPointGlobalRef(user->da_matstate, cell, mat,  &mcell); CHKERRQ(ierr);
           offset = NULL;
           ierr = DMPlexPointGlobalRef(user->da_output,cell, xout, &offset); CHKERRQ(ierr);
           phase = &offset[0];
           composition = &offset[1];
       
           PetscReal max = -LARGE, interpolant[slist[0]];
           EvalInterpolant(interpolant,pcell,slist[0]);
           phase[0] = (PetscScalar) slist[0];
           for (g=0; g<slist[0]; g++) {
               if (interpolant[g] > max) {
                  max = interpolant[g];
                  phase[0] = (PetscScalar) slist[g+1];
               }
           }
           for (c=0; c<user->ncp; c++) {
               composition[c] = 0.0;
               for (g=0; g<slist[0]; g++)
                   composition[c] += interpolant[g]*mcell[g*user->ncp+c];        
           }
       }
       ierr = VecRestoreArray(solution,&fdof); CHKERRQ(ierr);
       ierr = VecRestoreArray(user->matstate,&mat); CHKERRQ(ierr);
       ierr = VecRestoreArray(Xout,&xout); CHKERRQ(ierr);

       sprintf(name, "%s_%d.vtu",user->params.outfile,user->params.step);
       PetscPrintf(PETSC_COMM_WORLD,"writing output to %s\n",name);
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
  ierr = TSSetMaxTime(*ts,user->params.finaltime);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(*ts,TS_EXACTFINALTIME_STEPOVER);CHKERRQ(ierr);
  ierr = TSSetPostStep(*ts,PostStep); CHKERRQ(ierr);
  ierr = TSGetAdapt(*ts,&adapt); CHKERRQ(ierr);
  ierr = TSAdaptSetType(adapt,TSADAPTDSP); CHKERRQ(ierr);
  ierr = TSAdaptSetStepLimits(adapt,user->params.mintimestep,user->params.maxtimestep); CHKERRQ(ierr);
  ierr = TSSetTolerances(*ts,user->params.abstol,NULL,user->params.reltol,NULL); CHKERRQ(ierr);
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
  VecTagger      refineTag = NULL, coarsenTag = NULL;
  /* numerical parameters */
  PetscErrorCode ierr;
  /* MPI rank */
  PetscInt       mpirank;
      
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Initialize program
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = PetscInitialize(&argc,&args,(char*)0,help); CHKERRQ(ierr);
  MPI_Comm_rank(PETSC_COMM_WORLD,&mpirank);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Initialize problem parameters
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SetUpGeometry(&ctx); CHKERRQ(ierr);
  ierr = SetUpInterface(&ctx); CHKERRQ(ierr);
  utility_init(&ctx);
  material_init(&ctx);
  ctx.rejectstage = 0;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create distributed mesh (DMPLEX) to manage 
   parallel vectors for the multi-phase field PDE
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = PetscOptionsInsertString(NULL,ctx.params.petscoptions);

  /* Create and distribute mesh */
  ierr = DMCreate(PETSC_COMM_WORLD,&ctx.da_solforest);CHKERRQ(ierr);
  ierr = DMSetType(ctx.da_solforest,ctx.dim == 2 ? DMP4EST : DMP8EST);CHKERRQ(ierr);
  ierr = DMForestSetTopology(ctx.da_solforest,"brick");CHKERRQ(ierr);
  ierr = DMForestSetInitialRefinement(ctx.da_solforest,ctx.params.initrefine);CHKERRQ(ierr);
  ierr = DMForestSetMaximumRefinement(ctx.da_solforest,ctx.params.maxnrefine);CHKERRQ(ierr);
  ierr = DMForestSetPartitionOverlap(ctx.da_solforest,1);CHKERRQ(ierr);
  ierr = DMSetFromOptions(ctx.da_solforest);CHKERRQ(ierr);
  ierr = DMSetUp(ctx.da_solforest);CHKERRQ(ierr);
  ierr = DMForestTemplate(ctx.da_solforest,PETSC_COMM_WORLD,&ctx.da_matforest);CHKERRQ(ierr);
  ierr = DMSetUp(ctx.da_matforest);CHKERRQ(ierr);
  ierr = DMConvert(ctx.da_solforest,DMPLEX,&ctx.da_solution);CHKERRQ(ierr);
  ierr = DMLocalizeCoordinates(ctx.da_solution);CHKERRQ(ierr);
  ierr = DMConvert(ctx.da_matforest,DMPLEX,&ctx.da_matstate);CHKERRQ(ierr);
  ierr = DMLocalizeCoordinates(ctx.da_matstate);CHKERRQ(ierr);
  
  /* Create finite volume discretisation for solution */
  {
    PetscFE solution_fe;
    ierr = PetscFECreateDefault(PETSC_COMM_WORLD,ctx.dim,
                                AS_SIZE+PF_SIZE+DP_SIZE,
                                PETSC_FALSE,NULL,PETSC_DEFAULT,&solution_fe);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) solution_fe, "solution");CHKERRQ(ierr);
    ierr = DMSetField(ctx.da_solution,0, NULL, (PetscObject) solution_fe);CHKERRQ(ierr);
    ierr = DMCreateDS(ctx.da_solution);CHKERRQ(ierr);
    ierr = PetscFEDestroy(&solution_fe);CHKERRQ(ierr);
    ierr = DMCopyDisc(ctx.da_solution,ctx.da_solforest);CHKERRQ(ierr);
  }
  
  /* Create finite volume discretisation for material state */
  {
    PetscFE matstate_fe;
    ierr = PetscFECreateDefault(PETSC_COMM_WORLD,ctx.dim,
                                CP_SIZE,
                                PETSC_FALSE,NULL,PETSC_DEFAULT,&matstate_fe);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) matstate_fe, "material state");CHKERRQ(ierr);
    ierr = DMSetField(ctx.da_matstate,0, NULL, (PetscObject) matstate_fe);CHKERRQ(ierr);
    ierr = DMCreateDS(ctx.da_matstate);CHKERRQ(ierr);
    ierr = PetscFEDestroy(&matstate_fe);CHKERRQ(ierr);
    ierr = DMCopyDisc(ctx.da_matstate,ctx.da_matforest);CHKERRQ(ierr);
  }

  /* Create finite volume discretisation for output */
  {
    PetscFE output_fe;
    ierr = DMClone(ctx.da_solution,&ctx.da_output);CHKERRQ(ierr);
    ierr = PetscFECreateDefault(PETSC_COMM_WORLD,ctx.dim,
                                1+ctx.ncp,
                                PETSC_FALSE,NULL,PETSC_DEFAULT,&output_fe);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) output_fe, " output");CHKERRQ(ierr);
    ierr = DMSetField(ctx.da_output, 0, NULL, (PetscObject) output_fe);CHKERRQ(ierr);
    ierr = DMCreateDS(ctx.da_output);CHKERRQ(ierr);
    ierr = PetscFEDestroy(&output_fe);CHKERRQ(ierr);
  }

  /* Pre-calculate local cells and geometry */
  {
    PetscSection lsection, gsection;
    PetscInt point, pstart, pend, ngdof, nldof, nsupp, nchild;
    
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
    ierr = DMPlexGetHeightStratum(ctx.da_solution, 1, &pstart, &pend); CHKERRQ(ierr);
    ctx.localfaces = malloc((pend-pstart)*sizeof(PetscInt)); ctx.nlocalfaces = 0;
    for (point = pstart; point < pend; ++point) {
        ierr = DMPlexGetSupportSize(ctx.da_solution, point, &nsupp);
        ierr = DMPlexGetTreeChildren(ctx.da_solution, point, &nchild, NULL);
        if (nsupp != 2 || nchild > 0) continue;
        ctx.localfaces[(ctx.nlocalfaces)++] = point;
    }
    ctx.cellgeom = malloc(ctx.ninteriorcells*sizeof(PetscReal));
    for (PetscInt cell = 0; cell < ctx.ninteriorcells; ++cell) {
        PetscReal cvolume;
        ierr = DMPlexComputeCellGeometryFVM(ctx.da_solution, cell, &cvolume, NULL, NULL);
        ctx.cellgeom[cell] = ctx.dim == 2 ? sqrt(cvolume) : cbrt(cvolume);
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
    for (localcell = 0; localcell < ctx.nlocalcells; ++localcell) {
        cell = ctx.localcells[localcell];
        ierr = DMPlexPointLocalRead(celldm,cell,cgeom,&cg);
        for (dim=0; dim<ctx.dim; ++dim) {
            ioff[dim] = (PetscInt) (cg->centroid[dim]*ctx.resolution[dim]/ctx.size[dim]);
        }
        off  = (ioff[2]*ctx.resolution[1] + ioff[1])*ctx.resolution[0] + ioff[0];
        ierr = DMSetLabelValue(ctx.da_solution, "phase", cell, ctx.phasevoxelmapping[off]);CHKERRQ(ierr);
    }
    ierr = VecRestoreArrayRead(cellgeom,&cgeom);CHKERRQ(ierr);
  }

  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Set up problem
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMCreateGlobalVector(ctx.da_solution, &solution);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) solution, "solution");CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(ctx.da_matstate, &ctx.matstate);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) ctx.matstate, "material state");CHKERRQ(ierr);
  ierr = SetUpProblem(solution,&ctx);CHKERRQ(ierr);
  
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Set up and perform time integration 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  for (PetscInt initcoarseniter = 0; initcoarseniter < ctx.params.initcoarsen; ++initcoarseniter){
      PetscPrintf(PETSC_COMM_WORLD,"...initial refinement/coarsening (%d)...\n",initcoarseniter);
      /* Adapt mesh */
      PetscInt cell;
      PetscScalar *fdof, *offset;
      DMLabel adaptlabel;
      DM postsolforest, postmatforest;
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

      postsolforest = NULL;
      ierr = DMAdaptLabel(ctx.da_solforest,adaptlabel,&postsolforest);CHKERRQ(ierr);
      if (postsolforest != NULL) {
          ierr = DMAdaptLabel(ctx.da_matforest,adaptlabel,&postmatforest);CHKERRQ(ierr);
          
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

          ierr = DMCreateGlobalVector(postmatforest, &postvec);CHKERRQ(ierr);
          ierr = DMForestTransferVec(ctx.da_matforest, ctx.matstate, postmatforest, postvec, PETSC_TRUE, 0.0);CHKERRQ(ierr);
          ierr = VecDestroy(&ctx.matstate);CHKERRQ(ierr);
          ierr = DMDestroy(&ctx.da_matforest);CHKERRQ(ierr);
          ierr = DMDestroy(&ctx.da_matstate);CHKERRQ(ierr);
          ctx.da_matforest = postmatforest;
          ierr = DMForestSetAdaptivityForest(ctx.da_matforest,NULL);CHKERRQ(ierr);
          ierr = DMConvert(ctx.da_matforest,DMPLEX,&ctx.da_matstate);CHKERRQ(ierr);
          ierr = DMCopyDisc(ctx.da_matforest,ctx.da_matstate);CHKERRQ(ierr);
          ierr = DMCreateGlobalVector(ctx.da_matstate,&ctx.matstate);CHKERRQ(ierr);
          ierr = VecCopy(postvec,ctx.matstate);CHKERRQ(ierr);
          ierr = VecDestroy(&postvec);CHKERRQ(ierr);

          {
            PetscSection lsection, gsection;
            PetscInt point, pstart, pend, ngdof, nldof, nsupp, nchild;
            
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
                PetscReal cvolume;
                ierr = DMPlexComputeCellGeometryFVM(ctx.da_solution, cell, &cvolume, NULL, NULL);
                ctx.cellgeom[cell] = ctx.dim == 2 ? sqrt(cvolume) : cbrt(cvolume);
            }
          }
          
          {
            PetscFE output_fe;
            ierr = DMDestroy(&ctx.da_output);CHKERRQ(ierr);
            ierr = DMClone(ctx.da_solution,&ctx.da_output);CHKERRQ(ierr);
            ierr = PetscFECreateDefault(PETSC_COMM_WORLD,ctx.dim,
                                        1+ctx.ncp,
                                        PETSC_FALSE,NULL,PETSC_DEFAULT,&output_fe);CHKERRQ(ierr);
            ierr = PetscObjectSetName((PetscObject) output_fe, " output");CHKERRQ(ierr);
            ierr = DMSetField(ctx.da_output, 0, NULL, (PetscObject) output_fe);CHKERRQ(ierr);
            ierr = PetscFEDestroy(&output_fe);CHKERRQ(ierr);
          }
      }
      ierr = DMLabelDestroy(&adaptlabel);CHKERRQ(ierr);
  }    
  ierr = InitializeTS(ctx.da_solution, &ctx, &ts);CHKERRQ(ierr);

  PetscReal currenttime = 0.0;
  PetscInt  nsteps = 0;
  for (;currenttime < ctx.params.finaltime;) {
      ierr = TSSetStepNumber(ts,nsteps);CHKERRQ(ierr);
      ierr = TSSetTime(ts,currenttime);CHKERRQ(ierr);
      ierr = TSSetTimeStep(ts,ctx.params.timestep);CHKERRQ(ierr);
      ierr = TSSetMaxSteps(ts,nsteps+ctx.params.amrinterval);CHKERRQ(ierr);
      ierr = TSSolve(ts,solution);CHKERRQ(ierr);
      ierr = TSGetSolveTime(ts,&currenttime);CHKERRQ(ierr);
      ierr = TSGetStepNumber(ts,&nsteps);CHKERRQ(ierr);
      ierr = TSGetTimeStep(ts,&ctx.params.timestep);CHKERRQ(ierr);
      {
          PetscPrintf(PETSC_COMM_WORLD,"...remeshing...\n");
          /* Adapt mesh */
          PetscInt cell;
          PetscScalar *fdof, *offset;
          DMLabel adaptlabel;
          DM postsolforest, postmatforest;
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

          postsolforest = NULL;
          ierr = DMAdaptLabel(ctx.da_solforest,adaptlabel,&postsolforest);CHKERRQ(ierr);
          if (postsolforest != NULL) {
              ierr = DMAdaptLabel(ctx.da_matforest,adaptlabel,&postmatforest);CHKERRQ(ierr);
              
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

              ierr = DMCreateGlobalVector(postmatforest, &postvec);CHKERRQ(ierr);
              ierr = DMForestTransferVec(ctx.da_matforest, ctx.matstate, postmatforest, postvec, PETSC_TRUE, 0.0);CHKERRQ(ierr);
              ierr = VecDestroy(&ctx.matstate);CHKERRQ(ierr);
              ierr = DMDestroy(&ctx.da_matforest);CHKERRQ(ierr);
              ierr = DMDestroy(&ctx.da_matstate);CHKERRQ(ierr);
              ctx.da_matforest = postmatforest;
              ierr = DMForestSetAdaptivityForest(ctx.da_matforest,NULL);CHKERRQ(ierr);
              ierr = DMConvert(ctx.da_matforest,DMPLEX,&ctx.da_matstate);CHKERRQ(ierr);
              ierr = DMCopyDisc(ctx.da_matforest,ctx.da_matstate);CHKERRQ(ierr);
              ierr = DMCreateGlobalVector(ctx.da_matstate,&ctx.matstate);CHKERRQ(ierr);
              ierr = VecCopy(postvec,ctx.matstate);CHKERRQ(ierr);
              ierr = VecDestroy(&postvec);CHKERRQ(ierr);

              {
                PetscSection lsection, gsection;
                PetscInt point, pstart, pend, ngdof, nldof, nsupp, nchild;
            
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
                    PetscReal cvolume;
                    ierr = DMPlexComputeCellGeometryFVM(ctx.da_solution, cell, &cvolume, NULL, NULL);
                    ctx.cellgeom[cell] = ctx.dim == 2 ? sqrt(cvolume) : cbrt(cvolume);
                }
              }
          
              {
                PetscFE output_fe;
                ierr = DMDestroy(&ctx.da_output);CHKERRQ(ierr);
                ierr = DMClone(ctx.da_solution,&ctx.da_output);CHKERRQ(ierr);
                ierr = PetscFECreateDefault(PETSC_COMM_WORLD,ctx.dim,
                                            1+ctx.ncp,
                                            PETSC_FALSE,NULL,PETSC_DEFAULT,&output_fe);CHKERRQ(ierr);
                ierr = PetscObjectSetName((PetscObject) output_fe, " output");CHKERRQ(ierr);
                ierr = DMSetField(ctx.da_output, 0, NULL, (PetscObject) output_fe);CHKERRQ(ierr);
                ierr = PetscFEDestroy(&output_fe);CHKERRQ(ierr);
              }
              ierr = TSDestroy(&ts);CHKERRQ(ierr);
              ierr = InitializeTS(ctx.da_solution, &ctx, &ts);CHKERRQ(ierr);
          }
          ierr = DMLabelDestroy(&adaptlabel);CHKERRQ(ierr);
      }    
  }

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Free work space.  All PETSc objects should be destroyed when they
   are no longer needed.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = VecTaggerDestroy(&refineTag);CHKERRQ(ierr);
  ierr = VecTaggerDestroy(&coarsenTag);CHKERRQ(ierr);
  ierr = VecDestroy(&solution);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx.da_solution);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}
