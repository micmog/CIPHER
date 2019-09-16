
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
#include <sys/time.h>

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
    PFIELD            *pcell, *pcellL, *pcellR, *plapl, *plaplL, *plaplR, *pfrhs, fluxp;
    DFIELD            *dcell, *dcellL, *dcellR, *dlapl, *dlaplL, *dlaplR, *dprhs, fluxd;
    F2I               *slist, *slistL, *slistR;
    STATE             *mcell;
    const PetscInt    *scells;
    PetscReal         ffactor, cfactor, deltaL, deltaR;
    PetscInt          localcell, cell, localface, face; 
    PetscReal         interpolant[MAXAP], caplflux[MAXAP], caplsource[MAXAP];
    PetscReal         chemsource[MAXAP], rhs_unconstrained[MAXAP], cvgdot_phase[MAXAP];
    PetscBool         active[MAXAP];
    PetscReal         chemicalpotential[MAXCP], composition[CP_SIZE];
    PetscReal         mobilitycv[MAXCP], cavgdot[MAXCP];
    PetscReal         dcdm[MAXCP*MAXCP], dmdc[MAXCP*MAXCP];
    uint16_t          setintersect[MAXIP], injectionL[MAXIP], injectionR[MAXIP];
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
        slistL = (F2I    *) &offset[AS_OFFSET];
        pcellL = (PFIELD *) &offset[PF_OFFSET];
        dcellL = (DFIELD *) &offset[DP_OFFSET];
        offset = NULL;
        ierr = DMPlexPointLocalRef(user->da_solution, scells[1], fdof, &offset); CHKERRQ(ierr);
        slistR = (F2I    *) &offset[AS_OFFSET];
        pcellR = (PFIELD *) &offset[PF_OFFSET];
        dcellR = (DFIELD *) &offset[DP_OFFSET];
        if (slistL->i[0] <= 1 && slistR->i[0] <= 1 && user->nc <= 1) continue;

        /* get geometric data */
        deltaL = user->cellgeom[scells[0]]; deltaR = user->cellgeom[scells[0]];
        ffactor = 2.0*(deltaL < deltaR ? deltaL*deltaL : deltaR*deltaR)/(deltaL + deltaR);
        
        /* get common active phases */
        if (slistL->i[0] > 1 || slistR->i[0] > 1) {
            SetIntersection(setintersect,injectionL,injectionR,slistL->i,slistR->i);
            for (g=0; g<setintersect[0]; g++) 
                fluxp.p[g] = (pcellR->p[injectionR[g]] - pcellL->p[injectionL[g]]);
        }
        for (c=0; c<user->nc-1; c++) fluxd.d[c] = (dcellR->d[c] - dcellL->d[c]);
        
        {
            offset = NULL;
            ierr = DMPlexPointLocalRef(user->da_solution, scells[0], lap,  &offset); CHKERRQ(ierr);
            plaplL = (PFIELD *) &offset[PF_OFFSET];
            dlaplL = (DFIELD *) &offset[DP_OFFSET];
            cfactor = ffactor/FastPow(deltaL,3);
            if (slistL->i[0] > 1) {
                for (g=0; g<setintersect[0]; g++)
                    plaplL->p[injectionL[g]] += cfactor*fluxp.p[g];
            }
            for (c=0; c<user->nc-1; c++) dlaplL->d[c] += cfactor*fluxd.d[c];
        }
        {
            offset = NULL;
            ierr = DMPlexPointLocalRef(user->da_solution, scells[1], lap,  &offset); CHKERRQ(ierr);
            plaplR = (PFIELD *) &offset[PF_OFFSET];
            dlaplR = (DFIELD *) &offset[DP_OFFSET];
            cfactor = ffactor/FastPow(deltaR,3);
            if (slistR->i[0] > 1) {
                for (g=0; g<setintersect[0]; g++)
                    plaplR->p[injectionR[g]] -= cfactor*fluxp.p[g];
            }
            for (c=0; c<user->nc-1; c++) dlaplR->d[c] -= cfactor*fluxd.d[c];
        }
    }    
    /* Loop over cells and compute local contribution to the RHS */
    for (localcell = 0; localcell < user->nlocalcells; ++localcell) {
        cell = user->localcells[localcell];

        /* get fields */
        offset = NULL;
        ierr = DMPlexPointLocalRef(user->da_solution, cell, fdof, &offset); CHKERRQ(ierr);
        slist = (F2I    *) &offset[AS_OFFSET];
        pcell = (PFIELD *) &offset[PF_OFFSET];
        dcell = (DFIELD *) &offset[DP_OFFSET];
        offset = NULL;
        ierr = DMPlexPointLocalRef(user->da_solution, cell, lap,  &offset); CHKERRQ(ierr);
        plapl = (PFIELD *) &offset[PF_OFFSET];
        dlapl = (DFIELD *) &offset[DP_OFFSET];
        offset = NULL;
        ierr = DMPlexPointGlobalRef(user->da_solution,cell, rhs,  &offset); CHKERRQ(ierr);
        pfrhs = (PFIELD *) &offset[PF_OFFSET];
        dprhs = (DFIELD *) &offset[DP_OFFSET];
        mcell = NULL;
        ierr = DMPlexPointGlobalRef(user->da_matstate,cell, mat,  &mcell); CHKERRQ(ierr);
        
        /* calculate phase interpolants */
        EvalInterpolant(interpolant,pcell->p,slist->i[0]);
                
        /* update material state */
        memcpy(chemicalpotential,dcell->d,(user->nc-1)*sizeof(PetscScalar));
        
        /* get chemical part of diffusion potential */
        for (gk=0; gk<slist->i[0]; gk++) {
            for (gj=gk+1; gj<slist->i[0]; gj++) {
                interfacekj = (PetscInt) user->interfacelist[slist->i[gk+1]*user->np + slist->i[gj+1]];
                currentinterface = &user->interface[interfacekj];
                intval = sqrt(interpolant[gk]*interpolant[gj]);
                for (c=0; c<user->nc; c++) {
                    chemicalpotential[c] -= intval*currentinterface->potential[c];
                }
            }
        }    
        
        /* update composition for given chemical potential */
        memcpy(composition,mcell->c,CP_SIZE*sizeof(PetscScalar));
        user->rejectstage = Composition(composition,chemicalpotential,slist->i,user);
    
        if (slist->i[0] > 1) {
            /* phase capillary driving force */ 
            memset(caplflux  ,0,slist->i[0]*sizeof(PetscReal));
            memset(caplsource,0,slist->i[0]*sizeof(PetscReal));
            for (gk=0; gk<slist->i[0]; gk++) {
                for (gj=gk+1; gj<slist->i[0]; gj++) {
                    interfacekj = (PetscInt) user->interfacelist[slist->i[gk+1]*user->np + slist->i[gj+1]];
                    currentinterface = &user->interface[interfacekj];
                    caplflux  [gk] -= currentinterface->energy*plapl->p[gj];
                    caplflux  [gj] -= currentinterface->energy*plapl->p[gk];
                    caplsource[gk] -= currentinterface->energy*pcell->p[gj];
                    caplsource[gj] -= currentinterface->energy*pcell->p[gk];
                    for (gi=gj+1; gi<slist->i[0]; gi++) {
                        triplejunctionenergy = currentinterface->energy;
                        interfaceji = (PetscInt) user->interfacelist[slist->i[gj+1]*user->np + slist->i[gi+1]];
                        currentinterface = &user->interface[interfaceji];
                        triplejunctionenergy = triplejunctionenergy > currentinterface->energy ?
                                               triplejunctionenergy : currentinterface->energy;
                        interfaceki = (PetscInt) user->interfacelist[slist->i[gk+1]*user->np + slist->i[gi+1]];
                        currentinterface = &user->interface[interfaceki];
                        triplejunctionenergy = triplejunctionenergy > currentinterface->energy ?
                                               triplejunctionenergy : currentinterface->energy;
                        caplsource[gk] -= triplejunctionenergy*pcell->p[gj]*pcell->p[gi];
                        caplsource[gj] -= triplejunctionenergy*pcell->p[gk]*pcell->p[gi];
                        caplsource[gi] -= triplejunctionenergy*pcell->p[gk]*pcell->p[gj];
                    }
                }
                caplflux  [gk] *= 8.0*user->params.interfacewidth/PETSC_PI/PETSC_PI;
                caplsource[gk] *= 8.0/user->params.interfacewidth;
            }   

            /* phase chemical driving force */ 
            Chemenergy(chemsource,composition,dcell->d,slist->i,user);
            MatMulInterpolantDerivative(chemsource,pcell->p,slist->i[0]);

            /* build unconstrained RHS to calculate active set */ 
            nactivephases = 0.0;
            memset(rhs_unconstrained,0,slist->i[0]*sizeof(PetscReal));
            for (gk=0; gk<slist->i[0]; gk++) {
                for (gj=gk+1; gj<slist->i[0]; gj++) {
                    interfacekj = (PetscInt) user->interfacelist[slist->i[gk+1]*user->np + slist->i[gj+1]];
                    currentinterface = &user->interface[interfacekj];
                    rhsval = currentinterface->mobility*(  (caplflux  [gk] - caplflux  [gj])
                                                         + (caplsource[gk] - caplsource[gj])
                                                         + (chemsource[gj] - chemsource[gk])
                                                         * sqrt(interpolant[gk]*interpolant[gj]));
                    rhs_unconstrained[gk] += rhsval; rhs_unconstrained[gj] -= rhsval;
                }
                active[gk] =    (pcell->p[gk] >       TOL && pcell->p         [gk] < 1.0 - TOL)
                             || (pcell->p[gk] <       TOL && rhs_unconstrained[gk] > 0.0      )
                             || (pcell->p[gk] > 1.0 - TOL && rhs_unconstrained[gk] < 0.0      );
                if (active[gk]) nactivephases += 1.0;
            }

            /* build constrained RHS from active set*/ 
            for (gk=0; gk<slist->i[0]; gk++) {
                if (active[gk]) {
                    for (gj=gk+1; gj<slist->i[0]; gj++) {
                        if (active[gj]) {
                            interfacekj = (PetscInt) user->interfacelist[slist->i[gk+1]*user->np + slist->i[gj+1]];
                            currentinterface = &user->interface[interfacekj];
                            rhsval = currentinterface->mobility*(  (caplflux  [gk] - caplflux  [gj])
                                                                 + (caplsource[gk] - caplsource[gj])
                                                                 + (chemsource[gj] - chemsource[gk])
                                                                 * sqrt(interpolant[gk]*interpolant[gj]));
                            pfrhs->p[gk] += rhsval; pfrhs->p[gj] -= rhsval;
                        }
                    }
                    pfrhs->p[gk] /= nactivephases;
                }        
            }
        }    
            
        if (user->nc > 1) {
            /* calculate solute mobility matrix (volume-fixed frame of reference) */
            CompositionMobility(mobilitycv,composition,interpolant,slist->i,user);
  
            /* calculate avg composition rate */
            for (c=0; c<user->nc-1; c++) {
                cavgdot[c] = mobilitycv[c]*dlapl->d[c];
                for (g =0; g<slist->i[0];  g++) cvgdot_phase[g] = composition[g*user->nc+c];
                MatMulInterpolantDerivative(cvgdot_phase,pcell->p,slist->i[0]);
                for (g =0; g<slist->i[0];  g++) cavgdot[c] -= pfrhs->p[g]*cvgdot_phase[g];
            }
  
            /* calculate composition tangent wrt chemical potential */
            CompositionTangent(dcdm,composition,interpolant,slist->i,user);
            Invertmatrix(dmdc,dcdm);
  
            /* chemical potential RHS */
            for (cj=0; cj<user->nc-1; cj++) {
                for (ci=0; ci<user->nc-1; ci++) {
                    dprhs->d[cj] += dmdc[cj*(user->nc-1)+ci]*cavgdot[ci];
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
    PetscScalar    *fdof, *fadof, *mat, *offset;
    PFIELD         *pcell;
    DFIELD         *dcell;
    F2I            *slist, *glist, *nlist;
    STATE          *mcell;
    PetscInt       conesize, nsupp;
    const PetscInt *cone, *scells;
    PetscReal      intval;
    PetscReal      phiUP[MAXAP], interpolant[MAXAP], phiSS[MAXAP], compSS[MAXAP*MAXCP];
    PetscReal      interfacepotential[MAXCP], chemicalpotential[MAXCP];
    PetscInt       localcell, cell, face, supp;
    PetscInt       g, gj, gk, c, interfacekj;
    INTERFACE      *currentinterface;
    MATERIAL       *currentmaterial;
    uint16_t       superset[MAXIP], setunion[MAXIP], setintersection[MAXIP], injectionL[MAXIP], injectionR[MAXIP];
    
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
        slist = (F2I    *) &offset[AS_OFFSET];
        pcell = (PFIELD *) &offset[PF_OFFSET];
        dcell = (DFIELD *) &offset[DP_OFFSET];
        offset = NULL;
        ierr = DMPlexPointGlobalRef(user->da_solution,cell, fadof, &offset); CHKERRQ(ierr);
        glist = (F2I    *) &offset[AS_OFFSET];
        mcell = NULL;
        ierr = DMPlexPointGlobalRef(user->da_matstate,cell, mat,  &mcell); CHKERRQ(ierr);

        if (slist->i[0] > 1) {
            /* project phase fields back to Gibbs simplex */
            memcpy(phiUP,pcell->p,slist->i[0]*sizeof(PetscReal));
            SimplexProjection(pcell->p,phiUP,slist->i[0]);
            EvalInterpolant(interpolant,pcell->p,slist->i[0]);
    
            /* update material state */
            memset(interfacepotential,0,user->nc*sizeof(PetscScalar));
            memcpy(chemicalpotential,dcell->d,(user->nc-1)*sizeof(PetscScalar));
            /* get chemical part of diffusion potential */
            for (gk=0; gk<slist->i[0]; gk++) {
                for (gj=gk+1; gj<slist->i[0]; gj++) {
                    interfacekj = (PetscInt) user->interfacelist[slist->i[gk+1]*user->np + slist->i[gj+1]];
                    currentinterface = &user->interface[interfacekj];
                    intval = sqrt(interpolant[gk]*interpolant[gj]);
                    for (c=0; c<user->nc; c++) {
                        interfacepotential[c] += intval*currentinterface->potential[c];
                    }
                }
            }    
            for (c=0; c<user->nc-1; c++) 
                chemicalpotential[c] -= (interfacepotential[c] - interfacepotential[user->nc-1]);
        
            /* update composition for given chemical potential */
            Composition(mcell->c,chemicalpotential,slist->i,user);
        
            /* update active set */
            glist->i[0] = 0;
            for (g=0; g<slist->i[0];  g++) {
                if (interpolant[g] > TOL) {
                    glist->i[++(glist->i[0])] = slist->i[g+1];
                }
            }
        } else {
            /* project phase fields back to Gibbs simplex */
            pcell->p[0] = 1.0;
    
            /* update material state */
            Composition(mcell->c,dcell->d,slist->i,user);
        
            /* update active set */
            glist->i[0] = 1; glist->i[1] = slist->i[1];
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
        offset = NULL;
        ierr = DMPlexPointGlobalRef(user->da_solution, cell, fdof, &offset); CHKERRQ(ierr);
        slist = (F2I    *) &offset[AS_OFFSET];
        pcell = (PFIELD *) &offset[PF_OFFSET];
        mcell = NULL;
        ierr = DMPlexPointGlobalRef(user->da_matstate,cell, mat,  &mcell); CHKERRQ(ierr);

        /* add union of neighbouring cells */
        memcpy(superset,slist->i,MAXIP*sizeof(uint16_t));
        ierr = DMPlexGetConeSize(user->da_solution,cell,&conesize);
        ierr = DMPlexGetCone    (user->da_solution,cell,&cone    );
        for (face=0; face<conesize; face++) {
            ierr = DMPlexGetSupportSize(user->da_solution, cone[face], &nsupp);
            ierr = DMPlexGetSupport(user->da_solution, cone[face], &scells);
            for (supp=0; supp<nsupp; supp++) {
                if (scells[supp] != cell && scells[supp] < user->ninteriorcells) {
                    offset = NULL;
                    ierr = DMPlexPointLocalRef(user->da_solution, scells[supp], fadof, &offset); CHKERRQ(ierr);
                    nlist = (F2I    *) &offset[AS_OFFSET];
                    SetUnion(setunion,injectionL,injectionR,superset,nlist->i);
                    memcpy(superset,setunion,MAXIP*sizeof(uint16_t));
                }
            }
        }

        /* initialize new phase states */
        for (g=0; g<superset[0];  g++) {
            phiSS[g] = 0.0;
            currentmaterial = &user->material[user->phasematerialmapping[superset[g+1]]];
            memcpy(&compSS[g*user->nc],currentmaterial->c0,user->nc*sizeof(PetscReal));
        }
        
        /* reorder dofs to new active phase superset */
        SetIntersection(setintersection,injectionL,injectionR,slist->i,superset);
        for (g=0; g<setintersection[0];  g++) {
            phiSS[injectionR[g]] = pcell->p[injectionL[g]];
            memcpy(&compSS[injectionR[g]*user->nc],&mcell->c[injectionL[g]*user->nc],user->nc*sizeof(PetscReal));
        }
        memcpy(pcell->p,phiSS,superset[0]*sizeof(PetscReal));
        memcpy(mcell->c,compSS,superset[0]*user->nc*sizeof(PetscReal));
        memcpy(slist->i,superset,MAXIP*sizeof(uint16_t));
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
       PetscScalar       *xout;
       char              name[256];
       PetscViewer       viewer;
       O_DOFS            *ocell;
       
       ierr = DMGetGlobalVector(user->da_output,&Xout); CHKERRQ(ierr);
       ierr = PetscObjectSetName((PetscObject) Xout, "output");CHKERRQ(ierr);
       ierr = VecGetArray(solution,&fdof); CHKERRQ(ierr);
       ierr = VecGetArray(user->matstate,&mat); CHKERRQ(ierr);
       ierr = VecGetArray(Xout,&xout); CHKERRQ(ierr);
       for (localcell = 0; localcell < user->nlocalcells; ++localcell) {
           cell = user->localcells[localcell];
           /* get cell state */
           offset = NULL;
           ierr = DMPlexPointGlobalRef(user->da_solution,cell, fdof, &offset); CHKERRQ(ierr);
           slist = (F2I    *) &offset[AS_OFFSET];
           pcell = (PFIELD *) &offset[PF_OFFSET];
           dcell = (DFIELD *) &offset[DP_OFFSET];
           mcell = NULL;
           ierr = DMPlexPointGlobalRef(user->da_matstate,cell, mat,  &mcell); CHKERRQ(ierr);
           ocell = NULL;
           ierr = DMPlexPointGlobalRef(user->da_output  ,cell, xout, &ocell); CHKERRQ(ierr);
       
           PetscReal max = -LARGE, interpolant[slist->i[0]];
           EvalInterpolant(interpolant,pcell->p,slist->i[0]);
           for (g=0; g<slist->i[0]; g++) {
               if (interpolant[g] > max) {
                  max = interpolant[g];
                  ocell->p = (PetscReal) slist->i[g+1];
               }
           }
           for (c=0; c<user->nc; c++) {
               ocell->c[c] = 0.0;
               for (g=0; g<slist->i[0]; g++)
                   ocell->c[c] += interpolant[g]*mcell->c[g*user->nc+c];        
           }
       }
       ierr = VecRestoreArray(solution,&fdof); CHKERRQ(ierr);
       ierr = VecRestoreArray(user->matstate,&mat); CHKERRQ(ierr);
       ierr = VecRestoreArray(Xout,&xout); CHKERRQ(ierr);

       sprintf(name, "%s_%d.vtu",user->params.outfile,user->params.step);
       PetscPrintf(PETSC_COMM_WORLD,"writing output to %s\n",name);
       ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer);CHKERRQ(ierr);
       ierr = PetscViewerSetType(viewer, PETSCVIEWERVTK);CHKERRQ(ierr);
       ierr = PetscViewerFileSetName(viewer, name);CHKERRQ(ierr);
       ierr = VecView(Xout,viewer);
       ierr = DMRestoreGlobalVector(user->da_output,&Xout); CHKERRQ(ierr);
       ierr = PetscViewerDestroy(&viewer);
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
  TS             ts;
  VecTagger      refineTag = NULL, coarsenTag = NULL;
  /* numerical parameters */
  PetscErrorCode ierr;
  /* MPI rank */
  PetscInt       mpirank;
      
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Initialize program
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  assert(MAXIP == MAXFP*sizeof(PetscReal)/sizeof(uint16_t));
  assert(MAXAP <= MAXFP*sizeof(PetscReal)/sizeof(uint16_t));
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
  ierr = DMSetType(ctx.da_solforest,DMP8EST);CHKERRQ(ierr);
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
  
  /* Add ghost label */
  {
    DM da_ghosted = NULL;

    ierr = DMPlexConstructGhostCells(ctx.da_solution, NULL, NULL, &da_ghosted);CHKERRQ(ierr);
    if (da_ghosted) {
      ierr = DMDestroy(&ctx.da_solution);CHKERRQ(ierr);
      ctx.da_solution  = da_ghosted;
    }  
  }
  
  /* Pre-calculate ghost cell index */
  {
    DMLabel ghostlabel;
    PetscInt point, pstart, pend, ghost, nsupp, nchild;
    
    ierr = DMGetLabel(ctx.da_solution, "ghost", &ghostlabel);
    ierr = DMPlexGetHeightStratum(ctx.da_solution, 0, &pstart, &pend); CHKERRQ(ierr);
    ierr = DMPlexGetHybridBounds(ctx.da_solution,&ctx.ninteriorcells,NULL,NULL,NULL);
    ctx.localcells = malloc(ctx.ninteriorcells*sizeof(PetscInt)); ctx.nlocalcells = 0;
    for (point = pstart; point < ctx.ninteriorcells; ++point) {
        ierr = DMLabelGetValue(ghostlabel, point, &ghost);
        if (ghost >= 0) continue; 
        ctx.localcells[(ctx.nlocalcells)++] = point;
    }
    ierr = DMPlexGetHeightStratum(ctx.da_solution, 1, &pstart, &pend); CHKERRQ(ierr);
    ctx.localfaces = malloc((pend-pstart)*sizeof(PetscInt)); ctx.nlocalfaces = 0;
    for (point = pstart; point < pend; ++point) {
        ierr = DMLabelGetValue(ghostlabel, point, &ghost);
        ierr = DMPlexGetSupportSize(ctx.da_solution, point, &nsupp);
        ierr = DMPlexGetTreeChildren(ctx.da_solution, point, &nchild, NULL);
        if (ghost >= 0 || nsupp > 2 || nchild > 0) continue;
        ctx.localfaces[(ctx.nlocalfaces)++] = point;
    }
  }
  
  /* Add phase label */
  {
    DM celldm;
    Vec cellgeom;
    const PetscScalar *cgeom;
    PetscFVCellGeom *cg;
    PetscInt cell, localcell, off, ioff, joff, koff;

    ierr = DMPlexGetDataFVM(ctx.da_solution, NULL, &cellgeom, NULL, NULL);
    ierr = VecGetDM(cellgeom,&celldm);
    ierr = VecGetArrayRead(cellgeom,&cgeom);CHKERRQ(ierr);
    ierr = DMCreateLabel(ctx.da_solution, "phase");CHKERRQ(ierr);
    for (localcell = 0; localcell < ctx.nlocalcells; ++localcell) {
        cell = ctx.localcells[localcell];
        ierr = DMPlexPointLocalRead(celldm,cell,cgeom,&cg);
        ioff = (PetscInt) (cg->centroid[0]*ctx.resolution[0]/ctx.size[0]);
        joff = (PetscInt) (cg->centroid[1]*ctx.resolution[1]/ctx.size[1]);
        koff = (PetscInt) (cg->centroid[2]*ctx.resolution[2]/ctx.size[2]);
        off  = (koff*ctx.resolution[1] + joff)*ctx.resolution[0] + ioff;
        ierr = DMSetLabelValue(ctx.da_solution, "phase", cell, ctx.phasevoxelmapping[off]);CHKERRQ(ierr);
    }
    ierr = VecRestoreArrayRead(cellgeom,&cgeom);CHKERRQ(ierr);
  }
  
  /* Create finite volume discretisation for solution */
  PetscFE solution_fe;
  ierr = PetscFECreateDefault(PETSC_COMM_WORLD,3,
                              AS_SIZE+PF_SIZE+DP_SIZE,
                              PETSC_FALSE,NULL,PETSC_DEFAULT,&solution_fe);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) solution_fe, "solution");CHKERRQ(ierr);
  ierr = DMSetField(ctx.da_solution,0, NULL, (PetscObject) solution_fe);CHKERRQ(ierr);
  ierr = DMCreateDS(ctx.da_solution);CHKERRQ(ierr);
  ierr = DMPlexCreateClosureIndex(ctx.da_solution, NULL);
  ierr = PetscFEDestroy(&solution_fe);CHKERRQ(ierr);
  ierr = DMCopyDisc(ctx.da_solution,ctx.da_solforest);CHKERRQ(ierr);
  
  /* Create finite volume discretisation for material state */
  PetscFE matstate_fe;
  ierr = PetscFECreateDefault(PETSC_COMM_WORLD,3,
                              CP_SIZE,
                              PETSC_FALSE,NULL,PETSC_DEFAULT,&matstate_fe);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) matstate_fe, "material state");CHKERRQ(ierr);
  ierr = DMSetField(ctx.da_matstate,0, NULL, (PetscObject) matstate_fe);CHKERRQ(ierr);
  ierr = DMCreateDS(ctx.da_matstate);CHKERRQ(ierr);
  ierr = DMPlexCreateClosureIndex(ctx.da_matstate, NULL);
  ierr = PetscFEDestroy(&matstate_fe);CHKERRQ(ierr);
  ierr = DMCopyDisc(ctx.da_matstate,ctx.da_matforest);CHKERRQ(ierr);

  PetscFE output_fe;
  ierr = DMClone(ctx.da_solution,&ctx.da_output);
  ierr = PetscFECreateDefault(PETSC_COMM_WORLD,3,
                              1+DP_SIZE,
                              PETSC_FALSE,NULL,PETSC_DEFAULT,&output_fe);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) output_fe, "output");CHKERRQ(ierr);
  ierr = DMSetField(ctx.da_output,0, NULL, (PetscObject) output_fe);CHKERRQ(ierr);
  ierr = DMCreateDS(ctx.da_output);CHKERRQ(ierr);
  ierr = DMPlexCreateClosureIndex(ctx.da_output, NULL);
  ierr = PetscFEDestroy(&output_fe);CHKERRQ(ierr);

  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Set up problem
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMCreateGlobalVector(ctx.da_solution, &solution);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) solution, "solution");CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(ctx.da_matstate, &ctx.matstate);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) solution, "material state");CHKERRQ(ierr);
  ctx.cellgeom = malloc(ctx.ninteriorcells*sizeof(PetscReal));
  for (PetscInt cell = 0; cell < ctx.ninteriorcells; ++cell) {
      PetscReal cvolume;
      ierr = DMPlexComputeCellGeometryFVM(ctx.da_solution, cell, &cvolume, NULL, NULL);
      ctx.cellgeom[cell] = cbrt(cvolume);
  }
  ierr = SetUpProblem(solution,&ctx);CHKERRQ(ierr);
  
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Set up and perform time integration 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscReal currenttime = 0.0;
  PetscInt  nsteps = 0;
  for (;currenttime < ctx.params.finaltime;) {
      {
          PetscPrintf(PETSC_COMM_WORLD,"...remeshing...\n");
          /* Adapt mesh */
          PetscBool success;
          PetscInt cell;
          PetscScalar *fdof, *offset;
          F2I *slist;
          DMLabel adaptlabel;
          DM postsolforest, postmatforest;
          Vec lsolution, postvec;
          PetscSection lsection, gsection;
          PetscInt point, pstart, pend, ngdof, nldof, nsupp, nchild;
      
          ierr = DMLabelCreate(PETSC_COMM_SELF,"adapt",&adaptlabel);CHKERRQ(ierr);
          ierr = DMLabelSetDefaultValue(adaptlabel,DM_ADAPT_COARSEN);CHKERRQ(ierr);

          ierr = DMGetLocalVector(ctx.da_solution,&lsolution);CHKERRQ(ierr);
          ierr = DMGlobalToLocalBegin(ctx.da_solution,solution,INSERT_VALUES,lsolution); CHKERRQ(ierr);
          ierr = DMGlobalToLocalEnd(ctx.da_solution,solution,INSERT_VALUES,lsolution); CHKERRQ(ierr);
          ierr = VecGetArray(lsolution, &fdof); CHKERRQ(ierr);
          for (cell = 0; cell < ctx.ninteriorcells; ++cell) {
              ierr = DMPlexPointLocalRef(ctx.da_solution, cell, fdof, &offset); CHKERRQ(ierr);
              slist = (F2I *) &offset[AS_OFFSET];
              if (slist->i[0] > 1) {
                  ierr = DMLabelSetValue(adaptlabel, cell, DM_ADAPT_REFINE); CHKERRQ(ierr);
              }    
          }
          ierr = VecRestoreArray(lsolution, &fdof); CHKERRQ(ierr);
          ierr = DMRestoreLocalVector(ctx.da_solution,&lsolution);CHKERRQ(ierr);

          ierr = DMAdaptLabel(ctx.da_solforest,adaptlabel,&postsolforest);CHKERRQ(ierr);
          success = postsolforest != NULL;
          if (success) {
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
              
              PetscFE output_fe;
              ierr = DMDestroy(&ctx.da_output);
              ierr = DMClone(ctx.da_solution,&ctx.da_output);
              ierr = PetscFECreateDefault(PETSC_COMM_WORLD,3,
                                          1+DP_SIZE,
                                          PETSC_FALSE,NULL,PETSC_DEFAULT,&output_fe);CHKERRQ(ierr);
              ierr = PetscObjectSetName((PetscObject) output_fe, "output");CHKERRQ(ierr);
              ierr = DMSetField(ctx.da_output,0, NULL, (PetscObject) output_fe);CHKERRQ(ierr);
              ierr = DMCreateDS(ctx.da_output);CHKERRQ(ierr);
              ierr = DMPlexCreateClosureIndex(ctx.da_output, NULL);
              ierr = PetscFEDestroy(&output_fe);CHKERRQ(ierr);
              
              free(ctx.cellgeom);
              ctx.cellgeom = malloc(ctx.ninteriorcells*sizeof(PetscReal));
              for (PetscInt cell = 0; cell < ctx.ninteriorcells; ++cell) {
                  PetscReal cvolume;
                  ierr = DMPlexComputeCellGeometryFVM(ctx.da_solution, cell, &cvolume, NULL, NULL);
                  ctx.cellgeom[cell] = cbrt(cvolume);
              }
              
              if (ts) ierr = TSDestroy(&ts);
              ierr = InitializeTS(ctx.da_solution, &ctx, &ts);CHKERRQ(ierr);
          }
          ierr = DMLabelDestroy(&adaptlabel);CHKERRQ(ierr);
      }    
      ierr = TSSetStepNumber(ts,nsteps);CHKERRQ(ierr);
      ierr = TSSetTime(ts,currenttime);CHKERRQ(ierr);
      ierr = TSSetTimeStep(ts,ctx.params.timestep);CHKERRQ(ierr);
      ierr = TSSetMaxSteps(ts,nsteps+ctx.params.amrinterval);CHKERRQ(ierr);
      ierr = TSSolve(ts,solution);CHKERRQ(ierr);
      ierr = TSGetSolveTime(ts,&currenttime);CHKERRQ(ierr);
      ierr = TSGetStepNumber(ts,&nsteps);CHKERRQ(ierr);
      ierr = TSGetTimeStep(ts,&ctx.params.timestep);CHKERRQ(ierr);
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
