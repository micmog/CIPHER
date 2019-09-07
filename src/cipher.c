
static char help[] = "Solves multi phase field equations \n";

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <petscdm.h>
#include <petscdmplex.h>
#include <petscdmforest.h>
#include <petscds.h>
#include <petscts.h>
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
    Vec               localX, matstate0;
    PetscScalar       *fdof;
    PetscScalar       *rhs, *mat, *superset, *fvmgeom;
    FIELD             *ucell, *nucell, *fcell, lcell;
    F2I               *slist, *nslist;
    STATE             *mcell;
    PetscInt          conesize, nsupp;
    PetscScalar       *edgel, *nedgel, celledge, ncelledge, cellvol, cfactor;
    PetscReal         interpolant[MAXAP], caplflux[MAXAP], caplsource[MAXAP];
    PetscReal         chemsource[MAXAP], rhs_unconstrained[MAXAP], cvgdot_phase[MAXAP];
    PetscBool         active[MAXAP];
    PetscReal         interfacepotential[MAXCP], chemicalpotential[MAXCP];
    PetscReal         mobilitycv[MAXCP], cavgdot[MAXCP];
    PetscReal         dcdm[MAXCP*MAXCP], dmdc[MAXCP*MAXCP];
    const PetscInt    *cone, *scells;
    uint16_t          setintersect[MAXIP], injectionL[MAXIP], injectionR[MAXIP];
    PetscReal         nactivephases, intval, rhsval, triplejunctionenergy;
    PetscInt          localcell, cell, face, supp; 
    PetscInt          g, gi, gj, gk, c, ci, cj, interfacekj, interfaceji, interfaceki;
    INTERFACE         *currentinterface;
    
    /* Gather FVM residuals */
    ierr = DMGetLocalVector(user->da_solution,&localX); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(user->da_solution,X,INSERT_VALUES,localX); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->da_solution,X,INSERT_VALUES,localX); CHKERRQ(ierr);
    ierr = VecGetArray(localX, &fdof);
    ierr = VecZeroEntries(F); CHKERRQ(ierr);
    ierr = VecGetArray(F, &rhs);
    ierr = DMGetGlobalVector(user->da_matstate,&matstate0); CHKERRQ(ierr);
    ierr = VecCopy(user->matstate,matstate0); CHKERRQ(ierr);
    ierr = VecGetArray(matstate0,&mat); CHKERRQ(ierr);
    ierr = VecGetArray(user->activephasesuperset,&superset); CHKERRQ(ierr);
    ierr = VecGetArray(user->fvmgeom,&fvmgeom); CHKERRQ(ierr);

    /* Loop over cells and compute local contribution to the RHS */
    for (localcell = 0; localcell < user->nlocalcells; ++localcell) {
        cell = user->localcells[localcell];
        /* get fields */
        
        slist = NULL; ucell = NULL; fcell = NULL; mcell = NULL;
        ierr = DMPlexPointLocalRead(user->da_phaseID,cell,superset,&slist);
        ierr = DMPlexPointLocalRead(user->da_solution,cell,fdof,&ucell);
        ierr = DMPlexPointGlobalRef(user->da_solution,cell,rhs,&fcell);
        ierr = DMPlexPointGlobalRef(user->da_matstate,cell,mat,&mcell);
        
        /* calculate phase interpolants */
        EvalInterpolant(interpolant,ucell->p,slist->i[0]);
                
        /* update material state */
        memset(interfacepotential,0,user->nc*sizeof(PetscScalar));
        memcpy(chemicalpotential,ucell->m,(user->nc-1)*sizeof(PetscScalar));
        
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
        user->rejectstage = Composition(mcell->c,chemicalpotential,slist->i,user);
    
        if (slist->i[0] > 1 && user->nc > 1) {
            /* get laplacian from neighbouring cells */
            ierr = DMPlexGetConeSize(user->da_solution,cell,&conesize);
            ierr = DMPlexGetCone    (user->da_solution,cell,&cone    );
            edgel = NULL;
            ierr = DMPlexPointLocalRead(user->da_fvmgeom,cell,fvmgeom,&edgel);
            celledge = edgel[0]; cellvol = celledge*celledge*celledge;
            memset(lcell.p,0,MAXAP*sizeof(PetscReal));
            memset(lcell.m,0,MAXCP*sizeof(PetscReal));
            for (face=0; face<conesize; face++) {
                ierr = DMPlexGetSupportSize(user->da_solution, cone[face], &nsupp);
                ierr = DMPlexGetSupport(user->da_solution, cone[face], &scells);
                for (supp=0; supp<nsupp; supp++) {
                    if (cell != scells[supp]) {
                        nslist = NULL; nucell = NULL; nedgel = NULL;
                        ierr = DMPlexPointLocalRead(user->da_phaseID,scells[supp],superset,&nslist);
                        ierr = DMPlexPointLocalRead(user->da_solution,scells[supp],fdof,&nucell);
                        ierr = DMPlexPointLocalRead(user->da_fvmgeom,scells[supp],fvmgeom,&nedgel);
                        ncelledge = nedgel[0];
                        cfactor  = ncelledge < celledge ? ncelledge : celledge;
                        cfactor *= 2.0/(ncelledge + celledge)/cellvol;
            
                        /* get common active phases */
                        SetIntersection(setintersect,injectionL,injectionR,slist->i,nslist->i);
            
                        /* compute fluxes */
                        for (g=0; g<setintersect[0]; g++) 
                            lcell.p[injectionL[g]] += cfactor*(nucell->p[injectionR[g]] - ucell->p[injectionL[g]]);
                        for (c=0; c<user->nc-1; c++) 
                            lcell.m[           c ] += cfactor*(nucell->m[           c ] - ucell->m[           c ]);
                    }
                }
            }
    
            /* phase capillary driving force */ 
            memset(caplflux  ,0,slist->i[0]*sizeof(PetscReal));
            memset(caplsource,0,slist->i[0]*sizeof(PetscReal));
            for (gk=0; gk<slist->i[0]; gk++) {
                for (gj=gk+1; gj<slist->i[0]; gj++) {
                    interfacekj = (PetscInt) user->interfacelist[slist->i[gk+1]*user->np + slist->i[gj+1]];
                    currentinterface = &user->interface[interfacekj];
                    caplflux  [gk] -= currentinterface->energy*lcell.p[gj];
                    caplflux  [gj] -= currentinterface->energy*lcell.p[gk];
                    caplsource[gk] -= currentinterface->energy*ucell->p[gj];
                    caplsource[gj] -= currentinterface->energy*ucell->p[gk];
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
                        caplsource[gk] -= triplejunctionenergy*ucell->p[gj]*ucell->p[gi];
                        caplsource[gj] -= triplejunctionenergy*ucell->p[gk]*ucell->p[gi];
                        caplsource[gi] -= triplejunctionenergy*ucell->p[gk]*ucell->p[gj];
                    }
                }
                caplflux  [gk] *= 8.0*user->params.interfacewidth/PETSC_PI/PETSC_PI;
                caplsource[gk] *= 8.0/user->params.interfacewidth;
            }   

            /* phase chemical driving force */ 
            Chemenergy(chemsource,mcell->c,ucell->m,slist->i,user);
            MatMulInterpolantDerivative(chemsource,ucell->p,slist->i[0]);

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
                active[gk] =    (ucell->p[gk] >       TOL && ucell->p         [gk] < 1.0 - TOL)
                             || (ucell->p[gk] <       TOL && rhs_unconstrained[gk] > 0.0      )
                             || (ucell->p[gk] > 1.0 - TOL && rhs_unconstrained[gk] < 0.0      );
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
                            fcell->p[gk] += rhsval; fcell->p[gj] -= rhsval;
                        }
                    }
                    fcell->p[gk] /= nactivephases;
                }        
            }
            
            /* calculate solute mobility matrix (volume-fixed frame of reference) */
            CompositionMobility(mobilitycv,mcell->c,interpolant,slist->i,user);
  
            /* calculate avg composition rate */
            for (c=0; c<user->nc-1; c++) {
                cavgdot[c] = mobilitycv[c]*lcell.m[c];
                for (g =0; g<slist->i[0];  g++) cvgdot_phase[g] = mcell->c[g*user->nc+c];
                MatMulInterpolantDerivative(cvgdot_phase,ucell->p,slist->i[0]);
                for (g =0; g<slist->i[0];  g++) cavgdot[c] -= fcell->p[g]*cvgdot_phase[g];
            }
  
            /* calculate composition tangent wrt chemical potential */
            CompositionTangent(dcdm ,mcell->c,interpolant,slist->i,user);
            Invertmatrix(dmdc,dcdm);
  
            /* chemical potential RHS */
            for (cj=0; cj<user->nc-1; cj++) {
                for (ci=0; ci<user->nc-1; ci++) {
                    fcell->m[cj] += dmdc[cj*(user->nc-1)+ci]*cavgdot[ci];
                }
            }
        } else if (slist->i[0] > 1) {
            /* get laplacian from neighbouring cells */
            ierr = DMPlexGetConeSize(user->da_solution,cell,&conesize);
            ierr = DMPlexGetCone    (user->da_solution,cell,&cone    );
            edgel = NULL;
            ierr = DMPlexPointLocalRead(user->da_fvmgeom,cell,fvmgeom,&edgel);
            celledge = edgel[0]; cellvol = celledge*celledge*celledge;
            memset(lcell.p,0,MAXAP*sizeof(PetscReal));
            for (face=0; face<conesize; face++) {
                ierr = DMPlexGetSupportSize(user->da_solution, cone[face], &nsupp);
                ierr = DMPlexGetSupport(user->da_solution, cone[face], &scells);
                for (supp=0; supp<nsupp; supp++) {
                    if (cell != scells[supp]) {
                        nslist = NULL; nucell = NULL; nedgel = NULL;
                        ierr = DMPlexPointLocalRead(user->da_phaseID,scells[supp],superset,&nslist);
                        ierr = DMPlexPointLocalRead(user->da_solution,scells[supp],fdof,&nucell);
                        ierr = DMPlexPointLocalRead(user->da_fvmgeom,scells[supp],fvmgeom,&nedgel);
                        ncelledge = nedgel[0];
                        cfactor  = ncelledge < celledge ? ncelledge : celledge;
                        cfactor *= 2.0/(ncelledge + celledge)/cellvol;
            
                        /* get common active phases */
                        SetIntersection(setintersect,injectionL,injectionR,slist->i,nslist->i);
            
                        /* compute fluxes */
                        for (g=0; g<setintersect[0]; g++) 
                            lcell.p[injectionL[g]] += cfactor*(nucell->p[injectionR[g]] - ucell->p[injectionL[g]]);
                    }
                }
            }
    
            /* phase capillary driving force */ 
            memset(caplflux  ,0,slist->i[0]*sizeof(PetscReal));
            memset(caplsource,0,slist->i[0]*sizeof(PetscReal));
            for (gk=0; gk<slist->i[0]; gk++) {
                for (gj=gk+1; gj<slist->i[0]; gj++) {
                    interfacekj = (PetscInt) user->interfacelist[slist->i[gk+1]*user->np + slist->i[gj+1]];
                    currentinterface = &user->interface[interfacekj];
                    caplflux  [gk] -= currentinterface->energy*lcell.p[gj];
                    caplflux  [gj] -= currentinterface->energy*lcell.p[gk];
                    caplsource[gk] -= currentinterface->energy*ucell->p[gj];
                    caplsource[gj] -= currentinterface->energy*ucell->p[gk];
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
                        caplsource[gk] -= triplejunctionenergy*ucell->p[gj]*ucell->p[gi];
                        caplsource[gj] -= triplejunctionenergy*ucell->p[gk]*ucell->p[gi];
                        caplsource[gi] -= triplejunctionenergy*ucell->p[gk]*ucell->p[gj];
                    }
                }
                caplflux  [gk] *= 8.0*user->params.interfacewidth/PETSC_PI/PETSC_PI;
                caplsource[gk] *= 8.0/user->params.interfacewidth;
            }   

            /* phase chemical driving force */ 
            Chemenergy(chemsource,mcell->c,ucell->m,slist->i,user);
            MatMulInterpolantDerivative(chemsource,ucell->p,slist->i[0]);

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
                active[gk] =    (ucell->p[gk] >       TOL && ucell->p         [gk] < 1.0 - TOL)
                             || (ucell->p[gk] <       TOL && rhs_unconstrained[gk] > 0.0      )
                             || (ucell->p[gk] > 1.0 - TOL && rhs_unconstrained[gk] < 0.0      );
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
                            fcell->p[gk] += rhsval; fcell->p[gj] -= rhsval;
                        }
                    }
                    fcell->p[gk] /= nactivephases;
                }        
            }
        } else if (user->nc > 1) {
            /* get laplacian from neighbouring cells */
            ierr = DMPlexGetConeSize(user->da_solution,cell,&conesize);
            ierr = DMPlexGetCone    (user->da_solution,cell,&cone    );
            edgel = NULL;
            ierr = DMPlexPointLocalRead(user->da_fvmgeom,cell,fvmgeom,&edgel);
            celledge = edgel[0]; cellvol = celledge*celledge*celledge;
            memset(lcell.m,0,MAXCP*sizeof(PetscReal));
            for (face=0; face<conesize; face++) {
                ierr = DMPlexGetSupportSize(user->da_solution, cone[face], &nsupp);
                ierr = DMPlexGetSupport(user->da_solution, cone[face], &scells);
                for (supp=0; supp<nsupp; supp++) {
                    if (cell != scells[supp]) {
                        /* compute fluxes */
                        nucell = NULL; nedgel = NULL;
                        ierr = DMPlexPointLocalRead(user->da_solution,scells[supp],fdof,&nucell);
                        ierr = DMPlexPointLocalRead(user->da_fvmgeom,scells[supp],fvmgeom,&nedgel);
                        ncelledge = nedgel[0];
                        cfactor  = ncelledge < celledge ? ncelledge : celledge;
                        cfactor *= 2.0/(ncelledge + celledge)/cellvol;
            
                        /* compute fluxes */
                        for (c=0; c<user->nc-1; c++) 
                            lcell.m[c] += cfactor*(nucell->m[c] - ucell->m[c]);
                    }
                }
            }

            /* calculate solute mobility matrix (volume-fixed frame of reference) */
            CompositionMobility(mobilitycv,mcell->c,interpolant,slist->i,user);
  
            /* calculate avg composition rate */
            for (c=0; c<user->nc-1; c++) {
                cavgdot[c] = mobilitycv[c]*lcell.m[c];
                for (g =0; g<slist->i[0];  g++) cvgdot_phase[g] = mcell->c[g*user->nc+c];
                MatMulInterpolantDerivative(cvgdot_phase,ucell->p,slist->i[0]);
                for (g =0; g<slist->i[0];  g++) cavgdot[c] -= fcell->p[g]*cvgdot_phase[g];
            }
  
            /* calculate composition tangent wrt chemical potential */
            CompositionTangent(dcdm ,mcell->c,interpolant,slist->i,user);
            Invertmatrix(dmdc,dcdm);
  
            /* chemical potential RHS */
            for (cj=0; cj<user->nc-1; cj++) {
                for (ci=0; ci<user->nc-1; ci++) {
                    fcell->m[cj] += dmdc[cj*(user->nc-1)+ci]*cavgdot[ci];
                }
            }
        }
    }
    
    /* Restore FVM residuals */
    ierr = VecRestoreArray(localX, &fdof);
    ierr = DMRestoreLocalVector(user->da_solution,&localX); CHKERRQ(ierr);
    ierr = VecRestoreArray(F, &rhs);
    ierr = VecRestoreArray(matstate0,&mat); CHKERRQ(ierr);
    ierr = DMRestoreGlobalVector(user->da_matstate,&matstate0); CHKERRQ(ierr);
    ierr = VecRestoreArray(user->activephasesuperset,&superset); CHKERRQ(ierr);
    ierr = VecRestoreArray(user->fvmgeom,&fvmgeom); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}    

/*
 PostStep - Update active phases and write output
 */
PetscErrorCode PostStep(TS ts)
{
    PetscErrorCode ierr;
    AppCtx         *user;
    Vec            solution, gactivephaseset;
    PetscScalar    *fdof, *mat, *activeset, *superset, *gsuperset;
    FIELD          *ucell;
    F2I            *slist, *galist, *gslist, *nalist;
    STATE          *mcell;
    PetscInt       conesize, nsupp;
    const PetscInt *cone, *scells;
    PetscReal      intval;
    PetscReal      phiUP[MAXAP], interpolant[MAXAP], phiSS[MAXAP], compSS[MAXAP*MAXCP];
    PetscReal      interfacepotential[MAXCP], chemicalpotential[MAXCP];
    PetscInt       localcell, cell, face, supp, g, gj, gk, c, interfacekj;
    INTERFACE      *currentinterface;
    MATERIAL       *currentmaterial;
    uint16_t       setunion[MAXIP], setintersection[MAXIP], injectionL[MAXIP], injectionR[MAXIP];
    
    PetscFunctionBeginUser;    
    ierr = TSGetSolution(ts, &solution);CHKERRQ(ierr);    
    ierr = TSGetApplicationContext(ts,&user);CHKERRQ(ierr);
      
    /* Determine active phase set, update composition */
    ierr = DMGetGlobalVector(user->da_phaseID,&gactivephaseset); CHKERRQ(ierr);
    ierr = VecGetArray(gactivephaseset,&activeset); CHKERRQ(ierr);
    ierr = VecGetArray(user->activephasesuperset,&superset); CHKERRQ(ierr);
    ierr = VecGetArray(solution,&fdof);
    ierr = VecGetArray(user->matstate,&mat); CHKERRQ(ierr);
    for (localcell = 0; localcell < user->nlocalcells; ++localcell) {
        cell = user->localcells[localcell];
        /* get cell state */
        ucell = NULL; mcell = NULL; slist = NULL; galist = NULL;
        ierr = DMPlexPointGlobalRef(user->da_solution,cell,fdof,&ucell);
        ierr = DMPlexPointGlobalRef(user->da_matstate,cell,mat,&mcell);
        ierr = DMPlexPointGlobalRef(user->da_phaseID,cell,activeset,&galist);
        ierr = DMPlexPointLocalRead(user->da_phaseID,cell,superset,&slist);

        if (slist->i[0] > 1) {
            /* project phase fields back to Gibbs simplex */
            memcpy(phiUP,ucell->p,slist->i[0]*sizeof(PetscReal));
            SimplexProjection(ucell->p,phiUP,slist->i[0]);
            EvalInterpolant(interpolant,ucell->p,slist->i[0]);
    
            /* update material state */
            memset(interfacepotential,0,user->nc*sizeof(PetscScalar));
            memcpy(chemicalpotential,ucell->m,(user->nc-1)*sizeof(PetscScalar));
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
            galist->i[0] = 0;
            for (g=0; g<slist->i[0];  g++) {
                if (interpolant[g] > TOL) {
                    galist->i[++(galist->i[0])] = slist->i[g+1];
                }
            }
        } else {
            /* project phase fields back to Gibbs simplex */
            ucell->p[0] = 1.0;
    
            /* update material state */
            Composition(mcell->c,ucell->m,slist->i,user);
        
            /* update active set */
            galist->i[0] = 1; galist->i[1] = slist->i[1];
        }
    }
    ierr = VecRestoreArray(solution,&fdof);
    ierr = VecRestoreArray(user->matstate,&mat); CHKERRQ(ierr);
    ierr = VecRestoreArray(user->activephasesuperset,&superset); CHKERRQ(ierr);
    ierr = VecRestoreArray(gactivephaseset,&activeset); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(user->da_phaseID,gactivephaseset,INSERT_VALUES,user->activephaseset); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->da_phaseID,gactivephaseset,INSERT_VALUES,user->activephaseset); CHKERRQ(ierr);

    /* Determine active phase super set, reorder dofs in solution, initialize new phase frac & comp */
    ierr = VecGetArray(user->activephaseset,&activeset); CHKERRQ(ierr);
    ierr = VecGetArray(gactivephaseset,&gsuperset); CHKERRQ(ierr);
    ierr = VecGetArray(user->activephasesuperset,&superset); CHKERRQ(ierr);
    ierr = VecGetArray(solution,&fdof); CHKERRQ(ierr);
    ierr = VecGetArray(user->matstate,&mat); CHKERRQ(ierr);
    for (localcell = 0; localcell < user->nlocalcells; ++localcell) {
        cell = user->localcells[localcell];
        /* get cell state */
        ucell = NULL; mcell = NULL; slist = NULL; gslist = NULL;
        ierr = DMPlexPointGlobalRef(user->da_solution,cell,fdof,&ucell);
        ierr = DMPlexPointGlobalRef(user->da_matstate,cell,mat,&mcell);
        ierr = DMPlexPointGlobalRef(user->da_phaseID,cell,gsuperset,&gslist);
        ierr = DMPlexPointLocalRead(user->da_phaseID,cell,superset,&slist);

        /* add union of neighbouring cells */
        ierr = DMPlexGetConeSize(user->da_solution,cell,&conesize);
        ierr = DMPlexGetCone    (user->da_solution,cell,&cone    );
        for (face=0; face<conesize; face++) {
            ierr = DMPlexGetSupportSize(user->da_solution, cone[face], &nsupp);
            ierr = DMPlexGetSupport(user->da_solution, cone[face], &scells);
            for (supp=0; supp<nsupp; supp++) {
                if (cell != scells[supp]) {
                    nalist = NULL;
                    ierr = DMPlexPointLocalRead(user->da_phaseID,scells[supp],activeset,&nalist);
                    SetUnion(setunion,injectionL,injectionR,gslist->i,nalist->i);
                    memcpy(gslist->i,setunion,MAXIP*sizeof(uint16_t));
                }
            }
        }

        /* initialize new phase states */
        for (g=0; g<gslist->i[0];  g++) {
            phiSS[g] = 0.0;
            currentmaterial = &user->material[user->phasematerialmapping[gslist->i[g+1]]];
            memcpy(&compSS[g*user->nc],currentmaterial->c0,user->nc*sizeof(PetscReal));
        }
        
        /* reorder dofs to new active phase superset */
        SetIntersection(setintersection,injectionL,injectionR,slist->i,gslist->i);
        for (g=0; g<setintersection[0];  g++) {
            phiSS[injectionR[g]] = ucell->p[injectionL[g]];
            memcpy(&compSS[injectionR[g]*user->nc],&mcell->c[injectionL[g]*user->nc],user->nc*sizeof(PetscReal));
        }
        memcpy(ucell->p,phiSS,gslist->i[0]*sizeof(PetscReal));
        memcpy(mcell->c,compSS,gslist->i[0]*user->nc*sizeof(PetscReal));
    }    
    ierr = VecRestoreArray(solution,&fdof); CHKERRQ(ierr);
    ierr = VecRestoreArray(user->matstate,&mat); CHKERRQ(ierr);
    ierr = VecRestoreArray(gactivephaseset,&gsuperset); CHKERRQ(ierr);
    ierr = VecRestoreArray(user->activephaseset,&activeset); CHKERRQ(ierr);
    ierr = VecRestoreArray(user->activephasesuperset,&superset); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(user->da_phaseID,gactivephaseset,INSERT_VALUES,user->activephasesuperset); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->da_phaseID,gactivephaseset,INSERT_VALUES,user->activephasesuperset); CHKERRQ(ierr);
    ierr = DMRestoreGlobalVector(user->da_phaseID,&gactivephaseset); CHKERRQ(ierr);
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
       ierr = VecGetArray(user->activephasesuperset,&superset); CHKERRQ(ierr);
       ierr = VecGetArray(Xout,&xout); CHKERRQ(ierr);
       for (localcell = 0; localcell < user->nlocalcells; ++localcell) {
           cell = user->localcells[localcell];
           /* get cell state */
           ucell = NULL; mcell = NULL; slist = NULL; ocell = NULL;
           ierr = DMPlexPointGlobalRef(user->da_solution,cell,fdof,&ucell);
           ierr = DMPlexPointGlobalRef(user->da_matstate,cell,mat,&mcell);
           ierr = DMPlexPointLocalRead(user->da_phaseID,cell,superset,&slist);
           ierr = DMPlexPointGlobalRef(user->da_output,cell,xout,&ocell);
       
           PetscReal max = -LARGE, interpolant[slist->i[0]];
           EvalInterpolant(interpolant,ucell->p,slist->i[0]);
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

static PetscErrorCode AMRfunction(TS ts, Vec sol, VecTagger refineTag, VecTagger coarsenTag, AppCtx *user)
{
    PetscErrorCode ierr;
    DM                adapted_solution=NULL, adapted_matstate=NULL, adapted_output=NULL;
    DM                dforest_solution=NULL, dforest_matstate=NULL, dforest_output=NULL;
    Vec               localX, laplacian, errVec;
    const PetscScalar *fdof;
    PetscScalar       *lap;
    PetscInt          face, fstart, fend, cell, cstart, cend, nRefine, nCoarsen;
    DMLabel           ghostlabel = NULL;
    IS                refineIS, coarsenIS;
    PetscReal         time;
    
    PetscFunctionBegin;
    ierr = TSGetTime(ts,&time);CHKERRQ(ierr);

    /* Gather FVM residuals */
    ierr = DMGetLocalVector(user->da_solution,&localX); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(user->da_solution,sol,INSERT_VALUES,localX); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->da_solution,sol,INSERT_VALUES,localX); CHKERRQ(ierr);
    ierr = VecGetArrayRead(localX, &fdof);
    ierr = DMGetGlobalVector(user->da_solution,&laplacian); CHKERRQ(ierr);
    ierr = VecZeroEntries(laplacian); CHKERRQ(ierr);
    ierr = VecGetArray(laplacian, &lap);
    
    /* Loop over faces and add error indicator contribution to each cell */
    ierr = DMGetLabel(user->da_solution, "ghost", &ghostlabel);
    ierr = DMPlexGetHeightStratum(user->da_solution, 0, &cstart, &cend); CHKERRQ(ierr);
    ierr = DMPlexGetHeightStratum(user->da_solution, 1, &fstart, &fend); CHKERRQ(ierr);
    for (face = fstart; face < fend; ++face) {
        const PetscInt *scells;
        FIELD *fL, *fR, flux; 
        const FIELD *uL, *uR;
        PetscInt ghost, nsupp, nchild;

        ierr = DMLabelGetValue(ghostlabel, face, &ghost);
        ierr = DMPlexGetSupportSize(user->da_solution, face, &nsupp);
        ierr = DMPlexGetTreeChildren(user->da_solution, face, &nchild, NULL);
        
        /* Skip ghost, boundary or parent faces */
        if (ghost < 0 && nsupp < 3 && nchild == 0) {
            ierr = DMPlexGetSupport(user->da_solution, face, &scells);
            ierr = DMLabelGetValue(ghostlabel,scells[0],&ghost);
            
            /* get field */
            ierr = DMPlexPointLocalRead(user->da_solution, scells[0], fdof, &uL);
            ierr = DMPlexPointLocalRead(user->da_solution, scells[1], fdof, &uR);
            
            /* compute cel activity */
            for (PetscInt g=0; g<user->np  ; g++) flux.p[g] = (uL->p[g] - uR->p[g]);
            ierr = DMLabelGetValue(ghostlabel,scells[0],&ghost);
            if (ghost <= 0) {
                ierr = DMPlexPointGlobalRef(user->da_solution, scells[0], lap, &fL);
                for (PetscInt g=0; g<user->np  ; g++) fL->p[g] -= flux.p[g];
            }
            ierr = DMLabelGetValue(ghostlabel,scells[1],&ghost);
            if (ghost <= 0) {
                ierr = DMPlexPointGlobalRef(user->da_solution, scells[1], lap, &fR);
                for (PetscInt g=0; g<user->np  ; g++) fR->p[g] -= flux.p[g];
            }
        }    
    }
    
    /* set up error indicator vector */
    DMLabel adaptlabel = NULL;
    PetscScalar *errArray;
    ierr = DMLabelCreate(PETSC_COMM_SELF,"adapt",&adaptlabel);CHKERRQ(ierr);
    ierr = VecCreateMPI(PetscObjectComm((PetscObject)user->da_solution),cend-cstart,PETSC_DETERMINE,&errVec);CHKERRQ(ierr);
    ierr = VecSetUp(errVec);CHKERRQ(ierr);
    ierr = VecGetArray(errVec,&errArray);CHKERRQ(ierr);
    for (cell=cstart; cell<cend; cell++) {
        const FIELD *lcell = NULL;
        ierr = DMPlexPointGlobalRead(user->da_solution,cell,lap,&lcell);CHKERRQ(ierr);
        errArray[cell-cstart] = 0.0;
        for (PetscInt g=0; g<user->np  ; g++) errArray[cell-cstart] += fabs(lcell->p[g]);
    }
    ierr = VecRestoreArray(errVec,&errArray);CHKERRQ(ierr);

    /* Restore FVM residuals */
    ierr = VecRestoreArrayRead(localX, &fdof);
    ierr = DMRestoreLocalVector(user->da_solution,&localX); CHKERRQ(ierr);
    ierr = VecRestoreArray(laplacian, &lap);
    ierr = DMRestoreGlobalVector(user->da_solution,&laplacian); CHKERRQ(ierr);

    /* Destroy old DM */
    ierr = VecTaggerComputeIS(refineTag,errVec,&refineIS);CHKERRQ(ierr);
    ierr = VecTaggerComputeIS(coarsenTag,errVec,&coarsenIS);CHKERRQ(ierr);
    ierr = ISGetSize(refineIS ,&nRefine );CHKERRQ(ierr);
    ierr = ISGetSize(coarsenIS,&nCoarsen);CHKERRQ(ierr);
    if (nRefine ) {ierr = DMLabelSetStratumIS(adaptlabel,DM_ADAPT_REFINE ,refineIS );CHKERRQ(ierr);}
    if (nCoarsen) {ierr = DMLabelSetStratumIS(adaptlabel,DM_ADAPT_COARSEN,coarsenIS);CHKERRQ(ierr);}
    ierr = ISDestroy(&coarsenIS);CHKERRQ(ierr);
    ierr = ISDestroy(&refineIS );CHKERRQ(ierr);
    ierr = VecDestroy(&errVec);CHKERRQ(ierr);

    if (nRefine || nCoarsen) { /* at least one cell is over the refinement threshold */
        ierr = DMConvert(user->da_solution,DMP8EST,&dforest_solution);CHKERRQ(ierr);
        ierr = DMAdaptLabel(dforest_solution,adaptlabel,&adapted_solution);CHKERRQ(ierr);
        ierr = DMConvert(user->da_matstate,DMP8EST,&dforest_matstate);CHKERRQ(ierr);
        ierr = DMAdaptLabel(dforest_matstate,adaptlabel,&adapted_matstate);CHKERRQ(ierr);
        ierr = DMConvert(user->da_output  ,DMP8EST,&dforest_output  );CHKERRQ(ierr);
        ierr = DMAdaptLabel(dforest_output  ,adaptlabel,&adapted_output  );CHKERRQ(ierr);
    }
    ierr = DMLabelDestroy(&adaptlabel);CHKERRQ(ierr);
    if (adapted_solution) {
        Vec solnew = NULL;
        ierr = DMForestTransferVec(dforest_solution, sol, adapted_solution, solnew, PETSC_TRUE, time);CHKERRQ(ierr);
        ierr = DMDestroy(&dforest_solution);CHKERRQ(ierr);
        ierr = DMDestroy(&user->da_solution);CHKERRQ(ierr);
        ierr = DMConvert(adapted_solution,DMPLEX,&user->da_solution);CHKERRQ(ierr);
        ierr = PetscObjectReference((PetscObject) user->da_solution);CHKERRQ(ierr);
        ierr = VecDestroy(&sol);CHKERRQ(ierr);
        sol = solnew;
        ierr = PetscObjectSetName((PetscObject) sol, "solution");CHKERRQ(ierr);
        ierr = TSDestroy(&ts);CHKERRQ(ierr);
        ierr = InitializeTS(user->da_solution, user, &ts);CHKERRQ(ierr);
    }
    if (adapted_matstate) {
        Vec matnew = NULL;
        ierr = DMForestTransferVec(dforest_matstate, user->matstate, adapted_matstate, matnew, PETSC_TRUE, time);CHKERRQ(ierr);
        ierr = DMDestroy(&dforest_matstate);CHKERRQ(ierr);
        ierr = DMDestroy(&user->da_matstate);CHKERRQ(ierr);
        ierr = DMConvert(adapted_matstate,DMPLEX,&user->da_matstate);CHKERRQ(ierr);
        ierr = VecDestroy(&user->matstate);CHKERRQ(ierr);
        user->matstate = matnew;
        ierr = PetscObjectSetName((PetscObject) user->matstate, "material state");CHKERRQ(ierr);
    }
    if (adapted_output) {
        ierr = DMDestroy(&dforest_output);CHKERRQ(ierr);
        ierr = DMDestroy(&user->da_output);CHKERRQ(ierr);
        ierr = DMConvert(adapted_output,DMPLEX,&user->da_output);CHKERRQ(ierr);
    }
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
  PetscLimiter   nonelimiter = NULL;
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
  /* Step up AMR */
  {
    VecTaggerBox refineBox, coarsenBox;

    refineBox.min  = refineBox.max  = PETSC_MAX_REAL;
    coarsenBox.min = coarsenBox.max = PETSC_MIN_REAL;

    ierr = VecTaggerCreate(PETSC_COMM_WORLD,&refineTag);CHKERRQ(ierr);
    ierr = PetscObjectSetOptionsPrefix((PetscObject)refineTag,"refine_");CHKERRQ(ierr);
    ierr = VecTaggerSetType(refineTag,VECTAGGERABSOLUTE);CHKERRQ(ierr);
    ierr = VecTaggerAbsoluteSetBox(refineTag,&refineBox);CHKERRQ(ierr);
    ierr = VecTaggerSetFromOptions(refineTag);CHKERRQ(ierr);
    ierr = VecTaggerSetUp(refineTag);CHKERRQ(ierr);

    ierr = VecTaggerCreate(PETSC_COMM_WORLD,&coarsenTag);CHKERRQ(ierr);
    ierr = PetscObjectSetOptionsPrefix((PetscObject)coarsenTag,"coarsen_");CHKERRQ(ierr);
    ierr = VecTaggerSetType(coarsenTag,VECTAGGERABSOLUTE);CHKERRQ(ierr);
    ierr = VecTaggerAbsoluteSetBox(coarsenTag,&coarsenBox);CHKERRQ(ierr);
    ierr = VecTaggerSetFromOptions(coarsenTag);CHKERRQ(ierr);
    ierr = VecTaggerSetUp(coarsenTag);CHKERRQ(ierr);
  }

  /* Create and distribute mesh */
  DMBoundaryType boundarytype[3] = {DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC};
  PetscReal lower[3] = {0.0, 0.0, 0.0}, upper[3] = {(PetscReal) ctx.resolution[0], (PetscReal) ctx.resolution[1], (PetscReal) ctx.resolution[2]};
  ierr = DMPlexCreateBoxMesh(PETSC_COMM_WORLD, 3, PETSC_FALSE, ctx.resolution, lower, upper, boundarytype, PETSC_TRUE, &ctx.da_solution);
  {
    PetscInt pstart, pend;
    ierr = DMCreateLabel(ctx.da_solution, "phase");CHKERRQ(ierr);
    ierr = DMPlexGetHeightStratum(ctx.da_solution, 0, &pstart, &pend); CHKERRQ(ierr);
    for (PetscInt point = pstart; point < pend; point++) {ierr = DMSetLabelValue(ctx.da_solution, "phase", point, ctx.phasevoxelmapping[point]);CHKERRQ(ierr);}
  }
  {
    PetscPartitioner part;
    DM da_distributed = NULL;
    /* Distribute mesh over processes */
    ierr = DMSetBasicAdjacency(ctx.da_solution, PETSC_TRUE, PETSC_FALSE);CHKERRQ(ierr);
    ierr = DMPlexGetPartitioner(ctx.da_solution, &part);CHKERRQ(ierr);
    ierr = PetscPartitionerSetFromOptions(part);CHKERRQ(ierr);
    ierr = DMPlexDistribute(ctx.da_solution, 1, NULL, &da_distributed);CHKERRQ(ierr);
    if (da_distributed) {
      ierr = DMDestroy(&ctx.da_solution);CHKERRQ(ierr);
      ctx.da_solution  = da_distributed;
    }  
  }
  ierr = DMSetFromOptions(ctx.da_solution); CHKERRQ(ierr);
  ierr = DMClone(ctx.da_solution, &ctx.da_matstate); CHKERRQ(ierr);
  ierr = DMClone(ctx.da_solution, &ctx.da_output); CHKERRQ(ierr);
  ierr = DMLocalizeCoordinates(ctx.da_output);
  {
    DM da_ghosted = NULL;

    ierr = DMPlexConstructGhostCells(ctx.da_solution, NULL, NULL, &da_ghosted);CHKERRQ(ierr);
    if (da_ghosted) {
      ierr = DMDestroy(&ctx.da_solution);CHKERRQ(ierr);
      ctx.da_solution  = da_ghosted;
    }  
  }
  ierr = DMClone(ctx.da_solution, &ctx.da_fvmgeom); CHKERRQ(ierr);
  ierr = DMClone(ctx.da_solution, &ctx.da_phaseID); CHKERRQ(ierr);
  
  /* Create finite volume discretisation for solution */
  PetscFV solution_fvm;
  ierr = PetscFVCreate(PETSC_COMM_WORLD, &solution_fvm);CHKERRQ(ierr);
  ierr = PetscFVSetSpatialDimension(solution_fvm,3);CHKERRQ(ierr);
  ierr = PetscFVSetNumComponents(solution_fvm, MAXAP+MAXCP);CHKERRQ(ierr);
  ierr = PetscFVSetType(solution_fvm,PETSCFVLEASTSQUARES);CHKERRQ(ierr);
  ierr = PetscLimiterCreate(PETSC_COMM_WORLD, &nonelimiter);CHKERRQ(ierr);
  ierr = PetscLimiterSetType(nonelimiter, PETSCLIMITERNONE);CHKERRQ(ierr);
  ierr = PetscFVSetLimiter(solution_fvm, nonelimiter);CHKERRQ(ierr);
  ierr = PetscFVSetFromOptions(solution_fvm);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) solution_fvm,"solution");CHKERRQ(ierr);
  {
    for (PetscInt c=0; c < MAXAP; c++) {
        char compname[256]  = "unknown";
        ierr = PetscSNPrintf(compname,sizeof(compname),"order parameter %d",c);CHKERRQ(ierr);
        ierr = PetscFVSetComponentName(solution_fvm,c,compname);CHKERRQ(ierr);
    }
    for (PetscInt c=MAXAP; c < MAXAP+MAXCP; c++) {
        char compname[256]  = "unknown";
        if (c < MAXAP+ctx.nc)
        ierr = PetscSNPrintf(compname,sizeof(compname),"chemical potential %s",ctx.componentname[c-MAXAP]);CHKERRQ(ierr);
        ierr = PetscFVSetComponentName(solution_fvm,c,compname);CHKERRQ(ierr);
    }
  }
  ierr = DMSetField(ctx.da_solution, 0, NULL, (PetscObject) solution_fvm);CHKERRQ(ierr);
  ierr = DMCreateDS(ctx.da_solution);CHKERRQ(ierr);
  {
    PetscInt cstart, cend;
    DMLabel ghostlabel;

    ierr = DMGetLabel(ctx.da_solution, "ghost", &ghostlabel);
    ierr = DMPlexGetHeightStratum(ctx.da_solution, 0, &cstart, &cend); CHKERRQ(ierr);
    ctx.localcells = malloc(cend*sizeof(PetscInt)); ctx.nlocalcells = 0;
    for (PetscInt cell = cstart; cell < cend; ++cell) {
        PetscInt ghost;
        ierr = DMLabelGetValue(ghostlabel, cell, &ghost);
        if (ghost <= 0) ctx.localcells[(ctx.nlocalcells)++] = cell;
    }    
  }
  ierr = DMPlexCreateClosureIndex(ctx.da_solution, NULL);

  /* Precalculated finite volume geometry information */
  {
    DM da_coordinates;
    Vec coordinates;
    PetscSection geomsection, coordsection;
    PetscInt pstart, pend;
    
    ierr = PetscSectionCreate(PETSC_COMM_WORLD, &geomsection);
    ierr = DMPlexGetChart(ctx.da_fvmgeom, &pstart, &pend);
    ierr = PetscSectionSetChart(geomsection, pstart, pend);
    ierr = DMPlexGetHeightStratum(ctx.da_fvmgeom, 0, &pstart, &pend);
    for (PetscInt point=pstart; point<pend; ++point) ierr = PetscSectionSetDof(geomsection, point, 1);
    ierr = DMPlexGetHeightStratum(ctx.da_fvmgeom, 1, &pstart, &pend);
    for (PetscInt point=pstart; point<pend; ++point) ierr = PetscSectionSetDof(geomsection, point, 1);
    ierr = PetscSectionSetUp(geomsection);
    ierr = DMSetSection(ctx.da_fvmgeom, geomsection);
    PetscSectionDestroy(&geomsection);
    ierr = DMGetSection(ctx.da_fvmgeom, &geomsection);
    ierr = DMCreateLocalVector(ctx.da_fvmgeom,&ctx.fvmgeom);
    
    ierr = DMGetCoordinatesLocal(ctx.da_fvmgeom, &coordinates);
    ierr = DMGetCoordinateDM(ctx.da_fvmgeom, &da_coordinates);
    ierr = DMGetSection(da_coordinates, &coordsection);
    ierr = DMPlexGetHeightStratum(ctx.da_fvmgeom, 0, &pstart, &pend);
    for (PetscInt point=pstart; point<pend; ++point) {
        PetscInt coordsize;
        PetscReal *coords;
        PetscScalar edges;

        coords = NULL;
        DMPlexVecGetClosure(da_coordinates, coordsection, coordinates, point, &coordsize, &coords);
        edges = FastPow((  fabs(coords[0] - coords[3]) < (PetscReal) ctx.resolution[0] - fabs(coords[0] - coords[3]) 
                         ? fabs(coords[0] - coords[3]) : (PetscReal) ctx.resolution[0] - fabs(coords[0] - coords[3])),2)
              + FastPow((  fabs(coords[1] - coords[4]) < (PetscReal) ctx.resolution[1] - fabs(coords[1] - coords[4]) 
                         ? fabs(coords[1] - coords[4]) : (PetscReal) ctx.resolution[1] - fabs(coords[1] - coords[4])),2)           
              + FastPow((  fabs(coords[2] - coords[5]) < (PetscReal) ctx.resolution[2] - fabs(coords[2] - coords[5]) 
                         ? fabs(coords[2] - coords[5]) : (PetscReal) ctx.resolution[2] - fabs(coords[2] - coords[5])),2);        
        edges = sqrt(edges);
        DMPlexVecSetClosure(ctx.da_fvmgeom, geomsection , ctx.fvmgeom, point, &edges, INSERT_VALUES);
        DMPlexVecRestoreClosure(da_coordinates, coordsection, coordinates, point, &coordsize, &coords);
    }
    ierr = DMPlexGetHeightStratum(ctx.da_fvmgeom, 1, &pstart, &pend);
    for (PetscInt point=pstart; point<pend; ++point) {
        PetscInt coordsize;
        PetscReal *coords;
        PetscScalar edges;

        coords = NULL;
        DMPlexVecGetClosure(da_coordinates, coordsection, coordinates, point, &coordsize, &coords);
        edges = FastPow((  fabs(coords[0] - coords[3]) < (PetscReal) ctx.resolution[0] - fabs(coords[0] - coords[3]) 
                         ? fabs(coords[0] - coords[3]) : (PetscReal) ctx.resolution[0] - fabs(coords[0] - coords[3])),2)
              + FastPow((  fabs(coords[1] - coords[4]) < (PetscReal) ctx.resolution[1] - fabs(coords[1] - coords[4]) 
                         ? fabs(coords[1] - coords[4]) : (PetscReal) ctx.resolution[1] - fabs(coords[1] - coords[4])),2)           
              + FastPow((  fabs(coords[2] - coords[5]) < (PetscReal) ctx.resolution[2] - fabs(coords[2] - coords[5]) 
                         ? fabs(coords[2] - coords[5]) : (PetscReal) ctx.resolution[2] - fabs(coords[2] - coords[5])),2);        
        edges = sqrt(edges);
        DMPlexVecSetClosure(ctx.da_fvmgeom, geomsection , ctx.fvmgeom, point, &edges, INSERT_VALUES);
        DMPlexVecRestoreClosure(da_coordinates, coordsection, coordinates, point, &coordsize, &coords);
    }
  }
  ierr = DMPlexCreateClosureIndex(ctx.da_fvmgeom, NULL);

  /* Create finite volume discretisation for material state */
  PetscFV matstate_fvm;
  ierr = PetscFVCreate(PETSC_COMM_WORLD, &matstate_fvm);CHKERRQ(ierr);
  ierr = PetscFVSetSpatialDimension(matstate_fvm, 3);CHKERRQ(ierr);
  ierr = PetscFVSetNumComponents(matstate_fvm, MAXAP*MAXCP);CHKERRQ(ierr);
  ierr = PetscFVSetLimiter(matstate_fvm, nonelimiter);CHKERRQ(ierr);
  ierr = PetscFVSetFromOptions(matstate_fvm);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) matstate_fvm,"material state");CHKERRQ(ierr);
  ierr = DMSetField(ctx.da_matstate, 0, NULL, (PetscObject) matstate_fvm);CHKERRQ(ierr);
  ierr = DMCreateDS(ctx.da_matstate);CHKERRQ(ierr);
  ierr = DMPlexCreateClosureIndex(ctx.da_matstate, NULL);

  /* Create finite volume discretisation for material state */
  PetscFV phaseID_fvm;
  ierr = PetscFVCreate(PETSC_COMM_WORLD, &phaseID_fvm);CHKERRQ(ierr);
  ierr = PetscFVSetSpatialDimension(phaseID_fvm, 3);CHKERRQ(ierr);
  ierr = PetscFVSetNumComponents(phaseID_fvm, MAXFP);CHKERRQ(ierr);
  ierr = PetscFVSetLimiter(phaseID_fvm, nonelimiter);CHKERRQ(ierr);
  ierr = PetscFVSetFromOptions(phaseID_fvm);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) phaseID_fvm,"active phase set");CHKERRQ(ierr);
  ierr = DMSetField(ctx.da_phaseID, 0, NULL, (PetscObject) phaseID_fvm);CHKERRQ(ierr);
  ierr = DMCreateDS(ctx.da_phaseID);CHKERRQ(ierr);
  ierr = DMPlexCreateClosureIndex(ctx.da_phaseID, NULL);

  /* Create finite volume discretisation for material state */
  PetscFV output_fvm;
  ierr = PetscFVCreate(PETSC_COMM_WORLD, &output_fvm);CHKERRQ(ierr);
  ierr = PetscFVSetSpatialDimension(output_fvm, 3);CHKERRQ(ierr);
  ierr = PetscFVSetNumComponents(output_fvm, 1+MAXCP);CHKERRQ(ierr);
  ierr = PetscFVSetLimiter(output_fvm, nonelimiter);CHKERRQ(ierr);
  ierr = PetscFVSetFromOptions(output_fvm);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) output_fvm,"output");CHKERRQ(ierr);
  {
    ierr = PetscFVSetComponentName(output_fvm,0,"phase id");CHKERRQ(ierr);
    for (PetscInt c=1; c < 1+MAXCP; c++) {
        char compname[256];
        if (c <= ctx.nc) {
            ierr = PetscSNPrintf(compname,sizeof(compname),"composition %s",ctx.componentname[c-1]);CHKERRQ(ierr);
        } else {
            ierr = PetscSNPrintf(compname,sizeof(compname),"unknown field %d",c);CHKERRQ(ierr);
        }
        ierr = PetscFVSetComponentName(output_fvm,c,compname);CHKERRQ(ierr);
    }
  }
  ierr = DMSetField(ctx.da_output, 0, NULL, (PetscObject) output_fvm);CHKERRQ(ierr);
  ierr = DMCreateDS(ctx.da_output);CHKERRQ(ierr);
  ierr = DMPlexCreateClosureIndex(ctx.da_output, NULL);

  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Set up problem
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMCreateGlobalVector(ctx.da_solution, &solution);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) solution, "solution");CHKERRQ(ierr);
  ierr = DMCreateLocalVector(ctx.da_phaseID, &ctx.activephaseset);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) ctx.activephaseset, "phase active set");CHKERRQ(ierr);  
  ierr = DMCreateLocalVector(ctx.da_phaseID, &ctx.activephasesuperset);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) ctx.activephasesuperset, "phase active super set");CHKERRQ(ierr);  
  ierr = DMCreateGlobalVector(ctx.da_matstate, &ctx.matstate);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) ctx.matstate, "material state");CHKERRQ(ierr);  
  ierr = SetUpProblem(solution,&ctx);CHKERRQ(ierr);

  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create time stepping object
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = InitializeTS(ctx.da_solution, &ctx, &ts);CHKERRQ(ierr);

  /* initial AMR */
  PetscInt AMRiter;
  for (AMRiter=0; AMRiter<ctx.params.maxnrefine; ++AMRiter) {
      ierr = AMRfunction(ts, solution, refineTag, coarsenTag, &ctx);CHKERRQ(ierr);
  }
  
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Set up and perform time integration 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = TSSetTimeStep(ts,ctx.params.timestep);CHKERRQ(ierr);
  {
      PetscReal currenttime;
      PetscInt  adaptiter, nsteps;

      ierr   = TSSetMaxSteps(ts,ctx.params.amrinterval);CHKERRQ(ierr);
      ierr   = TSSolve(ts,solution);CHKERRQ(ierr);
      ierr   = TSGetSolveTime(ts,&currenttime);CHKERRQ(ierr);
      ierr   = TSGetStepNumber(ts,&nsteps);CHKERRQ(ierr);
      for (adaptiter=0; currenttime < ctx.params.finaltime; adaptiter++) {
          ierr = AMRfunction(ts, solution, refineTag, coarsenTag, &ctx);CHKERRQ(ierr);
          ierr = TSSetStepNumber(ts,nsteps);CHKERRQ(ierr);
          ierr = TSSetTime(ts,currenttime);CHKERRQ(ierr);
          ierr = TSSetTimeStep(ts,ctx.params.timestep);CHKERRQ(ierr);
          ierr = TSSetMaxSteps(ts,nsteps+ctx.params.amrinterval);CHKERRQ(ierr);
          ierr = TSSolve(ts,solution);CHKERRQ(ierr);
          ierr = TSGetSolveTime(ts,&currenttime);CHKERRQ(ierr);
          ierr = TSGetStepNumber(ts,&nsteps);CHKERRQ(ierr);
          ierr = TSGetTimeStep(ts,&ctx.params.timestep);CHKERRQ(ierr);
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
  ierr = VecDestroy(&ctx.matstate);CHKERRQ(ierr);
  ierr = PetscLimiterDestroy(&nonelimiter);CHKERRQ(ierr);
  ierr = PetscFVDestroy(&solution_fvm);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx.da_solution);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx.da_matstate);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx.da_output  );CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}
