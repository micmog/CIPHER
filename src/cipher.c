
static char help[] = "Solves multi phase field equations \n";

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include "init.h"
#include "material.h"
#include "utility.h"
#include "typedef.h"

/* User-defined routines */
extern PetscErrorCode RHSFunction(TS,PetscReal,Vec,Vec,void *);
extern PetscErrorCode PostStep(TS);

/*
 RHSFunction - Evaluates RHS function, F(x)
 */

PetscErrorCode RHSFunction(TS ts,PetscReal ftime,Vec X,Vec F,void *ptr)
{
    PetscErrorCode ierr;
    AppCtx         *user = (AppCtx*) ptr;
    DM             da_solution;
    Vec            localX, matstate0;
    FIELD          ***fdof, ***rhs; 
    STATE          ***matstate;
    F2I            ***slist;
    PetscInt       xs, ys, zs, xm, ym, zm;

    /* Isotropic FD stencil coefficients from https://arxiv.org/abs/1202.3299 */
    PetscReal    fdcoeff0 = -152.0/36.0, fdcoeff1 = 16.0/36.0, fdcoeff2 = 4.0/36.0, fdcoeff3 = 1.0/36.0; 
    PetscReal    fdcoeff[3][3][3] = {{{fdcoeff3, fdcoeff2, fdcoeff3},{fdcoeff2, fdcoeff1, fdcoeff2},{fdcoeff3, fdcoeff2, fdcoeff3}},
                                     {{fdcoeff2, fdcoeff1, fdcoeff2},{fdcoeff1, fdcoeff0, fdcoeff1},{fdcoeff2, fdcoeff1, fdcoeff2}},
                                     {{fdcoeff3, fdcoeff2, fdcoeff3},{fdcoeff2, fdcoeff1, fdcoeff2},{fdcoeff3, fdcoeff2, fdcoeff3}}};
    
    
    PetscFunctionBeginUser;

    ierr = TSGetDM(ts,&da_solution); CHKERRQ(ierr);
    ierr = DMGetLocalVector(da_solution,&localX); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(da_solution,X,INSERT_VALUES,localX); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da_solution,X,INSERT_VALUES,localX); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_solution,localX,&fdof); CHKERRQ(ierr);
    ierr = VecZeroEntries(F); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_solution,F,&rhs); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da_phaseID,user->activephasesuperset,&slist); CHKERRQ(ierr);
    ierr = DMGetGlobalVector(user->da_matstate,&matstate0); CHKERRQ(ierr);
    ierr = VecCopy(user->matstate,matstate0); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da_matstate,matstate0,&matstate); CHKERRQ(ierr);

    ierr = DMDAGetCorners(da_solution,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
    zm+=zs; ym+=ys; xm+=xs; 
    /* Compute function over RHS on active cells */
    for (PetscInt k=zs; k<zm; k++) {
        for (PetscInt j=ys; j<ym; j++) {
            for (PetscInt i=xs; i<xm; i++) {
                if (   ((user->resolution[2] - 1) == k) 
                    || ( 0                        == k)) {
                    memset(rhs[k][j][i].p,0.0,MAXAP*sizeof(PetscReal));
                    memset(rhs[k][j][i].m,0.0,MAXCP*sizeof(PetscReal));
                } else {
                    PetscInt ie=i+1,iw=i-1,jn=j+1,js=j-1,kt=k+1,kb=k-1;
                    
                    /* update material state */
                    Composition(matstate[k][j][i].c,fdof[k][j][i].m,slist[k][j][i].i,user);
                    
                    /* calculate phase interpolants */
                    PetscReal interpolant[slist[k][j][i].i[0]];
                    EvalInterpolant(interpolant,fdof[k][j][i].p,slist[k][j][i].i[0]);
                    
                    /* more than 1 phase active at this point neighbourhood */
                    memset(rhs[k][j][i].p,0.0,MAXAP*sizeof(PetscReal));
                    if (slist[k][j][i].i[0] > 1) {
                        PetscReal phidel[slist[k][j][i].i[0]];
                        uint16_t setintersect[MAXAP], injectionA[MAXAP], injectionB[MAXAP];
                        for (PetscInt kk=kb; kk<=kt; kk++) {
                            for (PetscInt jj=js; jj<=jn; jj++) {
                                for (PetscInt ii=iw; ii<=ie; ii++) {
                                    SetIntersection(setintersect,injectionA,injectionB,slist[k ][j ][i ].i,slist[kk][jj][ii].i);
                                    for (PetscInt g=0; g<setintersect[0]; g++)
                                        phidel[injectionA[g]] += fdcoeff[kk-kb][jj-js][ii-iw]*fdof[kk][jj][ii].p[injectionB[g]];
                                }
                            }
                        }
                    
                        /* capillary driving force */ 
                        PetscReal caplflux[slist[k][j][i].i[0]], caplsource[slist[k][j][i].i[0]];
                        memset(caplflux  ,0,slist[k][j][i].i[0]*sizeof(PetscReal));
                        memset(caplsource,0,slist[k][j][i].i[0]*sizeof(PetscReal));
                        for (PetscInt gk=0; gk<slist[k][j][i].i[0]; gk++) {
                            for (PetscInt gj=gk+1; gj<slist[k][j][i].i[0]; gj++) {
                                PetscInt interfacekj = (PetscInt) user->interfacelist[  slist[k][j][i].i[gk+1]*user->np
                                                                                      + slist[k][j][i].i[gj+1]         ];
                                INTERFACE *currentinterface = &user->interface[interfacekj];
                                caplflux  [gk] -= currentinterface->energy*phidel         [gj];
                                caplflux  [gj] -= currentinterface->energy*phidel         [gk];
                                caplsource[gk] -= currentinterface->energy*fdof[k][j][i].p[gj];
                                caplsource[gj] -= currentinterface->energy*fdof[k][j][i].p[gk];
                            }
                            caplflux  [gk] *= 8.0*user->len/PETSC_PI/PETSC_PI;
                            caplsource[gk] *= 8.0/user->len;
                        }   
                    
                        /* chemical driving force */ 
                        PetscReal chemsource[slist[k][j][i].i[0]];
                        Chemenergy(chemsource,matstate[k][j][i].c,fdof[k][j][i].m,slist[k][j][i].i,user);
                        MatMulInterpolantDerivative(chemsource,fdof[k][j][i].p,slist[k][j][i].i[0]);
                    
                        /* build total RHS and forward solution  */ 
                        PetscReal rhs_unconstrained[slist[k][j][i].i[0]], nactivephases = 0.0;
                        PetscBool active[slist[k][j][i].i[0]];
                        memset(rhs_unconstrained,0,slist[k][j][i].i[0]*sizeof(PetscReal));
                        for (PetscInt gk=0; gk<slist[k][j][i].i[0]; gk++) {
                            for (PetscInt gj=gk+1; gj<slist[k][j][i].i[0]; gj++) {
                                    PetscInt interfacekj = (PetscInt) user->interfacelist[  slist[k][j][i].i[gk+1]*user->np
                                                                                          + slist[k][j][i].i[gj+1]         ];
                                    INTERFACE *currentinterface = &user->interface[interfacekj];
                                    PetscReal val = currentinterface->mobility*(  (caplflux  [gk] - caplflux  [gj])
                                                                                  + (caplsource[gk] - caplsource[gj])
                                                                                  + (chemsource[gj] - chemsource[gk])
                                                                                  * sqrt(interpolant[gk]*interpolant[gj]));
                                    rhs_unconstrained[gk] += val; rhs_unconstrained[gj] -= val;
                            }
                            active[gk] =    (fdof[k][j][i].p[gk] >       TOL && fdof[k][j][i].p  [gk] < 1.0 - TOL)
                                         || (fdof[k][j][i].p[gk] <       TOL && rhs_unconstrained[gk] > 0.0      )
                                         || (fdof[k][j][i].p[gk] > 1.0 - TOL && rhs_unconstrained[gk] < 0.0      );
                            if (active[gk]) nactivephases += 1.0;
                        }

                        for (PetscInt gk=0; gk<slist[k][j][i].i[0]; gk++) {
                            if (active[gk]) {
                                for (PetscInt gj=gk+1; gj<slist[k][j][i].i[0]; gj++) {
                                    if (active[gj]) {
                                        PetscInt interfacekj = (PetscInt) user->interfacelist[  slist[k][j][i].i[gk+1]*user->np
                                                                                              + slist[k][j][i].i[gj+1]         ];
                                        INTERFACE *currentinterface = &user->interface[interfacekj];
                                        PetscReal val = currentinterface->mobility*(  (caplflux  [gk] - caplflux  [gj])
                                                                                      + (caplsource[gk] - caplsource[gj])
                                                                                      + (chemsource[gj] - chemsource[gk])
                                                                                      * sqrt(interpolant[gk]*interpolant[gj]));
                                        rhs[k][j][i].p[gk] += val; rhs[k][j][i].p[gj] -= val;
                                    }
                                }
                                rhs[k][j][i].p[gk] /= nactivephases;
                            }        
                        }
                    }
                    
                    /* more than 1 component active at this point */
                    memset(rhs[k][j][i].m,0.0,MAXCP*sizeof(PetscReal));
                    if (user->nc > 1) {
                        /* gather neighbourhood field values */
                        PetscReal chempotdel[user->nc-1];
                        memset(chempotdel,0,(user->nc-1)*sizeof(PetscReal));
                        for (PetscInt kk=kb; kk<=kt; kk++) {
                            for (PetscInt jj=js; jj<=jn; jj++) {
                                for (PetscInt ii=iw; ii<=ie; ii++) {
                                    for (PetscInt c=0; c<user->nc-1; c++) {
                                        chempotdel[c] += fdcoeff[kk-kb][jj-js][ii-iw]*fdof[kk][jj][ii].m[c];
                                    }
                                }
                            }
                        }
                      
                        /* calculate mobility matrix (volume-fixed frame of reference) */
                        PetscReal mobilitycv[user->nc-1];
                        CompositionMobility(mobilitycv,matstate[k][j][i].c,interpolant,slist[k][j][i].i,user);
                      
                        /* calculate avg composition rate */
                        PetscReal cavgdot[user->nc-1], cvgdot_phase[slist[k][j][i].i[0]];
                        for (PetscInt c=0; c<user->nc-1; c++) {
                            cavgdot[c] = mobilitycv[c]*chempotdel[c];
                            for (PetscInt g =0; g<slist[k][j][i].i[0];  g++) cvgdot_phase[g] = matstate[k][j][i].c[g*user->nc+c];
                            MatMulInterpolantDerivative(cvgdot_phase,fdof[k][j][i].p,slist[k][j][i].i[0]);
                            for (PetscInt g =0; g<slist[k][j][i].i[0];  g++) cavgdot[c] -= rhs[k][j][i].p[g]*cvgdot_phase[g];
                        }
                      
                        /* calculate composition tangent wrt chemical potential */
                        PetscReal dcdm[(user->nc-1)*(user->nc-1)], dmdc[(user->nc-1)*(user->nc-1)];
                        CompositionTangent(dcdm ,matstate[k][j][i].c,interpolant,slist[k][j][i].i,user);
                        Invertmatrix(dmdc,dcdm);
                      
                        /* update chemical potential */
                        for (PetscInt cj=0; cj<user->nc-1; cj++) {
                            for (PetscInt ci=0; ci<user->nc-1; ci++) {
                                rhs[k][j][i].m[cj] += dmdc[cj*(user->nc-1)+ci]*cavgdot[ci];
                            }
                        }
                    }
                }      
            }
        }
    }       
    ierr = DMDAVecRestoreArray(user->da_matstate,matstate0,&matstate); CHKERRQ(ierr);
    ierr = DMRestoreGlobalVector(user->da_matstate,&matstate0); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->da_phaseID,user->activephasesuperset,&slist); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_solution,F,&rhs); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_solution,localX,&fdof); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(da_solution,&localX); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

/*
 PostStep - Update active phases and write output
 */
PetscErrorCode PostStep(TS ts)
{
    PetscErrorCode    ierr;
    AppCtx            *user;
    DM                da_solution;
    Vec               solution, gactivephaseset, gactivephasesuperset;
    FIELD             ***fdof;
    STATE             ***matstate;
    F2I               ***alist, ***slist, ***galist, ***gslist;
    PetscInt          xs, ys, zs, xm, ym, zm;    
  
    PetscFunctionBeginUser;    
    
    ierr = TSGetDM(ts,&da_solution); CHKERRQ(ierr);
    ierr = TSGetSolution(ts, &solution);CHKERRQ(ierr);    
    ierr = TSGetApplicationContext(ts,&user);CHKERRQ(ierr);
    ierr = DMDAGetCorners(da_solution,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
    zm+=zs; ym+=ys; xm+=xs; 
    
    /* Determine active phase set */
    ierr = DMDAVecGetArray(da_solution,solution,&fdof); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da_matstate,user->matstate,&matstate); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da_phaseID,user->activephasesuperset,&slist); CHKERRQ(ierr);
    ierr = DMGetGlobalVector(user->da_phaseID,&gactivephaseset); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da_phaseID,gactivephaseset,&galist); CHKERRQ(ierr);
    for (PetscInt k=zs; k<zm; k++) {
        for (PetscInt j=ys; j<ym; j++) {
            for (PetscInt i=xs; i<xm; i++) {
                /* update material state */
                PetscReal phiUP[slist[k][j][i].i[0]], interpolant[slist[k][j][i].i[0]];
                memcpy(phiUP,fdof[k][j][i].p,slist[k][j][i].i[0]*sizeof(PetscReal));
                SimplexProjection(fdof[k][j][i].p,phiUP,slist[k][j][i].i[0]);
                EvalInterpolant(interpolant,fdof[k][j][i].p,slist[k][j][i].i[0]);
                Composition(matstate[k][j][i].c,fdof[k][j][i].m,slist[k][j][i].i,user);

                /* update active set */
                galist[k][j][i].i[0] = 0;
                for (PetscInt g=0; g<slist[k][j][i].i[0];  g++) {
                    if (interpolant[g] > TOL) {
                        galist[k][j][i].i[++galist[k][j][i].i[0]] = slist[k][j][i].i[g+1];
                    }
                }
            }
        }
    }        
    ierr = DMDAVecRestoreArray(user->da_phaseID,gactivephaseset,&galist); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(user->da_phaseID,gactivephaseset,INSERT_VALUES,user->activephaseset); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->da_phaseID,gactivephaseset,INSERT_VALUES,user->activephaseset); CHKERRQ(ierr);
    ierr = DMRestoreGlobalVector(user->da_phaseID,&gactivephaseset); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->da_phaseID,user->activephasesuperset,&slist); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da_matstate,user->matstate,&matstate); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_solution,solution,&fdof); CHKERRQ(ierr);

    /* Determine active phase super set, reorder dofs in solution, update composition, initialize new phase frac & comp */
    ierr = DMDAVecGetArray(da_solution,solution,&fdof); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da_matstate,user->matstate,&matstate); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da_phaseID,user->activephaseset,&alist); CHKERRQ(ierr);
    ierr = DMGetGlobalVector(user->da_phaseID,&gactivephasesuperset); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da_phaseID,gactivephasesuperset,&gslist); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da_phaseID,user->activephasesuperset,&slist); CHKERRQ(ierr);
    for (PetscInt k=zs; k<zm; k++) {
        for (PetscInt j=ys; j<ym; j++) {
            for (PetscInt i=xs; i<xm; i++) {
                /* update active phase superset */
                PetscInt ie=i+1,iw=i-1,jn=j+1,js=j-1,kt=k+1,kb=k-1;
                uint16_t setunion[MAXAP], setintersection[MAXAP], injectionA[MAXAP], injectionB[MAXAP];
                gslist[k][j][i].i[0] = 0;
                for (PetscInt kk=kb; kk<=kt; kk++) {
                    for (PetscInt jj=js; jj<=jn; jj++) {
                        for (PetscInt ii=iw; ii<=ie; ii++) {
                            SetUnion(setunion,injectionA,injectionB,alist[kk][jj][ii].i,gslist[k][j][i].i);
                            memcpy(gslist[k][j][i].i,setunion,2*MAXAP);          
                        }
                    }
                }
                
                /* initialize new phase states */
                PetscReal phiSS[gslist[k][j][i].i[0]], compSS[gslist[k][j][i].i[0]*user->nc];
                for (PetscInt g=0; g<gslist[k][j][i].i[0];  g++) {
                    phiSS[g] = 0.0;
                    MATERIAL *currentmaterial = &user->material[user->phasematerialmapping[gslist[k][j][i].i[g+1]]];
                    memcpy(&compSS[g*user->nc],currentmaterial->c0,user->nc*sizeof(PetscReal));
                }
                
                /* reorder dofs to new active phase superset */
                SetIntersection(setintersection,injectionA,injectionB,slist[k][j][i].i,gslist[k][j][i].i);
                for (PetscInt g=0; g<setintersection[0];  g++) {
                    phiSS[injectionB[g]] = fdof[k][j][i].p[injectionA[g]];
                    memcpy(&compSS[injectionB[g]*user->nc],&matstate[k][j][i].c[injectionA[g]*user->nc],user->nc*sizeof(PetscReal));
                }
                memcpy(fdof[k][j][i].p,phiSS,gslist[k][j][i].i[0]*sizeof(PetscReal));
                memcpy(matstate[k][j][i].c,compSS,gslist[k][j][i].i[0]*user->nc*sizeof(PetscReal));
            }
        }
    }        
    ierr = DMDAVecRestoreArray(user->da_phaseID,user->activephasesuperset,&slist); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da_phaseID,gactivephasesuperset,&gslist); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(user->da_phaseID,gactivephasesuperset,INSERT_VALUES,user->activephasesuperset); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->da_phaseID,gactivephasesuperset,INSERT_VALUES,user->activephasesuperset); CHKERRQ(ierr);
    ierr = DMRestoreGlobalVector(user->da_phaseID,&gactivephasesuperset); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->da_phaseID,user->activephaseset,&alist); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_solution,solution,&fdof); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da_matstate,user->matstate,&matstate); CHKERRQ(ierr);
    ierr = TSSetSolution(ts,solution); CHKERRQ(ierr);
    
    /* write output */
    ierr = TSGetStepNumber(ts, &user->step);CHKERRQ(ierr);
    if (user->step%user->outputfreq == 0) {
       Vec               Xout[user->nc+1];
       O_DOFS            xout[user->nc+1];
       char              name[256];
       PetscViewer       viewer;
       
       for (PetscInt c=0; c<user->nc+1; c++) {
           ierr = DMGetGlobalVector(user->da_output,&Xout[c]); CHKERRQ(ierr);
           ierr = DMDAVecGetArray(user->da_output,Xout[c],&xout[c].o); CHKERRQ(ierr);
       }
       ierr = DMDAVecGetArrayRead(da_solution,solution,&fdof); CHKERRQ(ierr);
       ierr = DMDAVecGetArrayRead(user->da_matstate,user->matstate,&matstate); CHKERRQ(ierr);
       ierr = DMDAVecGetArrayRead(user->da_phaseID,user->activephasesuperset,&slist); CHKERRQ(ierr);
       for (PetscInt k=zs; k<zm; k++) {
           for (PetscInt j=ys; j<ym; j++) {
               for (PetscInt i=xs; i<xm; i++) {
                   PetscReal max = -LARGE, interpolant[slist[k][j][i].i[0]];
                   EvalInterpolant(interpolant,fdof[k][j][i].p,slist[k][j][i].i[0]);
                   for (PetscInt g=0; g<slist[k][j][i].i[0]; g++) {
                       if (interpolant[g] > max) {
                          max = interpolant[g];
                          xout[0].o[k][j][i] = (PetscReal) slist[k][j][i].i[g+1];
                       }
                   }
                   for (PetscInt c=0; c<user->nc; c++) {
                       xout[c+1].o[k][j][i] = 0.0;
                       for (PetscInt g=0; g<slist[k][j][i].i[0]; g++)
                           xout[c+1].o[k][j][i] += interpolant[g]*matstate[k][j][i].c[g*user->nc+c];        
                   }
               }
           }
       }       
       ierr = DMDAVecRestoreArrayRead(user->da_phaseID,user->activephasesuperset,&slist); CHKERRQ(ierr);
       ierr = DMDAVecRestoreArrayRead(user->da_matstate,user->matstate,&matstate); CHKERRQ(ierr);
       ierr = DMDAVecRestoreArrayRead(da_solution,solution,&fdof); CHKERRQ(ierr);
       sprintf(name, "%s_%d.vtr",user->outfile,user->step);
       PetscPrintf(PETSC_COMM_WORLD,"writing output to %s\n",name);
       ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,name,FILE_MODE_WRITE,&viewer);
       ierr = DMDAVecRestoreArray(user->da_output,Xout[0],&xout[0].o); CHKERRQ(ierr);
       ierr = PetscObjectSetName((PetscObject) Xout[0],"phase");
       ierr = VecView(Xout[0],viewer);
       ierr = DMRestoreGlobalVector(user->da_output,&Xout[0]); CHKERRQ(ierr);
       for (PetscInt c=0; c<user->nc; c++) {
           ierr = DMDAVecRestoreArray(user->da_output,Xout[c+1],&xout[c+1].o); CHKERRQ(ierr);
           sprintf(name,"composition %s",user->componentname[c]);
           ierr = PetscObjectSetName((PetscObject) Xout[c+1],name);
           ierr = VecView(Xout[c+1],viewer);
           ierr = DMRestoreGlobalVector(user->da_output,&Xout[c+1]); CHKERRQ(ierr);
       }
       ierr = PetscViewerDestroy(&viewer);
    }
    PetscFunctionReturn(0);
}

int main(int argc,char **args)
{
  /* user-defined work context */
  AppCtx         user;         
  /* solution grid */
  DM             da_solution;
  Vec            solution;
  /* time stepping context */
  TS             ts;
  /* numerical parameters */
  PetscReal      finaltime;
  PetscErrorCode ierr;
      
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Initialize program
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  assert(MAXAP == MAXTP*sizeof(PetscReal)/sizeof(uint16_t));
  ierr = PetscInitialize(&argc,&args,(char*)0,help); CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Initialize problem parameters
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SetUpGeometry(&user); CHKERRQ(ierr);
  ierr = SetUpInterface(&user); CHKERRQ(ierr);
  utility_init(&user);
  material_init(&user);
  finaltime = 10.0;
  ierr = PetscOptionsGetReal(NULL,NULL,"-finaltime",&finaltime,NULL); CHKERRQ(ierr);
  user.dt = 0.5;
  ierr = PetscOptionsGetReal(NULL,NULL,"-mintimestep",&user.dt,NULL); CHKERRQ(ierr);
  user.dtmax = 10.0;
  ierr = PetscOptionsGetReal(NULL,NULL,"-maxtimestep",&user.dtmax,NULL); CHKERRQ(ierr);
  user.len = 4.0;
  ierr = PetscOptionsGetReal(NULL,NULL,"-l",&user.len,NULL); CHKERRQ(ierr);
  user.ptol = TOL;
  ierr = PetscOptionsGetReal(NULL,NULL,"-ptol",&user.ptol,NULL); CHKERRQ(ierr);
  user.ctol = TOL;
  ierr = PetscOptionsGetReal(NULL,NULL,"-ctol",&user.ctol,NULL); CHKERRQ(ierr);
  user.outputfreq = 1;
  ierr = PetscOptionsGetInt(NULL,NULL,"-outfreq",&user.outputfreq,NULL); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL,NULL,"-outfile",user.outfile,128,NULL);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create distributed array (DMDA) to manage parallel grid and vectors
   for the multi-phase field PDE
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMDACreate3d(PETSC_COMM_WORLD, 
                      DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,
                      DMDA_STENCIL_BOX, 
                      user.resolution[0],user.resolution[1],user.resolution[2],
                      PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
                      MAXAP+MAXCP,1,
                      NULL,NULL,NULL,
                      &da_solution); CHKERRQ(ierr);
  ierr = DMSetFromOptions(da_solution); CHKERRQ(ierr);
  ierr = DMSetUp(da_solution); CHKERRQ(ierr);
  ierr = DMSetApplicationContext(da_solution,&user); CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(da_solution,&solution); CHKERRQ(ierr);
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create aux distributed array (DMDA) to manage bitwise data
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMDACreate3d(PETSC_COMM_WORLD, 
                      DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,
                      DMDA_STENCIL_BOX, 
                      user.resolution[0],user.resolution[1],user.resolution[2],
                      PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
                      MAXTP,1,
                      NULL,NULL,NULL,
                      &user.da_phaseID); CHKERRQ(ierr);
  ierr = DMSetUp(user.da_phaseID); CHKERRQ(ierr);
  ierr = DMCreateLocalVector(user.da_phaseID,&user.activephaseset); CHKERRQ(ierr);
  ierr = DMCreateLocalVector(user.da_phaseID,&user.activephasesuperset); CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create aux distributed array (DMDA) to manage bitwise data
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMDACreate3d(PETSC_COMM_WORLD, 
                      DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,
                      DMDA_STENCIL_BOX, 
                      user.resolution[0],user.resolution[1],user.resolution[2],
                      PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
                      MAXAP*MAXCP,1,
                      NULL,NULL,NULL,
                      &user.da_matstate); CHKERRQ(ierr);
  ierr = DMSetUp(user.da_matstate); CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(user.da_matstate,&user.matstate); CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create aux distributed array (DMDA) to manage output
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMDACreate3d(PETSC_COMM_WORLD, 
                      DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,
                      DMDA_STENCIL_BOX, 
                      user.resolution[0],user.resolution[1],user.resolution[2],
                      PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
                      1,1,
                      NULL,NULL,NULL,
                      &user.da_output); CHKERRQ(ierr);
  ierr = DMSetUp(user.da_output); CHKERRQ(ierr);

  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Initialize solution
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = TSCreate(PETSC_COMM_WORLD,&ts); CHKERRQ(ierr);
  ierr = TSSetApplicationContext(ts, &user);CHKERRQ(ierr);
  ierr = TSSetProblemType(ts,TS_NONLINEAR); CHKERRQ(ierr);
  ierr = TSSetRHSFunction(ts,NULL,RHSFunction,&user); CHKERRQ(ierr);
  ierr = TSSetMaxTime(ts,finaltime); CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP); CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts,user.dt); CHKERRQ(ierr);
  ierr = TSSetPostStep(ts,PostStep); CHKERRQ(ierr);
  ierr = TSSetDM(ts,da_solution); CHKERRQ(ierr);
  ierr = TSSetSolution(ts,solution); CHKERRQ(ierr);
  ierr = TSSetType(ts,TSRK); CHKERRQ(ierr);
  ierr = TSSetTolerances(ts,user.ptol,NULL,user.ctol,NULL); CHKERRQ(ierr);
  
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Set up and perform time integration 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = TSSetFromOptions(ts); CHKERRQ(ierr);  
  ierr = SetUpProblem(da_solution,solution,&user); CHKERRQ(ierr);
  ierr = TSSolve(ts,solution); CHKERRQ(ierr);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Free work space.  All PETSc objects should be destroyed when they
   are no longer needed.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = VecDestroy(&solution); CHKERRQ(ierr);
  ierr = DMDestroy(&da_solution); CHKERRQ(ierr);
  ierr = VecDestroy(&user.activephasesuperset); CHKERRQ(ierr);
  ierr = VecDestroy(&user.activephaseset); CHKERRQ(ierr);
  ierr = DMDestroy(&user.da_phaseID); CHKERRQ(ierr);
  ierr = VecDestroy(&user.matstate); CHKERRQ(ierr);
  ierr = DMDestroy(&user.da_matstate); CHKERRQ(ierr);
  ierr = DMDestroy(&user.da_output); CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}
