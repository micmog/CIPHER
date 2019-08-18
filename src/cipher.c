
static char help[] = "Solves multi phase field equations \n";

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <petscdm.h>
#include <petscdmda.h>
#include "init.h"
#include "material.h"
#include "utility.h"
#include "typedef.h"

/* User-defined routines */
extern PetscErrorCode EulerStep(DM, Vec, AppCtx*);
extern PetscErrorCode PostStep(DM, Vec, AppCtx*);

/*
 EulerStep - Evaluates Euler step
 */

PetscErrorCode EulerStep(DM da,Vec X,AppCtx *user)
{
    PetscErrorCode ierr;
    Vec            localX, globalX, globalXe, globalXdot, gactivephases, globalComp;
    F_DOFS         ***fdof, ***gfdof, ***fdot, ***lte;
    C_DOFS         ***cdof; 
    F2I            ***glist, ***nalist, ***alist;
    PetscInt       xs, ys, zs, xm, ym, zm;
    PetscScalar    deltaphi[MAXAP]; 
    uint16_t       superset[MAXAP];
    PetscScalar    fdcoeff[3][3][3] = {{{1.0/36.0, 4.0/36.0, 1.0/36.0},{4.0/36.0, 16.0/36.0, 4.0/36.0},{1.0/36.0, 4.0/36.0, 1.0/36.0}},
                                       {{4.0/36.0, 16.0/36.0, 4.0/36.0},{16.0/36.0, -152.0/36.0, 16.0/36.0},{4.0/36.0, 16.0/36.0, 4.0/36.0}},
                                       {{1.0/36.0, 4.0/36.0, 1.0/36.0},{4.0/36.0, 16.0/36.0, 4.0/36.0},{1.0/36.0, 4.0/36.0, 1.0/36.0}}};
    
    PetscFunctionBeginUser;

    ierr = DMGetLocalVector(da,&localX); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(da,X,INSERT_VALUES,localX); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da,X,INSERT_VALUES,localX); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da,localX,&fdof); CHKERRQ(ierr);
    ierr = DMGetGlobalVector(user->daCmp,&globalComp); CHKERRQ(ierr);
    ierr = VecCopy(user->cvec,globalComp); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->daCmp,globalComp,&cdof); CHKERRQ(ierr);
    ierr = DMGetGlobalVector(user->daIdx,&gactivephases); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->daIdx,gactivephases,&glist); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->daIdx,user->activeneighbourphases,&nalist); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->daIdx,user->activephases,&alist); CHKERRQ(ierr);
    ierr = DMGetGlobalVector(da,&globalXe  ); CHKERRQ(ierr);
    ierr = DMGetGlobalVector(da,&globalX   ); CHKERRQ(ierr);
    ierr = DMGetGlobalVector(da,&globalXdot); CHKERRQ(ierr);
    ierr = VecCopy(user->rhs0,globalXdot); CHKERRQ(ierr);
    ierr = VecZeroEntries(globalXe); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,globalXe  ,&lte  ); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,globalX   ,&gfdof); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,globalXdot,&fdot ); CHKERRQ(ierr);

    ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
    zm+=zs; ym+=ys; xm+=xs; 
    /* Compute function over RHS on active cells */
    int err = 1; 
    for (PetscInt k=zs; k<zm && err; k++) {
        for (PetscInt j=ys; j<ym && err; j++) {
            for (PetscInt i=xs; i<xm && err; i++) {
                /* more than 1 phase active at this point neighbourhood */
                PetscInt ie=i+1,iw=i-1,jn=j+1,js=j-1,kt=k+1,kb=k-1;
                memset(deltaphi,0,MAXAP*sizeof(PetscScalar));
                if (   ((alist[k][j][i].i[0] + nalist[k][j][i].i[0]) <= 1) 
                    || ((user->resolution[2] - 1                   ) == k) 
                    || ( 0                                           == k)) {
                    memcpy(glist[k][j][i].i,alist[k][j][i].i,MAXAP*sizeof(uint16_t   ));
                    memcpy(gfdof[k][j][i].p, fdof[k][j][i].p,MAXAP*sizeof(PetscScalar));
                } else {   
                    /* gather neighbourhood field values */
                    superset[0] = union_uint16(&nalist[k][j][i].i[1],nalist[k][j][i].i[0],
                                               & alist[k][j][i].i[1], alist[k][j][i].i[0], 
                                               &superset[1]);
                    PetscScalar phi[user->np], phidel[user->np], phidot[user->np], phir[superset[0]];
                    memset(phi   ,0,user->np*sizeof(PetscScalar));
                    memset(phidot,0,user->np*sizeof(PetscScalar));
                    memset(phidel,0,user->np*sizeof(PetscScalar));
                    for (PetscInt kk=kb; kk<=kt; kk++) {
                        for (PetscInt jj=js; jj<=jn; jj++) {
                            for (PetscInt ii=iw; ii<=ie; ii++) {
                                for (PetscInt g=0; g<alist[kk][jj][ii].i[0]; g++) {
                                    phidel[alist[kk][jj][ii].i[g+1]] += 
                                        fdcoeff[kk-kb][jj-js][ii-iw]*fdof[kk][jj][ii].p[g];
                                }
                            }
                        }
                    }
                    for (PetscInt g=0; g<alist[k ][j ][i ].i[0]; g++) {
                        phi   [alist[k ][j ][i ].i[g+1]]  = fdof[k ][j ][i ].p[g];
                        phidot[alist[k ][j ][i ].i[g+1]]  = fdot[k ][j ][i ].p[g];
                    }
                    for (PetscInt g=0; g<superset[0]; g++) phir[g] = phi[superset[g+1]];
                    
                    PetscScalar compos[user->np*user->nc], compag[superset[0]*user->nc];
                    for (PetscInt g=0; g<nalist[k ][j ][i ].i[0]; g++)
                        memcpy(&compos[nalist[k][j][i].i[g+1]*user->nc],
                               user->material[user->phasematerialmapping[nalist[k][j][i].i[g+1]]].c0,
                               user->nc*sizeof(PetscScalar)); // ToDo: use neighbour composition instead of c0
                    for (PetscInt g=0; g< alist[k ][j ][i ].i[0]; g++)
                        memcpy(&compos[ alist[k][j][i].i[g+1]*user->nc],
                               &cdof[k][j][i].c[g*user->nc],user->nc*sizeof(PetscScalar));
                    for (PetscInt g=0; g<superset[0]; g++)
                        memcpy(&compag[g*user->nc],&compos[superset[g+1]*user->nc],user->nc*sizeof(PetscScalar));

                    /* calculate interpolants */
                    PetscScalar interpolant[superset[0]];
                    EvalInterpolant(interpolant,phir,superset[0]);
                    /* capillary driving force */ 
                    PetscScalar caplflux[superset[0]], caplsource[superset[0]];
                    memset(caplflux  ,0,superset[0]*sizeof(PetscScalar));
                    memset(caplsource,0,superset[0]*sizeof(PetscScalar));
                    for (PetscInt gk=0; gk<superset[0]; gk++) {
                        for (PetscInt gj=gk+1; gj<superset[0]; gj++) {
                            PetscInt interfacekj = (PetscInt) user->interfacelist[  superset[gk+1]*user->np
                                                                                  + superset[gj+1]         ];
                            INTERFACE *currentinterface = &user->interface[interfacekj];
                            caplflux  [gk] -= currentinterface->energy*phidel[superset[gj+1]];
                            caplflux  [gj] -= currentinterface->energy*phidel[superset[gk+1]];
                            caplsource[gk] -= currentinterface->energy*phir  [         gj   ];
                            caplsource[gj] -= currentinterface->energy*phir  [         gk   ];
                        }
                        caplflux  [gk] *= 8.0*user->len/PETSC_PI/PETSC_PI;
                        caplsource[gk] *= 8.0/user->len;
                    }   
                    /* chemical driving force */ 
                    PetscScalar chemsource[superset[0]], chemsourcer[superset[0]];
                    Chemenergy(chemsourcer,compag,fdof[k][j][i].m,superset,user);
                    MatMulInterpolantDerivative(chemsource,chemsourcer,phir,superset[0]);
                    /* build total RHS and forward solution  */ 
                    PetscScalar rhs[superset[0]], phiunprojected[superset[0]], phiprojected[superset[0]];
                    memset(rhs,0,superset[0]*sizeof(PetscScalar));
                    for (PetscInt gk=0; gk<superset[0]; gk++) {
                        for (PetscInt gj=gk+1; gj<superset[0]; gj++) {
                            PetscInt interfacekj = (PetscInt) user->interfacelist[  superset[gk+1]*user->np
                                                                                  + superset[gj+1]         ];
                            INTERFACE *currentinterface = &user->interface[interfacekj];
                            PetscScalar val = currentinterface->mobility*(  (caplflux  [gk] - caplflux  [gj])
                                                                          + (caplsource[gk] - caplsource[gj])
                                                                          + (chemsource[gj] - chemsource[gk])
                                                                          * sqrt(interpolant[gk]*interpolant[gj]));
                            rhs[gk] += val; rhs[gj] -= val;
                        }
                        rhs[gk] /= ((PetscScalar) superset[0]);
                        phiunprojected[gk] = phir[gk] + user->dt*rhs[gk];
                    }
                    /* project solution and update active list */
                    SimplexProjection(phiprojected ,phiunprojected ,superset[0]);
                    glist[k][j][i].i[0] = 0;
                    for (PetscInt g=0; g<superset[0];g++) {
                        if (phiprojected[g] > TOL || rhs[g] > TOL) {
                            lte  [k][j][i].p[  glist[k][j][i].i[0]] = 
                                (phiprojected[g] - phir[g] - user->dt*phidot[superset[g+1]])/user->ptol;
                            fdot [k][j][i].p[  glist[k][j][i].i[0]] = (phiprojected[g] - phir[g])/user->dt;
                            gfdof[k][j][i].p[  glist[k][j][i].i[0]] =  phiprojected[g];
                            deltaphi        [  glist[k][j][i].i[0]] =  phiprojected[g] - phir[g];
                            memcpy(&cdof[k][j][i].c[glist[k][j][i].i[0]*user->nc],
                                   &compag[g*user->nc],user->nc*sizeof(PetscScalar));
                            glist[k][j][i].i[++glist[k][j][i].i[0]] = superset[g+1];
                        }
                    }
                }
                /* update compositions */
                if (   ( user->nc               == 1)
                    || ((user->resolution[2]-1) == k) 
                    || ( 0                      == k)) {
                    memcpy(gfdof[k][j][i].m,fdof[k][j][i].m,(user->nc-1)*sizeof(PetscScalar));
                } else {   
                    /* gather neighbourhood field values */
                    PetscScalar chempotdel[user->nc-1];
                    memset(chempotdel,0,(user->nc-1)*sizeof(PetscScalar));
                    for (PetscInt kk=kb; kk<=kt; kk++) {
                        for (PetscInt jj=js; jj<=jn; jj++) {
                            for (PetscInt ii=iw; ii<=ie; ii++) {
                                for (PetscInt c=0; c<user->nc-1; c++) {
                                    chempotdel[c] += fdcoeff[kk-kb][jj-js][ii-iw]*fdof[kk][jj][ii].m[c];
                                }
                            }
                        }
                    }
                    /* calculate interpolants */
                    PetscScalar interpolant[glist[k][j][i].i[0]];
                    EvalInterpolant(interpolant,gfdof[k][j][i].p,glist[k][j][i].i[0]);
                    /* calculate mobility matrix (volume-fixed frame of reference) */
                    PetscScalar mobilitycv[user->nc-1];
                    CompositionMobility(mobilitycv,cdof[k][j][i].c,interpolant,glist[k][j][i].i,user);
                    /* calculate composition rate */
                    PetscScalar cavgdot[user->nc-1], deltac[glist[k][j][i].i[0]], ctemp[glist[k][j][i].i[0]];
                    for (PetscInt c=0; c<user->nc-1; c++) {
                        cavgdot[c] = mobilitycv[c]*chempotdel[c];
                        lte  [k][j][i].m[c] = user->dt*(cavgdot[c] - fdot[k][j][i].m[c])/user->ctol;
                        fdot [k][j][i].m[c] = cavgdot[c];
                        cavgdot[c] *= user->dt;
                        for (PetscInt g =0; g<glist[k][j][i].i[0];  g++) ctemp[g] = cdof[k][j][i].c[g*user->nc+c];
                        MatMulInterpolantDerivative(deltac,ctemp,gfdof[k][j][i].p,glist[k][j][i].i[0]);
                        for (PetscInt g =0; g<glist[k][j][i].i[0];  g++) cavgdot[c] -= deltac[g]*deltaphi[g];
                    }
                    /* calculate composition tangent wrt chemical potential */
                    PetscScalar dcdm[(user->nc-1)*(user->nc-1)], dmdc[(user->nc-1)*(user->nc-1)];
                    CompositionTangent(dcdm ,cdof[k][j][i].c,interpolant,glist[k][j][i].i,user);
                    Invertmatrix(dmdc,dcdm);
                    /* update chemical potential */
                    for (PetscInt cj=0; cj<user->nc; cj++) {
                        gfdof[k][j][i].m[cj] = fdof[k][j][i].m[cj];
                        for (PetscInt ci=0; ci<user->nc; ci++) {
                            gfdof[k][j][i].m[cj] += dmdc[cj*(user->nc-1)+ci]*cavgdot[ci];
                        }
                    }
                    /* get chemical part of diffusion potential */
                    PetscScalar interfacepotential[user->nc];
                    memset(interfacepotential,0,user->nc*sizeof(PetscScalar));
                    for (PetscInt c=0; c<user->nc; c++) {
                        for (PetscInt gk=0; gk<glist[k][j][i].i[0]; gk++) {
                            for (PetscInt gj=gk+1; gj<glist[k][j][i].i[0]; gj++) {
                                PetscInt interfacekj = (PetscInt) user->interfacelist[  glist[k][j][i].i[gk+1]*user->np
                                                                                      + glist[k][j][i].i[gj+1]         ];
                                INTERFACE *currentinterface = &user->interface[interfacekj];
                                interfacepotential[c] += sqrt(interpolant[gk]*interpolant[gj])*currentinterface->potential[c];
                            }
                        }
                    }    
                    PetscScalar chempot[user->nc-1];
                    memcpy(chempot,gfdof[k][j][i].m,(user->nc-1)*sizeof(PetscScalar));
                    for (PetscInt c=0; c<user->nc-1; c++) chempot[c] -= (interfacepotential[c] - interfacepotential[user->nc-1]);
                    /* update composition for given chemical potential */
                    err = Composition(cdof[k][j][i].c,chempot,glist[k][j][i].i,user);
                }
            }
        }
    }       
    ierr = DMDAVecRestoreArray(da,globalXdot,&fdot ); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,globalX   ,&gfdof); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,globalXe  ,&lte  ); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->daCmp,globalComp,&cdof); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->daIdx,user->activephases,&alist); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->daIdx,user->activeneighbourphases,&nalist); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->daIdx,gactivephases,&glist); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da,localX,&fdof); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(da,&localX); CHKERRQ(ierr);
    ierr = MPI_Allreduce(MPI_IN_PLACE,(void *) &err,1,MPI_INT,MPI_LAND,PETSC_COMM_WORLD);

    /* update time step */
    PetscScalar error;
    ierr = VecNorm(globalXe,NORM_INFINITY,&error); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"Step = %d, dt = %g, time = %g, LTE = %g\n",user->step,user->dt,user->time,error);
    if (error > 10.0) {
        /* cut back time step... */
        PetscPrintf(PETSC_COMM_WORLD,"Step rejected\n");
        user->dt /= 2.0;
    } else if (!err) {
        PetscPrintf(PETSC_COMM_WORLD,"Step rejected due to exception\n");
        user->dt /= 2.0;
    } else {
        /* communicate updated neighbourhood bits across processors... */
        ierr = DMGlobalToLocalBegin(user->daIdx,gactivephases,INSERT_VALUES,user->activephases); CHKERRQ(ierr);
        ierr = DMGlobalToLocalEnd(user->daIdx,gactivephases,INSERT_VALUES,user->activephases); CHKERRQ(ierr);
        ierr = DMDAVecGetArray(user->daIdx,user->activephases,&alist); CHKERRQ(ierr);
        ierr = DMDAVecGetArray(user->daIdx,user->activeneighbourphases,&nalist); CHKERRQ(ierr);
        for (PetscInt k=zs; k<zm; k++) {
            for (PetscInt j=ys; j<ym; j++) {
                for (PetscInt i=xs; i<xm; i++) {
                    nalist[k][j][i].i[0] = 0;
                    superset[0] = union_uint16(&nalist[k][j][i].i[1],nalist[k][j][i].i[0],&alist[k  ][j  ][i+1].i[1],alist[k  ][j  ][i+1].i[0],&superset[1]);
                    memcpy(nalist[k][j][i].i,superset,2*MAXAP);          
                    superset[0] = union_uint16(&nalist[k][j][i].i[1],nalist[k][j][i].i[0],&alist[k  ][j  ][i-1].i[1],alist[k  ][j  ][i-1].i[0],&superset[1]);
                    memcpy(nalist[k][j][i].i,superset,2*MAXAP);          
                    superset[0] = union_uint16(&nalist[k][j][i].i[1],nalist[k][j][i].i[0],&alist[k  ][j+1][i  ].i[1],alist[k  ][j+1][i  ].i[0],&superset[1]);
                    memcpy(nalist[k][j][i].i,superset,2*MAXAP);          
                    superset[0] = union_uint16(&nalist[k][j][i].i[1],nalist[k][j][i].i[0],&alist[k  ][j-1][i  ].i[1],alist[k  ][j-1][i  ].i[0],&superset[1]);
                    memcpy(nalist[k][j][i].i,superset,2*MAXAP);          
                    superset[0] = union_uint16(&nalist[k][j][i].i[1],nalist[k][j][i].i[0],&alist[k+1][j  ][i  ].i[1],alist[k+1][j  ][i  ].i[0],&superset[1]);
                    memcpy(nalist[k][j][i].i,superset,2*MAXAP);          
                    superset[0] = union_uint16(&nalist[k][j][i].i[1],nalist[k][j][i].i[0],&alist[k-1][j  ][i  ].i[1],alist[k-1][j  ][i  ].i[0],&superset[1]);
                    memcpy(nalist[k][j][i].i,superset,2*MAXAP);          
                    nalist[k][j][i].i[0] = difference_uint16(&superset[1],superset[0],&alist[k][j][i].i[1],alist[k][j][i].i[0],&nalist[k][j][i].i[1]);
                }
            }
        }       
        ierr = DMDAVecRestoreArray(user->daIdx,user->activeneighbourphases,&nalist); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(user->daIdx,user->activephases,&alist); CHKERRQ(ierr);
        /* forward solution... */
        user->step++; user->time += user->dt;
        user->error1 = user->error0; user->error0 = user->error; user->error = error;
        ierr = VecCopy(globalX,X); CHKERRQ(ierr);
        ierr = VecCopy(globalComp,user->cvec); CHKERRQ(ierr);
        ierr = VecCopy(globalXdot,user->rhs0 ); CHKERRQ(ierr);
        /* adapt time step... */
        if (error > 1.0) {
            user->dt /= 2.0;
        } else if (user->error1 && user->error0 && user->error) {
            user->dt *= pow(user->error ,-user->kI/4.0)
                     *  pow(user->error0,-user->kI/2.0)
                     *  pow(user->error1,-user->kI/4.0);
        }
        user->dt = user->dt < user->dtmax ? user->dt : user->dtmax;
    }
    ierr = DMRestoreGlobalVector(user->daIdx,&gactivephases); CHKERRQ(ierr);  
    ierr = DMRestoreGlobalVector(user->daCmp,&globalComp); CHKERRQ(ierr);
    ierr = DMRestoreGlobalVector(da,&globalX   ); CHKERRQ(ierr);
    ierr = DMRestoreGlobalVector(da,&globalXdot); CHKERRQ(ierr);
    ierr = DMRestoreGlobalVector(da,&globalXe  ); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

/*
 PostStep - Update active phases and write output
 */
PetscErrorCode PostStep(DM da,Vec X,AppCtx *user)
{
    PetscErrorCode    ierr;
    F_DOFS            ***gfdof;
    C_DOFS            ***cdof;
    F2I               ***alist;
    PetscInt          xs, ys, zs, xm, ym, zm;    
  
    PetscFunctionBeginUser;    
    /* write output */
    if (user->step%user->outputfreq == 0) {
       Vec               Xout[user->nc+1];
       O_DOFS            xout[user->nc+1];
       char              name[256];
       PetscViewer       viewer;
       
       for (PetscInt c=0; c<user->nc+1; c++) {
           ierr = DMGetGlobalVector(user->daOut,&Xout[c]); CHKERRQ(ierr);
           ierr = DMDAVecGetArray(user->daOut,Xout[c],&xout[c].o); CHKERRQ(ierr);
       }
       ierr = DMDAVecGetArray(da,X,&gfdof); CHKERRQ(ierr);
       ierr = DMDAVecGetArray(user->daCmp,user->cvec,&cdof); CHKERRQ(ierr);
       ierr = DMDAVecGetArray(user->daIdx,user->activephases,&alist); CHKERRQ(ierr);
       ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
       zm+=zs; ym+=ys; xm+=xs; 
       for (PetscInt k=zs; k<zm; k++) {
           for (PetscInt j=ys; j<ym; j++) {
               for (PetscInt i=xs; i<xm; i++) {
                   PetscScalar max = -LARGE, interpolant[alist[k][j][i].i[0]];
                   EvalInterpolant(interpolant,gfdof[k][j][i].p,alist[k][j][i].i[0]);
                   for (PetscInt g=0; g<alist[k][j][i].i[0]; g++) {
                       if (interpolant[g] > max) {
                          max = interpolant[g];
                          xout[0].o[k][j][i] = (PetscScalar) alist[k][j][i].i[g+1];
                       }
                   }
                   for (PetscInt c=0; c<user->nc; c++) {
                       xout[c+1].o[k][j][i] = 0.0;
                       for (PetscInt g=0; g<alist[k][j][i].i[0]; g++)
                           xout[c+1].o[k][j][i] += interpolant[g]*cdof[k][j][i].c[g*user->nc+c];        
                   }
               }
           }
       }       
       ierr = DMDAVecRestoreArray(user->daIdx,user->activephases,&alist); CHKERRQ(ierr);
       ierr = DMDAVecRestoreArray(user->daCmp,user->cvec,&cdof); CHKERRQ(ierr);
       ierr = DMDAVecRestoreArray(da,X,&gfdof); CHKERRQ(ierr);
       sprintf(name, "%s_%d.vtr",user->outfile,user->step);
       PetscPrintf(PETSC_COMM_WORLD,"writing output to %s\n",name);
       ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,name,FILE_MODE_WRITE,&viewer);
       ierr = DMDAVecRestoreArray(user->daOut,Xout[0],&xout[0].o); CHKERRQ(ierr);
       ierr = PetscObjectSetName((PetscObject) Xout[0],"phase");
       ierr = VecView(Xout[0],viewer);
       ierr = DMRestoreGlobalVector(user->daOut,&Xout[0]); CHKERRQ(ierr);
       for (PetscInt c=0; c<user->nc; c++) {
           ierr = DMDAVecRestoreArray(user->daOut,Xout[c+1],&xout[c+1].o); CHKERRQ(ierr);
           sprintf(name,"composition %s",user->componentname[c]);
           ierr = PetscObjectSetName((PetscObject) Xout[c+1],name);
           ierr = VecView(Xout[c+1],viewer);
           ierr = DMRestoreGlobalVector(user->daOut,&Xout[c+1]); CHKERRQ(ierr);
       }
       ierr = PetscViewerDestroy(&viewer);
    }
    PetscFunctionReturn(0);
}

int main(int argc,char **args)
{
  /* user-defined work context */
  AppCtx         *user = (AppCtx *) malloc(sizeof(struct AppCtx));         
  /* solution grid */
  DM             da;
  /* approx solution */
  Vec            x;
  /* numerical parameters */
  PetscReal      finaltime;
  PetscErrorCode ierr;
      
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Initialize program
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  assert(MAXAP == MAXTP*sizeof(PetscScalar)/sizeof(uint16_t));
  ierr = PetscInitialize(&argc,&args,(char*)0,help); CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Initialize problem parameters
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SetUpGeometry(user); CHKERRQ(ierr);
  ierr = SetUpInterface(user); CHKERRQ(ierr);
  utility_init(user);
  material_init(user);
  finaltime = 10.0;
  ierr = PetscOptionsGetReal(NULL,NULL,"-finaltime",&finaltime,NULL); CHKERRQ(ierr);
  user->dt = 0.5;
  ierr = PetscOptionsGetReal(NULL,NULL,"-mintimestep",&user->dt,NULL); CHKERRQ(ierr);
  user->dtmax = 10.0;
  ierr = PetscOptionsGetReal(NULL,NULL,"-maxtimestep",&user->dtmax,NULL); CHKERRQ(ierr);
  user->len = 4.0;
  ierr = PetscOptionsGetReal(NULL,NULL,"-l",&user->len,NULL); CHKERRQ(ierr);
  user->ptol = TOL;
  ierr = PetscOptionsGetReal(NULL,NULL,"-ptol",&user->ptol,NULL); CHKERRQ(ierr);
  user->ctol = TOL;
  ierr = PetscOptionsGetReal(NULL,NULL,"-ctol",&user->ctol,NULL); CHKERRQ(ierr);
  user->kI = 0.4;
  ierr = PetscOptionsGetReal(NULL,NULL,"-kI",&user->kI,NULL); CHKERRQ(ierr);
  user->outputfreq = 1;
  ierr = PetscOptionsGetInt(NULL,NULL,"-outfreq",&user->outputfreq,NULL); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL,NULL,"-outfile",user->outfile,128,NULL);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create distributed array (DMDA) to manage parallel grid and vectors
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMDACreate3d(PETSC_COMM_WORLD, 
                      DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,
                      DMDA_STENCIL_BOX, 
                      user->resolution[0],user->resolution[1],user->resolution[2],
                      PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
                      MAXAP+MAXCP,1,
                      NULL,NULL,NULL,
                      &da); CHKERRQ(ierr);
  ierr = DMSetFromOptions(da); CHKERRQ(ierr);
  ierr = DMSetUp(da); CHKERRQ(ierr);
  ierr = DMSetApplicationContext(da,user); CHKERRQ(ierr);
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create aux distributed array (DMDA) to manage output
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMDACreate3d(PETSC_COMM_WORLD, 
                      DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,
                      DMDA_STENCIL_BOX, 
                      user->resolution[0],user->resolution[1],user->resolution[2],
                      PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
                      1,1,
                      NULL,NULL,NULL,
                      &user->daOut); CHKERRQ(ierr);
  ierr = DMSetUp(user->daOut); CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create aux distributed array (DMDA) to manage bitwise data
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMDACreate3d(PETSC_COMM_WORLD, 
                      DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,
                      DMDA_STENCIL_BOX, 
                      user->resolution[0],user->resolution[1],user->resolution[2],
                      PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
                      MAXTP,1,
                      NULL,NULL,NULL,
                      &user->daIdx); CHKERRQ(ierr);
  ierr = DMSetUp(user->daIdx); CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create aux distributed array (DMDA) to manage bitwise data
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMDACreate3d(PETSC_COMM_WORLD, 
                      DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,
                      DMDA_STENCIL_BOX, 
                      user->resolution[0],user->resolution[1],user->resolution[2],
                      PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
                      MAXCP*MAXAP,1,
                      NULL,NULL,NULL,
                      &user->daCmp); CHKERRQ(ierr);
  ierr = DMSetUp(user->daCmp); CHKERRQ(ierr);

  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Extract global vectors from DMDA; then duplicate for remaining
   vectors that are the same types
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMCreateGlobalVector(da,&x); CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(da,&user->rhs0 ); CHKERRQ(ierr);
  ierr = VecZeroEntries(user->rhs0 ); CHKERRQ(ierr);
  ierr = DMCreateLocalVector(user->daIdx,&user->activeneighbourphases); CHKERRQ(ierr);
  ierr = DMCreateLocalVector(user->daIdx,&user->activephases); CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(user->daCmp,&user->cvec); CHKERRQ(ierr);
  ierr = SetUpProblem(da,user,x); CHKERRQ(ierr);
  
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Perform time integration 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  user->time = 0.0; user->error = 0.0; user->step = 0;
  for (;user->time < finaltime;) {
      ierr = PostStep(da,x,user);
      ierr = EulerStep(da,x,user);
  }

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Free work space.  All PETSc objects should be destroyed when they
   are no longer needed.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = VecDestroy(&x); CHKERRQ(ierr);
  ierr = VecDestroy(&user->rhs0); CHKERRQ(ierr);
  ierr = DMDestroy(&da); CHKERRQ(ierr);
  ierr = VecDestroy(&user->activeneighbourphases); CHKERRQ(ierr);
  ierr = VecDestroy(&user->activephases); CHKERRQ(ierr);
  ierr = DMDestroy(&user->daIdx); CHKERRQ(ierr);
  ierr = VecDestroy(&user->cvec); CHKERRQ(ierr);
  ierr = DMDestroy(&user->daCmp); CHKERRQ(ierr);
  ierr = DMDestroy(&user->daOut); CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}
