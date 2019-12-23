/**
 * Initialisation module file
 */
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <petscdm.h>
#include <petscdmplex.h>
#include <petscdmforest.h>
#include <petscviewerhdf5.h>
#include "material.h"
#include "utility.h"
#include "typedef.h"

/* Including my own header for checking by compiler */
#define INIT_IMPORT
#include "init.h"

/*
 SetUpGeometry - Import initial microstructure from geom file
 */
PetscErrorCode SetUpGeometry(AppCtx *user)
{
    char           geomfile[PETSC_MAX_PATH_LEN] = "";
    char           *buffer = 0, *substr = 0;
    MATERIAL       *currentmaterial;
    PetscMPIInt    rank;
    
    PetscFunctionBeginUser;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    PetscOptionsGetString(NULL,NULL,"-geomfile",geomfile,PETSC_MAX_PATH_LEN,NULL);
    FILE *infile=NULL;
    if (rank == 0) {
        infile = fopen (geomfile, "r");
        if (infile==NULL) {
            printf("Error: input file can't be opened. \n");
            printf("Error: Please make sure the file %s is in the CWD.\n",geomfile);
            return 1;
        }
    }
    PetscInt filesize;
    if (infile) {
        fseek(infile, 0, SEEK_END);
        filesize = ftell(infile);
        fseek(infile, 0, SEEK_SET);
    }    
    MPI_Bcast(&filesize,1,MPIU_INT ,0,PETSC_COMM_WORLD);
    buffer = malloc(filesize);
    if (infile) {
        fread(buffer, 1, filesize, infile);
        fclose(infile);
    }
    MPI_Bcast(buffer,filesize,MPI_CHAR,0,PETSC_COMM_WORLD);
    
    char *tok, *savetok;
    substr = Extract(buffer, "<header>", "</header>");
    tok = strtok_r(substr, "\n", &savetok);
    /* initialise header information */
    user->dim = 3;  //why doesn't this need to be &user->dim ?
    user->ncp = 1;
    user->npf = 1;
    user->nmat = 1;
    user->interpolation = LINEAR_INTERPOLATION;
    user->resolution = malloc(user->dim*sizeof(PetscInt));
    user->size = malloc(user->dim*sizeof(PetscReal));
    for (PetscInt dim=0; dim<user->dim; ++dim) {
        user->resolution[dim] = 1;
        user->size[dim] = 1.0;
    }
    while (tok !=NULL) {
        // process the line
        if (strstr(tok, "dimension ") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%d", &user->dim);
            free(user->resolution); user->resolution = malloc(user->dim*sizeof(PetscInt));
            free(user->size); user->size = malloc(user->dim*sizeof(PetscReal));
            for (PetscInt dim=0; dim<user->dim; ++dim) {
                user->resolution[dim] = 1;
                user->size[dim] = 1.0;
            }
        }    
        if (strstr(tok, "grid ") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            for (PetscInt dim=0; dim<user->dim; ++dim) {
                mtok = strtok_r(NULL, " ", &savemtok);
                user->resolution[dim] = atoi(strtok_r(NULL, " ", &savemtok));
            }
            user->phasevoxelmapping = malloc(user->resolution[2]*user->resolution[1]*user->resolution[0]*sizeof(uint16_t));
        }
        if (strstr(tok, "size ") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            for (PetscInt dim=0; dim<user->dim; ++dim) {
                mtok = strtok_r(NULL, " ", &savemtok);
                user->size[dim] = atof(strtok_r(NULL, " ", &savemtok));
            }
        }
        if (strstr(tok, "n_phases ") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%d", &user->npf);
            user->phasematerialmapping = malloc(user->npf*sizeof(uint16_t));
        }    
        if (strstr(tok, "n_materials ") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%d", &user->nmat);
            user->material = (MATERIAL *) malloc(user->nmat*sizeof(struct MATERIAL));
        }
        if (strstr(tok, "n_components ") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%d", &user->ncp);
            user->componentname = (char **) malloc(user->ncp*sizeof(char *));
            for (PetscInt c=0; c<user->ncp; c++) {
                user->componentname[c] = (char *) malloc(2*sizeof(char));
                sprintf(user->componentname[c],"%d",c+1); 
            }
            user->ndp = user->ncp-1;
        }    
        if (strstr(tok, "componentnames ") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            for (PetscInt c=0; c<user->ncp; c++) {
                sscanf(strtok_r(NULL, " ", &savemtok), "%s", user->componentname[c]);
            }
        }
        /* interpolation type */
        if (strstr(tok, "interpolation_type ") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            mtok = strtok_r(NULL, " ", &savemtok);
            if        (strstr(mtok, "linear"   ) != NULL) {
                user->interpolation = LINEAR_INTERPOLATION;
            } else if (strstr(mtok, "quadratic") != NULL) {
                user->interpolation = QUADRATIC_INTERPOLATION;
            } else if (strstr(mtok, "cubic"    ) != NULL) {
                user->interpolation = CUBIC_INTERPOLATION;
            } else {
                user->interpolation = NONE_INTERPOLATION;
            }
        }
        //advance the token
        tok = strtok_r(NULL, "\n", &savetok);
    }
    /* header information sanity checks */
    for (PetscInt dim=0; dim<user->dim; ++dim) {
        assert(user->resolution[dim] > 0);
        assert(user->size[dim] > 0);
    }
    assert(user->ncp  >  0    );
    assert(user->ncp  <= MAXCP);
    assert(user->npf  >  0    );
    assert(user->nmat >  0    );
    assert(user->interpolation != NONE_INTERPOLATION);
    free(substr);
    
    /* field offsets */
    AS_OFFSET = 0;
    AS_SIZE   = (MAXAP < user->npf ? MAXAP : user->npf) + 1;
    PF_OFFSET = AS_OFFSET+AS_SIZE;
    PF_SIZE   = (MAXAP < user->npf ? MAXAP : user->npf);
    DP_OFFSET = PF_OFFSET+PF_SIZE;
    DP_SIZE   = PF_SIZE*user->ndp;
    
    CP_OFFSET = 0;
    CP_SIZE   = PF_SIZE*user->ncp;
    /* initialise material information */
    currentmaterial = &user->material[0];
    for (PetscInt m=0; m<user->nmat; m++,currentmaterial++) {
        currentmaterial->model = NONE_CHEMENERGY;
        currentmaterial->c0 = malloc(user->ncp*sizeof(PetscReal));
        memset(currentmaterial->c0,0,user->ncp*sizeof(PetscReal));
        currentmaterial->mobilityc = (MOBILITY *) malloc(user->ncp*sizeof(MOBILITY));
        MOBILITY *currentmobility = &currentmaterial->mobilityc[0];
        for (PetscInt c=0; c<user->ncp; c++, currentmobility++) {
            currentmobility->m0 = 0.0;
            currentmobility->unary = (TSeries *) malloc(user->ncp*sizeof(TSeries));
            currentmobility->binary = (RK *) malloc(user->ncp*(user->ncp-1)/2*sizeof(struct RK));
            TSeries *currentunary = &currentmobility->unary[0];
            RK *currentbinary = &currentmobility->binary[0];
            for (PetscInt ck=0; ck<user->ncp; ck++,currentunary++) {
                currentunary->nTser = 0;
                currentunary->logCoeff = 0.0;
                for (PetscInt cj=ck+1; cj<user->ncp; cj++,currentbinary++) {
                    currentbinary->n = 0;
                }
            } 
        }
    }
    currentmaterial = &user->material[0];
    for (PetscInt m=0; m<user->nmat; m++,currentmaterial++) {
        QUAD *currentquad = &currentmaterial->energy.quad;
        CALPHAD *currentcalphad = &currentmaterial->energy.calphad;
        char starttag[128], endtag[128];
        sprintf(starttag, "<material %d>",m+1);        
        sprintf(endtag  , "</material %d>",m+1);        
        substr = Extract(buffer,starttag,endtag);
        assert(substr != NULL);
        tok = strtok_r(substr, "\n", &savetok);
        while (tok !=NULL) {
            /* model type */
            if (strstr(tok, "chemicalenergy_type ") != NULL) {
                char *mtok, *savemtok;
                mtok = strtok_r(tok, " ", &savemtok);
                mtok = strtok_r(NULL, " ", &savemtok);
                if (strstr(mtok, "quadratic") != NULL) {
                    currentmaterial->model = QUADRATIC_CHEMENERGY;
                    currentquad->ceq = (TSeries *) malloc(user->ncp*sizeof(TSeries));
                    currentquad->ref.nTser = 0;
                    currentquad->ref.logCoeff = 0.0;
                    currentquad->unary = (TSeries *) malloc(user->ncp*sizeof(TSeries));
                    currentquad->binary = (TSeries *) malloc(user->ncp*sizeof(TSeries));
                    currentquad->mobilityc = malloc(user->ncp*sizeof(PetscReal));
                    memset(currentquad->mobilityc,0,user->ncp*sizeof(PetscReal));
                    for (PetscInt c=0; c<user->ncp; c++) {
                        currentquad->ceq[c].nTser = 0;
                        currentquad->unary[c].nTser = 0;
                        currentquad->binary[c].nTser = 0;
                        currentquad->ceq[c].logCoeff = 0.0;
                        currentquad->unary[c].logCoeff = 0.0;
                        currentquad->binary[c].logCoeff = 0.0;
                    }
                } else if (strstr(mtok, "calphaddis") != NULL) {
                    currentmaterial->model = CALPHAD_CHEMENERGY;
                    currentcalphad->ref.nTser = 0;
                    currentcalphad->ref.logCoeff = 0.0;
                    currentcalphad->unary = (TSeries *) malloc(user->ncp*sizeof(TSeries));
                    currentcalphad->binary = (RK *) malloc(user->ncp*(user->ndp)*sizeof(struct RK)/2);
                    currentcalphad->ternary = (RK *) malloc(user->ncp*(user->ndp)*(user->ncp-2)*sizeof(struct RK)/6);
                    currentcalphad->mobilityc = malloc(user->ncp*sizeof(PetscReal));
                    memset(currentcalphad->mobilityc,0,user->ncp*sizeof(PetscReal));
                    TSeries *currentunary = &currentcalphad->unary[0];
                    RK *currentbinary = &currentcalphad->binary[0];
                    RK *currentternary = &currentcalphad->ternary[0];
                    for (PetscInt ck=0; ck<user->ncp; ck++,currentunary++) {
                        currentunary->nTser = 0;
                        currentunary->logCoeff = 0.0;
                        for (PetscInt cj=ck+1; cj<user->ncp; cj++,currentbinary++) {
                            currentbinary->n = 0;
                            for (PetscInt ci=cj+1; ci<user->ncp; ci++,currentternary++) { 
                                currentternary->n = 0;
                            }    
                        }
                    } 
                } else {
                    currentmaterial->model = NONE_CHEMENERGY;
                    currentmaterial->c0[0] = 1.0;
                }
            }
            /* molar volume */
            if (strstr(tok, "molarvolume ") != NULL) {
                char *mtok, *savemtok;
                mtok = strtok_r(tok, " ", &savemtok);
                sscanf(strtok_r(NULL, " ", &savemtok), "%lf", &currentmaterial->molarvolume);
            }
            /* component mobility for each material */
            if (strstr(tok, "mobilityc_") != NULL) {
                char *mtok, *savemtok;
                for (PetscInt c=0; c<user->ncp; c++) {
                    sprintf(starttag, "mobilityc_%s ",user->componentname[c]);
                    if (strstr(tok, starttag) != NULL){
                        mtok = strtok_r(tok, " ", &savemtok);
                        if (currentmaterial->model == QUADRATIC_CHEMENERGY) {
                            sscanf(strtok_r(NULL, " ", &savemtok), "%lf", &currentquad->mobilityc[c]);
                        } else if (currentmaterial->model == CALPHAD_CHEMENERGY) {
                            sscanf(strtok_r(NULL, " ", &savemtok), "%lf", &currentcalphad->mobilityc[c]);
                        }
                        c = user->ncp;
                    }
                }
            }
            /* component mobility for each material */
            if (strstr(tok, "mobilityc0_") != NULL) {
                char *mtok, *savemtok;
                MOBILITY *currentmobility = &currentmaterial->mobilityc[0];
                for (PetscInt c=0; c<user->ncp; c++, currentmobility++) {
                    sprintf(starttag, "mobilityc0_%s ",user->componentname[c]);
                    if (strstr(tok, starttag) != NULL){
                        mtok = strtok_r(tok, " ", &savemtok);
                        sscanf(strtok_r(NULL, " ", &savemtok), "%lf", &currentmobility->m0);
                        c = user->ncp;
                    }
                }
            }
            /* component migration energy for each material */
            if (strstr(tok, "migration_unary_coeff_") != NULL) {
                char *mtok, *savemtok;
                MOBILITY *currentmobility = &currentmaterial->mobilityc[0];
                for (PetscInt ck=0; ck<user->ncp; ck++, currentmobility++) {
                    TSeries *currentmigration = &currentmobility->unary[0];
                    for (PetscInt cj=0; cj<user->ncp; cj++, currentmigration++) {
                        sprintf(starttag, "migration_unary_coeff_%s_%s ",user->componentname[ck],user->componentname[cj]);
                        if (strstr(tok, starttag) != NULL){
                            mtok = strtok_r(tok, " ", &savemtok);
                            mtok = strtok_r(NULL, " ", &savemtok);
                            currentmigration->nTser = 0;
                            while (mtok != NULL) {
                                sscanf(mtok, "%lf", &currentmigration->coeff[currentmigration->nTser]);
                                currentmigration->nTser++;
                                mtok = strtok_r(NULL, " ", &savemtok);
                            }
                            ck=user->ncp;cj=user->ncp;
                        }
                    }
                }
            }
            if (strstr(tok, "migration_unary_exp_") != NULL) {
                char *mtok, *savemtok;
                MOBILITY *currentmobility = &currentmaterial->mobilityc[0];
                for (PetscInt ck=0; ck<user->ncp; ck++, currentmobility++) {
                    TSeries *currentmigration = &currentmobility->unary[0];
                    for (PetscInt cj=0; cj<user->ncp; cj++, currentmigration++) {
                        sprintf(starttag, "migration_unary_exp_%s_%s ",user->componentname[ck],user->componentname[cj]);
                        if (strstr(tok, starttag) != NULL){
                            mtok = strtok_r(tok, " ", &savemtok);
                            mtok = strtok_r(NULL, " ", &savemtok);
                            currentmigration->nTser = 0;
                            while (mtok != NULL) {
                                sscanf(mtok, "%d", &currentmigration->exp[currentmigration->nTser]);
                                currentmigration->nTser++;
                                mtok = strtok_r(NULL, " ", &savemtok);
                            }
                            ck=user->ncp;cj=user->ncp;
                        }
                    }
                }
            }
            if (strstr(tok, "migration_nbinary_") != NULL) {
                char *mtok, *savemtok;
                MOBILITY *currentmobility = &currentmaterial->mobilityc[0];
                for (PetscInt ck=0; ck<user->ncp; ck++, currentmobility++) {
                    RK *currentbinary = &currentmobility->binary[0];
                    for (PetscInt cj=0; cj<user->ncp; cj++) {
                        for (PetscInt ci=cj+1; ci<user->ncp; ci++, currentbinary++) {
                            sprintf(starttag, "migration_binary_%s_%s_%s ",user->componentname[ck],user->componentname[cj],user->componentname[ci]);
                            if (strstr(tok, starttag) != NULL){
                                mtok = strtok_r(tok, " ", &savemtok);
                                mtok = strtok_r(NULL, " ", &savemtok);
                                sscanf(mtok, "%d", &currentbinary->n);
                                currentbinary->enthalpy = (TSeries *) malloc(currentbinary->n*sizeof(TSeries));
                                TSeries *currentmigration = &currentbinary->enthalpy[0];
                                for (PetscInt nrk=0; nrk < currentbinary->n; nrk++, currentmigration++) {
                                    currentmigration->nTser = 0;
                                    currentmigration->logCoeff = 0.0;
                                }
                                ci=user->ncp;cj=user->ncp;ck=user->ncp;
                            }    
                        }    
                    }
                }
            }
            if (strstr(tok, "migration_binary_coeff_") != NULL) {
                char *mtok, *savemtok;
                MOBILITY *currentmobility = &currentmaterial->mobilityc[0];
                for (PetscInt ck=0; ck<user->ncp; ck++, currentmobility++) {
                    RK *currentbinary = &currentmobility->binary[0];
                    for (PetscInt cj=0; cj<user->ncp; cj++) {
                        for (PetscInt ci=cj+1; ci<user->ncp; ci++, currentbinary++) {
                            TSeries *currentmigration = &currentbinary->enthalpy[0];
                            for (PetscInt nrk=0; nrk < currentbinary->n; nrk++, currentmigration++) {
                                sprintf(starttag, "migration_binary_coeff_%s_%s_%s_%d ",user->componentname[ck],user->componentname[cj],user->componentname[ci],nrk+1);
                                if (strstr(tok, starttag) != NULL){
                                    mtok = strtok_r(tok, " ", &savemtok);
                                    mtok = strtok_r(NULL, " ", &savemtok);
                                    currentmigration->nTser = 0;
                                    while (mtok != NULL) {
                                        sscanf(mtok, "%lf", &currentmigration->coeff[currentmigration->nTser]);
                                        currentmigration->nTser++;
                                        mtok = strtok_r(NULL, " ", &savemtok);
                                    }
                                    ck=user->ncp;cj=user->ncp;ci=user->ncp;nrk=currentbinary->n;
                                }
                            }
                        }    
                    }
                }
            }
            if (strstr(tok, "migration_binary_exp_") != NULL) {
                char *mtok, *savemtok;
                MOBILITY *currentmobility = &currentmaterial->mobilityc[0];
                for (PetscInt ck=0; ck<user->ncp; ck++, currentmobility++) {
                    RK *currentbinary = &currentmobility->binary[0];
                    for (PetscInt cj=0; cj<user->ncp; cj++) {
                        for (PetscInt ci=cj+1; ci<user->ncp; ci++, currentbinary++) {
                            TSeries *currentmigration = &currentbinary->enthalpy[0];
                            for (PetscInt nrk=0; nrk < currentbinary->n; nrk++, currentmigration++) {
                                sprintf(starttag, "migration_binary_exp_%s_%s_%s_%d ",user->componentname[ck],user->componentname[cj],user->componentname[ci],nrk+1);
                                if (strstr(tok, starttag) != NULL){
                                    mtok = strtok_r(tok, " ", &savemtok);
                                    mtok = strtok_r(NULL, " ", &savemtok);
                                    currentmigration->nTser = 0;
                                    while (mtok != NULL) {
                                        sscanf(mtok, "%d", &currentmigration->exp[currentmigration->nTser]);
                                        currentmigration->nTser++;
                                        mtok = strtok_r(NULL, " ", &savemtok);
                                    }
                                    ck=user->ncp;cj=user->ncp;ci=user->ncp;nrk=currentbinary->n;
                                }
                            }
                        }    
                    }
                }
            }
            /* initial composition */
            if (strstr(tok, "c0 ") != NULL) {
                char *mtok, *savemtok;
                mtok = strtok_r(tok, " ", &savemtok);
                mtok = strtok_r(NULL, " ", &savemtok);
                PetscInt c=0;
                while (mtok != NULL) {
                    PetscReal c0;
                    sscanf(mtok, "%lf", &c0);
                    currentmaterial->c0[c] = c0;
                    mtok = strtok_r(NULL, " ", &savemtok);
                    c++;
                }
                assert(c == user->ncp);
            }
            /* initial composition */
            if (strstr(tok, "statekineticcoeff ") != NULL) {
                char *mtok, *savemtok;
                mtok = strtok_r(tok, " ", &savemtok);
                sscanf(strtok_r(NULL, " ", &savemtok), "%lf", &currentmaterial->statekineticcoeff);
            }
            /* equillibrium composition */
            if (strstr(tok, "quad_ceq_coeff_") != NULL) {
                char *mtok, *savemtok;
                TSeries *currentceq = &currentquad->ceq[0];
                for (PetscInt c=0; c<user->ncp; c++, currentceq++) {
                    sprintf(starttag, "quad_ceq_coeff_%s ",user->componentname[c]);
                    if (strstr(tok, starttag) != NULL){
                        mtok = strtok_r(tok, " ", &savemtok);
                        mtok = strtok_r(NULL, " ", &savemtok);
                        currentceq->nTser = 0;
                        while (mtok != NULL) {
                            sscanf(mtok, "%lf", &currentceq->coeff[currentceq->nTser]);
                            currentceq->nTser++;
                            mtok = strtok_r(NULL, " ", &savemtok);
                        }
                        c=user->ncp;
                    }
                }
            }
            if (strstr(tok, "quad_ceq_exp_") != NULL) {
                char *mtok, *savemtok;
                TSeries *currentceq = &currentquad->ceq[0];
                for (PetscInt c=0; c<user->ncp; c++, currentceq++) {
                    sprintf(starttag, "quad_ceq_exp_%s ",user->componentname[c]);
                    if (strstr(tok, starttag) != NULL){
                        mtok = strtok_r(tok, " ", &savemtok);
                        mtok = strtok_r(NULL, " ", &savemtok);
                        currentceq->nTser = 0;
                        while (mtok != NULL) {
                            sscanf(mtok, "%d", &currentceq->exp[currentceq->nTser]);
                            currentceq->nTser++;
                            mtok = strtok_r(NULL, " ", &savemtok);
                        }
                        c=user->ncp;
                    }
                }
            }
            /* F_chem parameters for each material */
            if (strstr(tok, "quad_refenthalpy_coeff ") != NULL) {
                char *mtok, *savemtok;
                mtok = strtok_r(tok, " ", &savemtok);
                mtok = strtok_r(NULL, " ", &savemtok);
                currentquad->ref.nTser = 0;
                while (mtok != NULL) {
                    sscanf(mtok, "%lf", &currentquad->ref.coeff[currentquad->ref.nTser]);
                    currentquad->ref.nTser++;
                    mtok = strtok_r(NULL, " ", &savemtok);
                }
            }
            if (strstr(tok, "quad_refenthalpy_exp ") != NULL) {
                char *mtok, *savemtok;
                mtok = strtok_r(tok, " ", &savemtok);
                mtok = strtok_r(NULL, " ", &savemtok);
                currentquad->ref.nTser = 0;
                while (mtok != NULL) {
                    sscanf(mtok, "%d", &currentquad->ref.exp[currentquad->ref.nTser]);
                    currentquad->ref.nTser++;
                    mtok = strtok_r(NULL, " ", &savemtok);
                }
            }
            /* F_chem parameters for each material */
            if (strstr(tok, "quad_unaryenthalpy_coeff_") != NULL) {
                char *mtok, *savemtok;
                TSeries *currentunary = &currentquad->unary[0];
                for (PetscInt c=0; c<user->ncp; c++, currentunary++) {
                    sprintf(starttag, "quad_unaryenthalpy_coeff_%s ",user->componentname[c]);
                    if (strstr(tok, starttag) != NULL){
                        mtok = strtok_r(tok, " ", &savemtok);
                        mtok = strtok_r(NULL, " ", &savemtok);
                        currentunary->nTser = 0;
                        while (mtok != NULL) {
                            sscanf(mtok, "%lf", &currentunary->coeff[currentunary->nTser]);
                            currentunary->nTser++;
                            mtok = strtok_r(NULL, " ", &savemtok);
                        }
                        c=user->ncp;
                    }
                }
            }
            if (strstr(tok, "quad_unaryenthalpy_exp_") != NULL) {
                char *mtok, *savemtok;
                TSeries *currentunary = &currentquad->unary[0];
                for (PetscInt c=0; c<user->ncp; c++, currentunary++) {
                    sprintf(starttag, "quad_unaryenthalpy_exp_%s ",user->componentname[c]);
                    if (strstr(tok, starttag) != NULL){
                        mtok = strtok_r(tok, " ", &savemtok);
                        mtok = strtok_r(NULL, " ", &savemtok);
                        currentunary->nTser = 0;
                        while (mtok != NULL) {
                            sscanf(mtok, "%d", &currentunary->exp[currentunary->nTser]);
                            currentunary->nTser++;
                            mtok = strtok_r(NULL, " ", &savemtok);
                        }
                        c=user->ncp;
                    }
                }
            }
            /* F_chem parameters for each material */
            if (strstr(tok, "quad_binaryenthalpy_coeff_") != NULL) {
                char *mtok, *savemtok;
                TSeries *currentbinary = &currentquad->binary[0];
                for (PetscInt c=0; c<user->ncp; c++, currentbinary++) {
                    sprintf(starttag, "quad_binaryenthalpy_coeff_%s ",user->componentname[c]);
                    if (strstr(tok, starttag) != NULL){
                        mtok = strtok_r(tok, " ", &savemtok);
                        mtok = strtok_r(NULL, " ", &savemtok);
                        currentbinary->nTser = 0;
                        while (mtok != NULL) {
                            sscanf(mtok, "%lf", &currentbinary->coeff[currentbinary->nTser]);
                            currentbinary->nTser++;
                            mtok = strtok_r(NULL, " ", &savemtok);
                        }
                        c=user->ncp;
                    }
                }
            }
            if (strstr(tok, "quad_binaryenthalpy_exp_") != NULL) {
                char *mtok, *savemtok;
                TSeries *currentbinary = &currentquad->binary[0];
                for (PetscInt c=0; c<user->ncp; c++, currentbinary++) {
                    sprintf(starttag, "quad_binaryenthalpy_exp_%s ",user->componentname[c]);
                    if (strstr(tok, starttag) != NULL){
                        mtok = strtok_r(tok, " ", &savemtok);
                        mtok = strtok_r(NULL, " ", &savemtok);
                        currentbinary->nTser = 0;
                        while (mtok != NULL) {
                            sscanf(mtok, "%d", &currentbinary->exp[currentbinary->nTser]);
                            currentbinary->nTser++;
                            mtok = strtok_r(NULL, " ", &savemtok);
                        }
                        c=user->ncp;
                    }
                }
            }
            /* F_chem parameters for each material */
            if (strstr(tok, "calphad_refenthalpy_coeff ") != NULL) {
                char *mtok, *savemtok;
                mtok = strtok_r(tok, " ", &savemtok);
                mtok = strtok_r(NULL, " ", &savemtok);
                currentcalphad->ref.nTser = 0;
                while (mtok != NULL) {
                    sscanf(mtok, "%lf", &currentcalphad->ref.coeff[currentcalphad->ref.nTser]);
                    currentcalphad->ref.nTser++;
                    mtok = strtok_r(NULL, " ", &savemtok);
                }
            }
            if (strstr(tok, "calphad_refenthalpy_exp ") != NULL) {
                char *mtok, *savemtok;
                mtok = strtok_r(tok, " ", &savemtok);
                mtok = strtok_r(NULL, " ", &savemtok);
                currentcalphad->ref.nTser = 0;
                while (mtok != NULL) {
                    sscanf(mtok, "%d", &currentcalphad->ref.exp[currentcalphad->ref.nTser]);
                    currentcalphad->ref.nTser++;
                    mtok = strtok_r(NULL, " ", &savemtok);
                }
            }
            /* F_chem parameters for each material */
            if (strstr(tok, "calphad_unaryenthalpy_coeff_") != NULL) {
                char *mtok, *savemtok;
                TSeries *currentunary = &currentcalphad->unary[0];
                for (PetscInt c=0; c<user->ncp; c++, currentunary++) {
                    sprintf(starttag, "calphad_unaryenthalpy_coeff_%s ",user->componentname[c]);
                    if (strstr(tok, starttag) != NULL){
                        mtok = strtok_r(tok, " ", &savemtok);
                        mtok = strtok_r(NULL, " ", &savemtok);
                        currentunary->nTser = 0;
                        while (mtok != NULL) {
                            sscanf(mtok, "%lf", &currentunary->coeff[currentunary->nTser]);
                            currentunary->nTser++;
                            mtok = strtok_r(NULL, " ", &savemtok);
                        }
                        c=user->ncp;
                    }
                }
            }
            if (strstr(tok, "calphad_unaryenthalpy_exp_") != NULL) {
                char *mtok, *savemtok;
                TSeries *currentunary = &currentcalphad->unary[0];
                for (PetscInt c=0; c<user->ncp; c++, currentunary++) {
                    sprintf(starttag, "calphad_unaryenthalpy_exp_%s ",user->componentname[c]);
                    if (strstr(tok, starttag) != NULL){
                        mtok = strtok_r(tok, " ", &savemtok);
                        mtok = strtok_r(NULL, " ", &savemtok);
                        currentunary->nTser = 0;
                        while (mtok != NULL) {
                            sscanf(mtok, "%d", &currentunary->exp[currentunary->nTser]);
                            currentunary->nTser++;
                            mtok = strtok_r(NULL, " ", &savemtok);
                        }
                        c=user->ncp;
                    }
                }
            }
            if (strstr(tok, "calphad_unaryenthalpy_logCoeff_") != NULL) { //if this works put in for binary and ternary also
                char *mtok, *savemtok;
                TSeries *currentunary = &currentcalphad->unary[0];
                for (PetscInt c=0; c<user->ncp; c++, currentunary++) {
                    sprintf(starttag, "calphad_unaryenthalpy_logCoeff_%s ",user->componentname[c]);
                    if (strstr(tok, starttag) != NULL){
                        mtok = strtok_r(tok, " ", &savemtok);
                        mtok = strtok_r(NULL, " ", &savemtok);
                        sscanf(mtok, "%lf", &currentunary->logCoeff);
                        c=user->ncp;
                    }
                }
            }
            if (strstr(tok, "calphad_nbinaryenthalpy_") != NULL) {
                char *mtok, *savemtok;
                RK *currentbinary = &currentcalphad->binary[0];
                for (PetscInt cj=0;cj<user->ncp;cj++) {
                    for (PetscInt ci=cj+1;ci<user->ncp;ci++, currentbinary++) {
                        sprintf(starttag, "calphad_nbinaryenthalpy_%s_%s ",user->componentname[cj],user->componentname[ci]);
                        if (strstr(tok, starttag) != NULL){
                            mtok = strtok_r(tok, " ", &savemtok);
                            mtok = strtok_r(NULL, " ", &savemtok);
                            sscanf(mtok, "%d", &currentbinary->n);
                            currentbinary->enthalpy = (TSeries *) malloc(currentbinary->n*sizeof(TSeries));
                            TSeries *currententhalpy = &currentbinary->enthalpy[0];
                            for (PetscInt nrk=0; nrk < currentbinary->n; nrk++, currententhalpy++) {
                                currententhalpy->nTser = 0;
                                currententhalpy->logCoeff = 0.0;
                            }
                            ci=user->ncp;cj=user->ncp;
                        }
                    }
                }
            }
            if (strstr(tok, "calphad_binaryenthalpy_coeff_") != NULL) {
                char *mtok, *savemtok;
                RK *currentbinary = &currentcalphad->binary[0];
                for (PetscInt cj=0;cj<user->ncp;cj++) {
                    for (PetscInt ci=cj+1;ci<user->ncp;ci++, currentbinary++) {
                        TSeries *currententhalpy = &currentbinary->enthalpy[0];
                        for (PetscInt nrk=0; nrk < currentbinary->n; nrk++, currententhalpy++) {
                            sprintf(starttag, "calphad_binaryenthalpy_coeff_%s_%s_%d ",user->componentname[cj],user->componentname[ci],nrk+1);
                            if (strstr(tok, starttag) != NULL){
                                mtok = strtok_r(tok, " ", &savemtok);
                                mtok = strtok_r(NULL, " ", &savemtok);
                                currententhalpy->nTser = 0;
                                while (mtok != NULL) {
                                    sscanf(mtok, "%lf", &currententhalpy->coeff[currententhalpy->nTser]);
                                    currententhalpy->nTser++;
                                    mtok = strtok_r(NULL, " ", &savemtok);
                                }
                                ci=user->ncp;cj=user->ncp;nrk=currentbinary->n;
                            }
                        }
                    }
                }
            }
            if (strstr(tok, "calphad_binaryenthalpy_exp_") != NULL) {
                char *mtok, *savemtok;
                RK *currentbinary = &currentcalphad->binary[0];
                for (PetscInt cj=0;cj<user->ncp;cj++) {
                    for (PetscInt ci=cj+1;ci<user->ncp;ci++, currentbinary++) {
                        TSeries *currententhalpy = &currentbinary->enthalpy[0];
                        for (PetscInt nrk=0; nrk < currentbinary->n; nrk++, currententhalpy++) {
                            sprintf(starttag, "calphad_binaryenthalpy_exp_%s_%s_%d ",user->componentname[cj],user->componentname[ci],nrk+1);
                            if (strstr(tok, starttag) != NULL){
                                mtok = strtok_r(tok, " ", &savemtok);
                                mtok = strtok_r(NULL, " ", &savemtok);
                                currententhalpy->nTser = 0;
                                while (mtok != NULL) {
                                    sscanf(mtok, "%d", &currententhalpy->exp[currententhalpy->nTser]);
                                    currententhalpy->nTser++;
                                    mtok = strtok_r(NULL, " ", &savemtok);
                                }
                                ci=user->ncp;cj=user->ncp;nrk=currentbinary->n;
                            }
                        }
                    }
                }
            }
            if (strstr(tok, "calphad_binaryenthalpy_logCoeff_") != NULL) { //if this works put in for binary and ternary also
                char *mtok, *savemtok;
                RK *currentbinary = &currentcalphad->binary[0];
                for (PetscInt cj=0;cj<user->ncp;cj++) {
                    for (PetscInt ci=cj+1;ci<user->ncp;ci++, currentbinary++) {
                        TSeries *currententhalpy = &currentbinary->enthalpy[0];
                        for (PetscInt nrk=0; nrk < currentbinary->n; nrk++, currententhalpy++) {
                            sprintf(starttag, "calphad_binaryenthalpy_logCoeff_%s_%s_%d ",user->componentname[cj],user->componentname[ci],nrk+1);
                            if (strstr(tok, starttag) != NULL){
                                mtok = strtok_r(tok, " ", &savemtok);
                                mtok = strtok_r(NULL, " ", &savemtok);
                                sscanf(mtok, "%lf", &currententhalpy->logCoeff);
                                ci=user->ncp;cj=user->ncp;nrk=currentbinary->n;
                            }
                        }
                    }
                }
            }                
            /* F_chem parameters for each material */
            if (strstr(tok, "calphad_nternaryenthalpy_") != NULL) {
                char *mtok, *savemtok;
                RK *currentternary = &currentcalphad->ternary[0];
                for (PetscInt ck=0;ck<user->ncp;ck++) {
                    for (PetscInt cj=ck+1;cj<user->ncp;cj++) {
                        for (PetscInt ci=cj+1;ci<user->ncp;ci++, currentternary++) {
                            sprintf(starttag, "calphad_nternaryenthalpy_%s_%s_%s ",user->componentname[ck],user->componentname[cj],user->componentname[ci]);
                            if (strstr(tok, starttag) != NULL){
                                mtok = strtok_r(tok, " ", &savemtok);
                                mtok = strtok_r(NULL, " ", &savemtok);
                                sscanf(mtok, "%d", &currentternary->n);
                                currentternary->i = malloc(currentternary->n*sizeof(PetscInt));
                                memset(currentternary->i,0,currentternary->n*sizeof(PetscInt));
                                currentternary->enthalpy = (TSeries *) malloc(currentternary->n*sizeof(TSeries));
                                TSeries *currententhalpy = &currentternary->enthalpy[0];
                                for (PetscInt nrk=0; nrk < currentternary->n; nrk++, currententhalpy++) {
                                    currententhalpy->nTser = 0;
                                    currententhalpy->logCoeff = 0.0;
                                }
                                ci=user->ncp;cj=user->ncp;ck=user->ncp;
                            }
                        }
                    }
                }
            }
            if (strstr(tok, "calphad_iternaryenthalpy_") != NULL) {
                char *mtok, *savemtok;
                RK *currentternary = &currentcalphad->ternary[0];
                for (PetscInt ck=0;ck<user->ncp;ck++) {
                    for (PetscInt cj=ck+1;cj<user->ncp;cj++) {
                        for (PetscInt ci=cj+1;ci<user->ncp;ci++, currentternary++) {
                            sprintf(starttag, "calphad_iternaryenthalpy_%s_%s_%s ",user->componentname[ck],user->componentname[cj],user->componentname[ci]);
                            if (strstr(tok, starttag) != NULL){
                                mtok = strtok_r(tok, " ", &savemtok);
                                mtok = strtok_r(NULL, " ", &savemtok);
                                PetscInt nrk = 0;
                                while (mtok != NULL) {
                                    PetscInt iternaryenthalpy;
                                    sscanf(mtok, "%d", &iternaryenthalpy);
                                    currentternary->i[nrk] = iternaryenthalpy-1;
                                    mtok = strtok_r(NULL, " ", &savemtok);
                                    nrk++;
                                }
                                ci=user->ncp;cj=user->ncp;ck=user->ncp;
                            }
                        }
                    }
                }
            }
            if (strstr(tok, "calphad_ternaryenthalpy_coeff_") != NULL) {
                char *mtok, *savemtok;
                RK *currentternary = &currentcalphad->ternary[0];
                for (PetscInt ck=0;ck<user->ncp;ck++) {
                    for (PetscInt cj=ck+1;cj<user->ncp;cj++) {
                        for (PetscInt ci=cj+1;ci<user->ncp;ci++, currentternary++) {
                            TSeries *currententhalpy = &currentternary->enthalpy[0];
                            for (PetscInt nrk=0; nrk < currentternary->n; nrk++, currententhalpy++) {
                                sprintf(starttag, "calphad_ternaryenthalpy_coeff_%s_%s_%s_%d ",user->componentname[ck],user->componentname[cj],user->componentname[ci],nrk+1);
                                if (strstr(tok, starttag) != NULL){
                                    mtok = strtok_r(tok, " ", &savemtok);
                                    mtok = strtok_r(NULL, " ", &savemtok);
                                    currententhalpy->nTser = 0;
                                    while (mtok != NULL) {
                                        sscanf(mtok, "%lf", &currententhalpy->coeff[currententhalpy->nTser]);
                                        currententhalpy->nTser++;
                                        mtok = strtok_r(NULL, " ", &savemtok);
                                    }
                                    ci=user->ncp;cj=user->ncp;ck=user->ncp;nrk=currentternary->n;
								}
							}
						}
					}
				}
			}
            if (strstr(tok, "calphad_ternaryenthalpy_exp_") != NULL) {
                char *mtok, *savemtok;
                RK *currentternary = &currentcalphad->ternary[0];
                for (PetscInt ck=0;ck<user->ncp;ck++) {
                    for (PetscInt cj=ck+1;cj<user->ncp;cj++) {
                        for (PetscInt ci=cj+1;ci<user->ncp;ci++, currentternary++) {
                            TSeries *currententhalpy = &currentternary->enthalpy[0];
                            for (PetscInt nrk=0; nrk < currentternary->n; nrk++, currententhalpy++) {
								sprintf(starttag, "calphad_ternaryenthalpy_exp_%s_%s_%s_%d ",user->componentname[ck],user->componentname[cj],user->componentname[ci],nrk+1);
								if (strstr(tok, starttag) != NULL){
									mtok = strtok_r(tok, " ", &savemtok);
									mtok = strtok_r(NULL, " ", &savemtok);
                                    currententhalpy->nTser = 0;
									while (mtok != NULL) {
                                        sscanf(mtok, "%d", &currententhalpy->exp[currententhalpy->nTser]);
                                        currententhalpy->nTser++;
										mtok = strtok_r(NULL, " ", &savemtok);
									}
									ci=user->ncp;cj=user->ncp;ck=user->ncp;nrk=currentternary->n;
	                            }
    	                    }
        	            }
            	    }
            	}
			}
            if (strstr(tok, "calphad_ternaryenthalpy_logCoeff_") != NULL) { //if this works put in for binary and ternary also
                char *mtok, *savemtok;
                RK *currentternary = &currentcalphad->ternary[0];
                for (PetscInt ck=0;ck<user->ncp;ck++) {
                    for (PetscInt cj=ck+1;cj<user->ncp;cj++) {
                        for (PetscInt ci=cj+1;ci<user->ncp;ci++, currentternary++) {
                            TSeries *currententhalpy = &currentternary->enthalpy[0];
                            for (PetscInt nrk=0; nrk < currentternary->n; nrk++, currententhalpy++) {
								sprintf(starttag, "calphad_ternaryenthalpy_logCoeff_%s_%s_%s_%d ",user->componentname[ck],user->componentname[cj],user->componentname[ci],nrk+1);
								if (strstr(tok, starttag) != NULL){
									mtok = strtok_r(tok, " ", &savemtok);
									mtok = strtok_r(NULL, " ", &savemtok);
                                    sscanf(mtok, "%lf", &currententhalpy->logCoeff);
									ci=user->ncp;cj=user->ncp;ck=user->ncp;nrk=currentternary->n;
                                }
                            }
                        }
                    }
                }
            }                
            //advance the token
            tok = strtok_r(NULL, "\n", &savetok);
        }
        /* material information sanity checks */
        assert(currentmaterial->molarvolume > 0.0);
        for (PetscInt c=0;c<user->ncp;c++) {
            assert(currentmaterial->c0[c] > 0.0);
        }
        if        (currentmaterial->model == QUADRATIC_CHEMENERGY) {
            QUAD *currentquad = &currentmaterial->energy.quad;
            for (PetscInt c=0;c<user->ncp;c++) {
                assert(currentquad->mobilityc[c] >= 0.0);
            }
        } else if (currentmaterial->model ==   CALPHAD_CHEMENERGY) {
            CALPHAD *currentcalphad = &currentmaterial->energy.calphad;
            for (PetscInt c=0;c<user->ncp;c++) {
                assert(currentcalphad->mobilityc[c] >= 0.0);
            }
        } else if (currentmaterial->model ==      NONE_CHEMENERGY) {
            assert(user->ncp == 1);
        }
        free(substr);
    }
    
    substr = Extract(buffer, "<phase_material_mapping>", "</phase_material_mapping>");
    tok = strtok_r(substr, "\n", &savetok);
    PetscInt ctrm=0;
    while (tok != NULL && ctrm < user->npf) {
        // process the line
        if (strstr(tok, "of") != NULL) {
            char stra[128], strb[128];
            PetscInt sof,m;
            sscanf(tok, "%s of %s", stra, strb);
            sof = atoi(stra); m = atoi(strb);
            assert(m <= user->nmat);
            for (PetscInt j=ctrm;j<ctrm+sof;j++) 
                user->phasematerialmapping[j] = m-1;
            ctrm += sof;
        } else if (strstr(tok, "to") != NULL) {
            char stra[128], strb[128];
            PetscInt m, mfrom, mto;
            sscanf(tok, "%s to %s", stra, strb);
            mfrom = atoi(stra); mto = atoi(strb);
            assert(mfrom <= user->nmat);
            assert(mto   <= user->nmat);
            if (mfrom < mto) {
                for (m=mfrom;m<=mto;m++) 
                    user->phasematerialmapping[ctrm++] = m-1;
            } else {
                for (m=mfrom;m>=mto;m--)
                    user->phasematerialmapping[ctrm++] = m-1;
            }   
        }
        else if (atoi(tok)) {
            assert(atoi(tok) <= user->nmat);
            user->phasematerialmapping[ctrm++] = atoi(tok)-1;
        }
        //advance the token
        tok = strtok_r(NULL, "\n", &savetok);
    }
    assert(ctrm == user->npf && tok == NULL);
    free(substr);    

    substr = Extract(buffer, "<voxel_phase_mapping>", "</voxel_phase_mapping>");
    tok = strtok_r(substr, "\n", &savetok);
    PetscInt ctrv=0;
    while (tok != NULL && ctrv < user->resolution[2]*user->resolution[1]*user->resolution[0]) {
        // process the line
        if (strstr(tok, "of") != NULL) {
            char stra[128], strb[128];
            PetscInt sof,p;
            sscanf(tok, "%s of %s", stra, strb);
            sof = atoi(stra); p = atoi(strb);
            for (PetscInt j=ctrv;j<ctrv+sof;j++) 
                user->phasevoxelmapping[j] = p-1;;
            ctrv+=sof;
        } else if (strstr(tok, "to") != NULL) {
            char stra[128], strb[128];
            PetscInt p, pfrom, pto;
            sscanf(tok, "%s to %s", stra, strb);
            pfrom = atoi(stra); pto = atoi(strb);
            if (pfrom < pto) {
                for (p=pfrom;p<=pto;p++) 
                    user->phasevoxelmapping[ctrv++] = p-1;
            } else {
                for (p=pfrom;p>=pto;p--)
                    user->phasevoxelmapping[ctrv++] = p-1;
            }   
        }
        else if (atoi(tok)) {
            user->phasevoxelmapping[ctrv++] = atoi(tok)-1;
        }
        //advance the token
        tok = strtok_r(NULL, "\n", &savetok);
    }
    assert(ctrv == user->resolution[0]*user->resolution[1]*user->resolution[2] && tok == NULL);
    free(substr);    

    substr = Extract(buffer, "<solution_parameters>", "</solution_parameters>");
    tok = strtok_r(substr, "\n", &savetok);
    /* default AMR parameters */
    user->amrparams.initrefine = 1;
    user->amrparams.initcoarsen = 0;
    user->amrparams.maxnrefine = 1;
    user->amrparams.minnrefine = 1;
    user->amrparams.initblocksize = malloc(user->dim*sizeof(PetscInt));
    for (PetscInt dim=0; dim<user->dim; ++dim) {
        user->amrparams.initblocksize[dim] = 2;
    }
    /* default solution parameters */
    user->solparams.finaltime = 1.0;
    user->solparams.timestep = TOL;
    user->solparams.mintimestep = TOL;
    user->solparams.maxtimestep = LARGE;
    user->solparams.step = 0;    
    user->solparams.interfacewidth = 4.0;
    user->solparams.temperature = 300.0;
    user->solparams.reltol = 1e-6;
    user->solparams.abstol = 1e-6;
    user->solparams.outputfreq = 1;
    strncpy(user->solparams.outfile, "output", 128);
    strncpy(user->solparams.petscoptions, "", PETSC_MAX_PATH_LEN);
    /* initialise solution parameters */
    while (tok !=NULL) {
        // process the line
        if (strstr(tok, "finaltime ") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%lf", &user->solparams.finaltime);
        }
        if (strstr(tok, "timestep0 ") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%lf", &user->solparams.timestep);
        }
        if (strstr(tok, "timestepmin ") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%lf", &user->solparams.mintimestep);
        }
        if (strstr(tok, "timestepmax ") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%lf", &user->solparams.maxtimestep);
        }
        if (strstr(tok, "interfacewidth ") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%lf", &user->solparams.interfacewidth);
        }    
        if (strstr(tok, "temperature ") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%lf", &user->solparams.temperature);
        }    
        if (strstr(tok, "initblocksize ") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            for (PetscInt dim=0; dim<user->dim; ++dim) {
                sscanf(strtok_r(NULL, " ", &savemtok), "%d", &user->amrparams.initblocksize[dim]);
            }
        }    
        if (strstr(tok, "initrefine ") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%d", &user->amrparams.initrefine);
        }    
        if (strstr(tok, "initcoarsen ") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%d", &user->amrparams.initcoarsen);
        }    
        if (strstr(tok, "maxnrefine ") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%d", &user->amrparams.maxnrefine);
        }    
        if (strstr(tok, "minnrefine ") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%d", &user->amrparams.minnrefine);
        }    
        if (strstr(tok, "amrinterval ") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%d", &user->amrparams.amrinterval);
        }    
        if (strstr(tok, "reltol ") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%lf", &user->solparams.reltol);
        }
        if (strstr(tok, "abstol ") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%lf", &user->solparams.abstol);
        }
        if (strstr(tok, "outputfreq ") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%d", &user->solparams.outputfreq);
        }    
        if (strstr(tok, "outfile ") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%s", user->solparams.outfile);
        }    
        if (strstr(tok, "petscoptions ") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            mtok = strtok_r(NULL, " ", &savemtok);
            while (mtok != NULL) {
                strcat(user->solparams.petscoptions,mtok);
                strcat(user->solparams.petscoptions," ");
                mtok = strtok_r(NULL, " ", &savemtok);
            }
        }    
        //advance the token
        tok = strtok_r(NULL, "\n", &savetok);
    }
    strcat(user->solparams.petscoptions," -dm_p4est_brick_bounds ");
    for (PetscInt dim=0; dim<user->dim; ++dim) {
        char strval[128];
        sprintf(strval,"0.0,%lf",user->size[dim]);
        strcat(user->solparams.petscoptions,strval);
        if (dim <user->dim-1) strcat(user->solparams.petscoptions,",");
    }
    strcat(user->solparams.petscoptions," -dm_p4est_brick_size ");
    for (PetscInt dim=0; dim<user->dim; ++dim) {
        char strval[128];
        sprintf(strval,"%d",user->amrparams.initblocksize[dim]);
        strcat(user->solparams.petscoptions,strval);
        if (dim <user->dim-1) strcat(user->solparams.petscoptions,",");
    }
    for (PetscInt dim=0; dim<user->dim; ++dim) {
        assert(user->amrparams.initblocksize[dim]*FastPow(2,user->amrparams.initrefine) == user->resolution[dim]);
    }
    free(substr);    
    free(buffer);   
    PetscFunctionReturn(0);
}

/*
 SetUpInterface - Import interface types from interface file
 */
PetscErrorCode SetUpInterface(AppCtx *user)
{
    char           interfacefile[PETSC_MAX_PATH_LEN] = "";
    char           *buffer = 0, *substr = 0;
    PetscMPIInt    rank;
    
    PetscFunctionBeginUser;  
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    PetscOptionsGetString(NULL,NULL,"-interfacefile",interfacefile,PETSC_MAX_PATH_LEN,NULL);
    FILE *infile=NULL;
    if (rank == 0) {
        infile = fopen (interfacefile, "r");
        if (infile==NULL) {
            printf("Error: input file can't be opened. \n");
            printf("Error: Please make sure the file %s is in the CWD.\n",interfacefile);
            return 1;
        }
    }
    PetscInt filesize;
    if (infile) {
        fseek(infile, 0, SEEK_END);
        filesize = ftell(infile);
        fseek(infile, 0, SEEK_SET);
    }    
    MPI_Bcast(&filesize,1,MPIU_INT ,0,PETSC_COMM_WORLD);
    buffer = malloc(filesize);
    if (rank == 0) {
        fread(buffer, 1, filesize, infile);
        fclose(infile);
    }
    MPI_Bcast(buffer,filesize,MPI_CHAR,0,PETSC_COMM_WORLD);
    
    char *tok, *savetok;
    substr = Extract(buffer, "<header>", "</header>");
    tok = strtok_r(substr, "\n", &savetok);
    while (tok !=NULL) {
        /* number of interface families */
        if (strstr(tok, "n_interfaces") != NULL) {
            char stra[128], strb[128];
            sscanf(tok, "%s %s", stra,strb); user->nf = atoi(strb);
            user->interface = (INTERFACE *) malloc(user->nf*sizeof(struct INTERFACE));
            INTERFACE *currentinterface = &user->interface[0];
            for (PetscInt interface=0; interface<user->nf; interface++,currentinterface++) {
                currentinterface->energy = 0.0;
                currentinterface->mobility = 0.0;
                currentinterface->potential = malloc(user->ncp*sizeof(PetscReal));
                currentinterface->mobilityc = malloc(user->ncp*sizeof(PetscReal));
                memset(currentinterface->potential,0,user->ncp*sizeof(PetscReal));
                memset(currentinterface->mobilityc,0,user->ncp*sizeof(PetscReal));
            }    
                
        }
        //advance the token
        tok = strtok_r(NULL, "\n", &savetok);
    }
    assert(user->nf >= 0);
    free(substr);
    
    INTERFACE *currentinterface = &user->interface[0];
    for (PetscInt interface=0; interface<user->nf; interface++,currentinterface++) {
        char starttag[128], endtag[128];
        sprintf(starttag, "<interface %d>",interface+1);        
        sprintf(endtag  , "</interface %d>",interface+1);        
        substr = Extract(buffer,starttag,endtag);
        assert(substr != NULL);
        tok = strtok_r(substr, "\n", &savetok);
        while (tok !=NULL) {
            char stra[128], strb[128];
            if (strstr(tok, "energy") != NULL) {
                sscanf(tok, "%s %s", stra,strb); currentinterface->energy = atof(strb);
                assert(currentinterface->energy >= 0.0);
            }
            if (strstr(tok, "mobility") != NULL) {
                sscanf(tok, "%s %s", stra,strb); currentinterface->mobility = atof(strb);
                assert(currentinterface->mobility >= 0.0);
            }
            if (strstr(tok, "potential") != NULL) {
                char *mtok, *savemtok;
                mtok = strtok_r(tok, " ", &savemtok);
                mtok = strtok_r(NULL, " ", &savemtok);
                PetscInt c = 0;
                while (mtok != NULL) {
                    PetscReal potential;
                    sscanf(mtok, "%lf", &potential);
                    currentinterface->potential[c] = potential;
                    mtok = strtok_r(NULL, " ", &savemtok);
                    c++;
                }
                assert(c == user->ncp);
            }
            if (strstr(tok, "mobilityc") != NULL) {
                char *mtok, *savemtok;
                mtok = strtok_r(tok, " ", &savemtok);
                mtok = strtok_r(NULL, " ", &savemtok);
                PetscInt c = 0;
                while (mtok != NULL) {
                    PetscReal mobilityc;
                    sscanf(mtok, "%lf", &mobilityc);
                    currentinterface->mobilityc[c] = mobilityc;
                    mtok = strtok_r(NULL, " ", &savemtok);
                    c++;
                }
                assert(c == user->ncp);
            }
            //advance the token
            tok = strtok_r(NULL, "\n", &savetok);
        }
        free(substr);
    }
    user->interfacelist = malloc(user->npf*user->npf*sizeof(uint16_t));
    
    substr = Extract(buffer, "<interface_mapping>", "</interface_mapping>");
    tok = strtok_r(substr, "\n", &savetok);
    PetscInt ctr=0;
    while (tok != NULL && ctr < user->npf*user->npf) {
        // process the line
        if (strstr(tok, "of") != NULL) {
            char stra[128], strb[128];
            PetscInt sof, interface;
            sscanf(tok, "%s of %s", stra, strb);
            sof = atoi(stra); interface = atoi(strb);
            for (PetscInt j=ctr;j<ctr+sof;j++) {
                PetscInt row = j/user->npf, col = j%user->npf;
                if (col > row) user->interfacelist[j] = (unsigned char) interface-1;
            }
            ctr+=sof;
        } else if (strstr(tok, "to") != NULL) {
            char stra[128], strb[128];
            PetscInt interface, interfacefrom, interfaceto;
            sscanf(tok, "%s to %s", stra, strb);
            interfacefrom = atoi(stra); interfaceto = atoi(strb);
            if (interfacefrom < interfaceto) {
                for (interface=interfacefrom;interface<=interfaceto;interface++) {
                    PetscInt row = ctr/user->npf, col = ctr%user->npf;
                    if (col > row) user->interfacelist[ctr] = (unsigned char) interface-1;
                    ctr++;
                }    
            } else {
                for (interface=interfacefrom;interface>=interfaceto;interface--) {
                    PetscInt row = ctr/user->npf, col = ctr%user->npf;
                    if (col > row) user->interfacelist[ctr] = (unsigned char) interface-1;
                    ctr++;
                }    
            }   
        } else if (atoi(tok)) {
            PetscInt interface = atoi(tok);
            PetscInt row = ctr/user->npf, col = ctr%user->npf;
            if (col > row) user->interfacelist[ctr] = (unsigned char) interface-1;
            ctr++;
        }
        //advance the token
        tok = strtok_r(NULL, "\n", &savetok);
    }
    assert(ctr == user->npf*user->npf);
    free(substr);
    free(buffer);    
    PetscFunctionReturn(0);
}

/*
 SetUpProblem - initializes solution
 */
PetscErrorCode SetUpProblem(Vec solution, AppCtx *user)
{
    PetscErrorCode    ierr;
    Vec               solutionl;
    PetscScalar       *fdof, *fdofl, *offset, *offsetg;
    PetscInt          localcell, cell, face, phase, g;
    PetscInt          conesize, supp, nsupp;
    const PetscInt    *cone, *scells;
    DMLabel           plabel = NULL;
    uint16_t          gslist[AS_SIZE], nalist[AS_SIZE];
    PetscScalar       *pcell, *dcell;
    uint16_t          setunion[AS_SIZE], injectionL[AS_SIZE], injectionR[AS_SIZE];
    MATERIAL          *currentmaterial;

    PetscFunctionBeginUser;  
    ierr = DMGetLabel(user->da_solution, "phase", &plabel);
    
    /* Determine active phase set */
    ierr = VecGetArray(solution, &fdof); CHKERRQ(ierr);
    for (localcell = 0; localcell < user->nlocalcells; ++localcell) {
        cell = user->localcells[localcell];

        /* set cell fields */
        ierr = DMLabelGetValue(plabel, cell, &phase);
        offset = NULL;
        ierr = DMPlexPointGlobalRef(user->da_solution, cell, fdof, &offset); CHKERRQ(ierr);
        offset[AS_OFFSET+0] = 1.0; offset[AS_OFFSET+1] = (PetscScalar) phase;
    }        
    ierr = VecRestoreArray(solution, &fdof);
    ierr = DMGetLocalVector(user->da_solution,&solutionl); CHKERRQ(ierr);  
    ierr = DMGlobalToLocalBegin(user->da_solution,solution,INSERT_VALUES,solutionl); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->da_solution,solution,INSERT_VALUES,solutionl); CHKERRQ(ierr);

    /* Determine superset with neighbouring active phase sets */
    ierr = VecGetArray(solution,&fdof); CHKERRQ(ierr);
    ierr = VecGetArray(solutionl,&fdofl); CHKERRQ(ierr);
    for (localcell = 0; localcell < user->nlocalcells; ++localcell) {
        cell = user->localcells[localcell];

        /* get cell state */
        offsetg = NULL;
        ierr = DMPlexPointGlobalRef(user->da_solution, cell, fdof, &offsetg); CHKERRQ(ierr);
        F2IFUNC(gslist,&offsetg[AS_OFFSET]);

        /* add union of neighbouring cells */
        ierr = DMPlexGetConeSize(user->da_solution,cell,&conesize);
        ierr = DMPlexGetCone    (user->da_solution,cell,&cone    );
        for (face=0; face<conesize; face++) {
            ierr = DMPlexGetSupportSize(user->da_solution, cone[face], &nsupp);
            ierr = DMPlexGetSupport(user->da_solution, cone[face], &scells);
            for (supp=0; supp<nsupp; supp++) {
                if (scells[supp] != cell && scells[supp] < user->ninteriorcells) {
                    offset = NULL;
                    ierr = DMPlexPointLocalRef(user->da_solution, scells[supp], fdofl, &offset); CHKERRQ(ierr);
                    F2IFUNC(nalist,&offset[AS_OFFSET]);
                    SetUnion(setunion,injectionL,injectionR,gslist,nalist);
                    memcpy(gslist,setunion,(setunion[0]+1)*sizeof(uint16_t));
                }
            }
        }
        I2FFUNC(&offsetg[AS_OFFSET],gslist);
    }    
    ierr = VecRestoreArray(solutionl,&fdofl); CHKERRQ(ierr);
    ierr = VecRestoreArray(solution,&fdof); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(user->da_solution,&solutionl); CHKERRQ(ierr);  

    /* Set initial phase field values */
    ierr = VecGetArray(solution, &fdof);
    for (localcell = 0; localcell < user->nlocalcells; ++localcell) {
        cell = user->localcells[localcell];

        /* get cell state */
        offset = NULL;
        ierr = DMPlexPointGlobalRef(user->da_solution, cell, fdof, &offset); CHKERRQ(ierr);
        F2IFUNC(gslist,&offset[AS_OFFSET]);
        pcell  = &offset[PF_OFFSET];
        dcell  = &offset[DP_OFFSET];
        PetscInt phase;
        ierr = DMLabelGetValue(plabel, cell, &phase);
    
        /* set initial conditions */
        for (g=0; g<gslist[0]; g++) {
            currentmaterial = &user->material[user->phasematerialmapping[gslist[g+1]]];
            ChemicalpotentialImplicit(&dcell[g*user->ndp],currentmaterial->c0,user->solparams.temperature,gslist[g+1],user);
            if (gslist[g+1] == phase) {
                pcell[g] = 1.0;
            } else {
                pcell[g] = 0.0;
            }
        }
    }        
    ierr = VecRestoreArray(solution, &fdof);
    PetscFunctionReturn(0);
}

