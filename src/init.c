/**
 * Initialisation module file
 */
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <petscdm.h>
#include <petscdmda.h>
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
    user->resolution[0] = 1;
    user->resolution[1] = 1;
    user->resolution[2] = 1;
    user->nc = 1;
    user->np = 1;
    user->nmat = 1;
    user->interpolation = LINEAR_INTERPOLATION;
    while (tok !=NULL) {
        // process the line
        if (strstr(tok, "grid") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            mtok = strtok_r(NULL, " ", &savemtok);
            user->resolution[0] = atoi(strtok_r(NULL, " ", &savemtok));
            mtok = strtok_r(NULL, " ", &savemtok);
            user->resolution[1] = atoi(strtok_r(NULL, " ", &savemtok));
            mtok = strtok_r(NULL, " ", &savemtok);
            user->resolution[2] = atoi(strtok_r(NULL, " ", &savemtok));
            user->phasevoxelmapping = malloc(user->resolution[2]*user->resolution[1]*user->resolution[0]*sizeof(uint16_t));
        }
        if (strstr(tok, "n_phases") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%d", &user->np);
            user->phasematerialmapping = malloc(user->np*sizeof(uint16_t));
        }    
        if (strstr(tok, "n_materials") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%d", &user->nmat);
            user->material = (MATERIAL *) malloc(user->nmat*sizeof(struct MATERIAL));
        }
        if (strstr(tok, "n_components") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%d", &user->nc);
            user->componentname = (char **) malloc(user->nc*sizeof(char *));
            for (PetscInt c=0; c<user->nc; c++) {
                user->componentname[c] = (char *) malloc(2*sizeof(char));
                sprintf(user->componentname[c],"%d",c+1); 
            }
        }    
        if (strstr(tok, "componentnames") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            for (PetscInt c=0; c<user->nc; c++) {
                sscanf(strtok_r(NULL, " ", &savemtok), "%s", user->componentname[c]);
            }
        }
        /* interpolation type */
        if (strstr(tok, "interpolation_type") != NULL) {
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
    assert(user->resolution[0] > 0  );
    assert(user->resolution[1] > 0  );
    assert(user->resolution[2] > 0  );
    assert(user->nc         >  0    );
    assert(user->nc         <= MAXCP);
    assert(user->np         >  0    );
    assert(user->nmat       >  0    );
    assert(user->interpolation != NONE_INTERPOLATION);
    free(substr);
    
    /* initialise material information */
    currentmaterial = &user->material[0];
    for (PetscInt m=0; m<user->nmat; m++,currentmaterial++) {
        currentmaterial->model = NONE_CHEMENERGY;
        currentmaterial->c0 = malloc(user->nc*sizeof(PetscReal));
        memset(currentmaterial->c0,0,user->nc*sizeof(PetscReal));
        QUAD *currentquad = &currentmaterial->energy.quad;
        currentquad->ceq = malloc(user->nc*sizeof(PetscReal));
        currentquad->unary = malloc(user->nc*sizeof(PetscReal));
        currentquad->binary = malloc(user->nc*sizeof(PetscReal));
        currentquad->mobilityc = malloc(user->nc*sizeof(PetscReal));
        CALPHAD *currentcalphad = &currentmaterial->energy.calphad;
        currentcalphad->unary = malloc(user->nc*sizeof(PetscReal));
        currentcalphad->binary = (RK *) malloc(user->nc*(user->nc-1)*sizeof(struct RK)/2);
        currentcalphad->ternary = (RK *) malloc(user->nc*(user->nc-1)*(user->nc-2)*sizeof(struct RK)/6);
        currentcalphad->mobilityc = malloc(user->nc*sizeof(PetscReal));
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
            if (strstr(tok, "chemicalenergy_type") != NULL) {
                char *mtok, *savemtok;
                mtok = strtok_r(tok, " ", &savemtok);
                mtok = strtok_r(NULL, " ", &savemtok);
                if (strstr(mtok, "quadratic") != NULL) {
                    currentmaterial->model = QUADRATIC_CHEMENERGY;
                    currentquad->ref = 0.0;
                    memset(currentquad->ceq,0,user->nc*sizeof(PetscReal));
                    memset(currentquad->unary,0,user->nc*sizeof(PetscReal));
                    memset(currentquad->binary,0,user->nc*sizeof(PetscReal));
                    memset(currentquad->mobilityc,0,user->nc*sizeof(PetscReal));
                } else if (strstr(mtok, "calphaddis") != NULL) {
                    currentmaterial->model = CALPHAD_CHEMENERGY;
                    currentcalphad->ref = 0.0;
                    currentcalphad->RT = 2436.002; // default is room temperature
                    memset(currentcalphad->unary,0,user->nc*sizeof(PetscReal));
                    memset(currentcalphad->mobilityc,0,user->nc*sizeof(PetscReal));
                    RK *currentbinary = &currentcalphad->binary[0];
                    RK *currentternary = &currentcalphad->ternary[0];
                    for (PetscInt ck=0; ck<user->nc; ck++) {
                        for (PetscInt cj=ck+1; cj<user->nc; cj++,currentbinary++) {
                            currentbinary->n = 0;
                            for (PetscInt ci=cj+1; ci<user->nc; ci++,currentternary++) { 
                                currentternary->n = 0;
                            }    
                        }
                    } 
                } else {
                    currentmaterial->model = NONE_CHEMENERGY;
                    currentmaterial->c0[0] = 1.0;
                }
            }
            /* temperature */
            if (strstr(tok, "RT") != NULL) {
                char *mtok, *savemtok;
                mtok = strtok_r(tok, " ", &savemtok);
                sscanf(strtok_r(NULL, " ", &savemtok), "%lf", &currentcalphad->RT);
            }
            /* molar volume */
            if (strstr(tok, "molarvolume") != NULL) {
                char *mtok, *savemtok;
                mtok = strtok_r(tok, " ", &savemtok);
                sscanf(strtok_r(NULL, " ", &savemtok), "%lf", &currentmaterial->molarvolume);
            }
            /* component mobility for each material */
            if (strstr(tok, "mobilityc") != NULL) {
                char *mtok, *savemtok;
                mtok = strtok_r(tok, " ", &savemtok);
                mtok = strtok_r(NULL, " ", &savemtok);
                PetscInt c = 0;
                while (mtok != NULL) {
                    PetscReal mobilityc;
                    sscanf(mtok, "%lf", &mobilityc);
                    if (currentmaterial->model == QUADRATIC_CHEMENERGY) {
                        currentquad->mobilityc[c] = mobilityc;
                    } else if (currentmaterial->model == CALPHAD_CHEMENERGY) {
                        currentcalphad->mobilityc[c] = mobilityc;
                    }   
                    mtok = strtok_r(NULL, " ", &savemtok);
                    c++;
                }
                assert(c == user->nc);
            }
            /* initial composition */
            if (strstr(tok, "c0") != NULL) {
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
                assert(c == user->nc);
            }
            /* equillibrium composition */
            if (strstr(tok, "quad_ceq") != NULL) {
                char *mtok, *savemtok;
                mtok = strtok_r(tok, " ", &savemtok);
                mtok = strtok_r(NULL, " ", &savemtok);
                PetscInt c=0;
                while (mtok != NULL) {
                    PetscReal ceq;
                    sscanf(mtok, "%lf", &ceq);
                    assert(ceq > 0.0);
                    currentquad->ceq[c] = ceq;
                    mtok = strtok_r(NULL, " ", &savemtok);
                    c++;
                }
                assert(c == user->nc);
            }
            /* F_chem parameters for each material */
            if (strstr(tok, "quad_refenthalpy") != NULL) {
                char *mtok, *savemtok;
                mtok = strtok_r(tok, " ", &savemtok);
                mtok = strtok_r(NULL, " ", &savemtok);
                while (mtok != NULL) {
                    PetscReal refenthalpy;
                    sscanf(mtok, "%lf", &refenthalpy);
                    currentquad->ref = refenthalpy;
                    mtok = strtok_r(NULL, " ", &savemtok);
                }
            }
            /* F_chem parameters for each material */
            if (strstr(tok, "quad_unaryenthalpy") != NULL) {
                char *mtok, *savemtok;
                mtok = strtok_r(tok, " ", &savemtok);
                mtok = strtok_r(NULL, " ", &savemtok);
                PetscInt c = 0;
                while (mtok != NULL) {
                    PetscReal unaryenthalpy;
                    sscanf(mtok, "%lf", &unaryenthalpy);
                    currentquad->unary[c] = unaryenthalpy;
                    mtok = strtok_r(NULL, " ", &savemtok);
                    c++;
                }
                assert(c == user->nc);
            }
            /* F_chem parameters for each material */
            if (strstr(tok, "quad_binaryenthalpy") != NULL) {
                char *mtok, *savemtok;
                mtok = strtok_r(tok, " ", &savemtok);
                mtok = strtok_r(NULL, " ", &savemtok);
                PetscInt c = 0;
                while (mtok != NULL) {
                    PetscReal binaryenthalpy;
                    sscanf(mtok, "%lf", &binaryenthalpy);
                    currentquad->binary[c] = binaryenthalpy;
                    mtok = strtok_r(NULL, " ", &savemtok);
                    c++;
                }
                assert(c == user->nc);
            }
            /* F_chem parameters for each material */
            if (strstr(tok, "calphad_refenthalpy") != NULL) {
                char *mtok, *savemtok;
                mtok = strtok_r(tok, " ", &savemtok);
                mtok = strtok_r(NULL, " ", &savemtok);
                while (mtok != NULL) {
                    PetscReal refenthalpy;
                    sscanf(mtok, "%lf", &refenthalpy);
                    currentcalphad->ref = refenthalpy;
                    mtok = strtok_r(NULL, " ", &savemtok);
                }
            }
            /* F_chem parameters for each material */
            if (strstr(tok, "calphad_unaryenthalpy") != NULL) {
                char *mtok, *savemtok;
                mtok = strtok_r(tok, " ", &savemtok);
                mtok = strtok_r(NULL, " ", &savemtok);
                PetscInt c = 0;
                while (mtok != NULL) {
                    PetscReal unaryenthalpy;
                    sscanf(mtok, "%lf", &unaryenthalpy);
                    currentcalphad->unary[c] = unaryenthalpy;
                    mtok = strtok_r(NULL, " ", &savemtok);
                    c++;
                }
                assert(c == user->nc);
            }
            /* F_chem parameters for each material */
            if (strstr(tok, "calphad_nbinaryenthalpy") != NULL) {
                char *mtok, *savemtok;
                mtok = strtok_r( tok, " ", &savemtok);
                PetscInt cj = atoi(strtok_r(NULL, " ", &savemtok))-1;
                PetscInt ci = atoi(strtok_r(NULL, " ", &savemtok))-1;
                assert(cj < ci && ci < user->nc);
                PetscInt offset = 0;
                for (PetscInt cjj=0;cjj<user->nc;cjj++) {
                    for (PetscInt cii=cjj+1;cii<user->nc;cii++,offset++) {
                        if (cjj == cj && cii == ci) break;
                    }
                }
                RK *currentbinary = &currentcalphad->binary[offset];
                mtok = strtok_r(NULL, " ", &savemtok);
                sscanf(mtok, "%d", &currentbinary->n);
                assert(currentbinary->n >= 0);
                currentbinary->enthalpy = malloc(currentbinary->n*sizeof(PetscReal));
                memset(currentbinary->enthalpy,0,currentbinary->n*sizeof(PetscReal));
            }
            /* F_chem parameters for each material */
            if (strstr(tok, "calphad_binaryenthalpy") != NULL) {
                char *mtok, *savemtok;
                mtok = strtok_r( tok, " ", &savemtok);
                PetscInt cj = atoi(strtok_r(NULL, " ", &savemtok))-1;
                PetscInt ci = atoi(strtok_r(NULL, " ", &savemtok))-1;
                assert(cj < ci && ci < user->nc);
                PetscInt offset = 0;
                for (PetscInt cjj=0;cjj<user->nc;cjj++) {
                    for (PetscInt cii=cjj+1;cii<user->nc;cii++,offset++) {
                        if (cjj == cj && cii == ci) break;
                    }
                }
                RK *currentbinary = &currentcalphad->binary[offset];
                mtok = strtok_r(NULL, " ", &savemtok);
                PetscInt nrk = 0;
                while (mtok != NULL) {
                    PetscReal binaryenthalpy;
                    sscanf(mtok, "%lf", &binaryenthalpy);
                    currentbinary->enthalpy[nrk] = binaryenthalpy;
                    mtok = strtok_r(NULL, " ", &savemtok);
                    nrk++;
                }
                assert(nrk == currentbinary->n);
            }
            /* F_chem parameters for each material */
            if (strstr(tok, "calphad_nternaryenthalpy") != NULL) {
                char *mtok, *savemtok;
                mtok = strtok_r( tok, " ", &savemtok);
                PetscInt ck = atoi(strtok_r(NULL, " ", &savemtok))-1;
                PetscInt cj = atoi(strtok_r(NULL, " ", &savemtok))-1;
                PetscInt ci = atoi(strtok_r(NULL, " ", &savemtok))-1;
                assert(ck < cj && cj < ci && ci < user->nc);
                PetscInt offset = 0;
                for (PetscInt ckk=0;ckk<user->nc;ckk++) {
                    for (PetscInt cjj=ckk+1;cjj<user->nc;cjj++) {
                        for (PetscInt cii=cjj+1;cii<user->nc;cii++,offset++) {
                            if (ckk == ck && cjj == cj && cii == ci) break;
                        }
                    }
                }
                RK *currentternary = &currentcalphad->ternary[offset];
                mtok = strtok_r(NULL, " ", &savemtok);
                sscanf(mtok, "%d", &currentternary->n);
                assert(currentternary->n >= 0 && currentternary->n <=3);
                currentternary->enthalpy = malloc(currentternary->n*sizeof(PetscReal));
                memset(currentternary->enthalpy,0,currentternary->n*sizeof(PetscReal));
                currentternary->i = malloc(currentternary->n*sizeof(PetscInt));
                memset(currentternary->i,0,currentternary->n*sizeof(PetscInt));
            }
            /* F_chem parameters for each material */
            if (strstr(tok, "calphad_iternaryenthalpy") != NULL) {
                char *mtok, *savemtok;
                mtok = strtok_r( tok, " ", &savemtok);
                PetscInt ck = atoi(strtok_r(NULL, " ", &savemtok))-1;
                PetscInt cj = atoi(strtok_r(NULL, " ", &savemtok))-1;
                PetscInt ci = atoi(strtok_r(NULL, " ", &savemtok))-1;
                assert(ck < cj && cj < ci && ci < user->nc);
                PetscInt offset = 0;
                for (PetscInt ckk=0;ckk<user->nc;ckk++) {
                    for (PetscInt cjj=ckk+1;cjj<user->nc;cjj++) {
                        for (PetscInt cii=cjj+1;cii<user->nc;cii++,offset++) {
                            if (ckk == ck && cjj == cj && cii == ci) break;
                        }
                    }
                }            
                RK *currentternary = &currentcalphad->ternary[offset];
                mtok = strtok_r(NULL, " ", &savemtok);
                PetscInt nrk = 0;
                while (mtok != NULL) {
                    PetscInt iternaryenthalpy;
                    sscanf(mtok, "%d", &iternaryenthalpy);
                    currentternary->i[nrk] = iternaryenthalpy-1;
                    assert(   currentternary->i[nrk] == ci 
                           || currentternary->i[nrk] == cj
                           || currentternary->i[nrk] == ck);
                    mtok = strtok_r(NULL, " ", &savemtok);
                    nrk++;
                }
                assert(nrk == currentternary->n);
            }
            /* F_chem parameters for each material */
            if (strstr(tok, "calphad_ternaryenthalpy") != NULL) {
                char *mtok, *savemtok;
                mtok = strtok_r( tok, " ", &savemtok);
                PetscInt ck = atoi(strtok_r(NULL, " ", &savemtok))-1;
                PetscInt cj = atoi(strtok_r(NULL, " ", &savemtok))-1;
                PetscInt ci = atoi(strtok_r(NULL, " ", &savemtok))-1;
                assert(ck < cj && cj < ci && ci < user->nc);
                PetscInt offset = 0;
                for (PetscInt ckk=0;ckk<user->nc;ckk++) {
                    for (PetscInt cjj=ckk+1;cjj<user->nc;cjj++) {
                        for (PetscInt cii=cjj+1;cii<user->nc;cii++,offset++) {
                            if (ckk == ck && cjj == cj && cii == ci) break;
                        }
                    }
                }            
                RK *currentternary = &currentcalphad->ternary[offset];
                mtok = strtok_r(NULL, " ", &savemtok);
                PetscInt nrk = 0;
                while (mtok != NULL) {
                    PetscReal ternaryenthalpy;
                    sscanf(mtok, "%lf", &ternaryenthalpy);
                    currentternary->enthalpy[nrk] = ternaryenthalpy;
                    mtok = strtok_r(NULL, " ", &savemtok);
                    nrk++;
                }
                assert(nrk == currentternary->n);
            }
            //advance the token
            tok = strtok_r(NULL, "\n", &savetok);
        }
        /* material information sanity checks */
        assert(currentmaterial->molarvolume > 0.0);
        for (PetscInt c=0;c<user->nc;c++) {
            assert(currentmaterial->c0[c] > 0.0);
        }
        if        (currentmaterial->model == QUADRATIC_CHEMENERGY) {
            QUAD *currentquad = &currentmaterial->energy.quad;
            for (PetscInt c=0;c<user->nc;c++) {
                assert(currentquad->ceq[c] > 0.0);
                assert(currentquad->mobilityc[c] > 0.0);
            }
        } else if (currentmaterial->model ==   CALPHAD_CHEMENERGY) {
            CALPHAD *currentcalphad = &currentmaterial->energy.calphad;
            assert(currentcalphad->RT > 0.0);
            for (PetscInt c=0;c<user->nc;c++) {
                assert(currentcalphad->mobilityc[c] > 0.0);
            }
        } else if (currentmaterial->model ==      NONE_CHEMENERGY) {
            assert(user->nc == 1);
        }
        free(substr);
    }
    
    substr = Extract(buffer, "<phase_material_mapping>", "</phase_material_mapping>");
    tok = strtok_r(substr, "\n", &savetok);
    PetscInt ctrm=0;
    while (tok != NULL && ctrm < user->np) {
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
        else {
            assert(atoi(tok) <= user->nmat);
            user->phasematerialmapping[ctrm++] = atoi(tok)-1;
        }
        //advance the token
        tok = strtok_r(NULL, "\n", &savetok);
    }
    assert(ctrm == user->np && tok == NULL);
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
        else {
            user->phasevoxelmapping[ctrv++] = atoi(tok)-1;
        }
        //advance the token
        tok = strtok_r(NULL, "\n", &savetok);
    }
    assert(ctrv == user->resolution[0]*user->resolution[1]*user->resolution[2] && tok == NULL);
    free(substr);    

    substr = Extract(buffer, "<solution_parameters>", "</solution_parameters>");
    tok = strtok_r(substr, "\n", &savetok);
    /* default solution parameters */
    user->params.finaltime = 1.0;
    user->params.timestep = TOL;
    user->params.mintimestep = TOL;
    user->params.maxtimestep = LARGE;
    user->params.step = 0;    
    user->params.interfacewidth = 5;
    user->params.feorder_phase = 1;
    user->params.feorder_chem = 1;
    user->params.reltol = 1e-6;
    user->params.abstol = 1e-6;
    user->params.outputfreq = 1;
    strncpy(user->params.outfile, "output", 128);
    /* initialise solution parameters */
    while (tok !=NULL) {
        // process the line
        if (strstr(tok, "finaltime") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%lf", &user->params.finaltime);
        }
        if (strstr(tok, "timestep0") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%lf", &user->params.timestep);
        }
        if (strstr(tok, "timestepmin") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%lf", &user->params.mintimestep);
        }
        if (strstr(tok, "timestepmax") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%lf", &user->params.maxtimestep);
        }
        if (strstr(tok, "interfacewidth") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%lf", &user->params.interfacewidth);
        }    
        if (strstr(tok, "feorder_phase") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%d", &user->params.feorder_phase);
        }    
        if (strstr(tok, "feorder_chem") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%d", &user->params.feorder_chem);
        }    
        if (strstr(tok, "reltol") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%lf", &user->params.reltol);
        }
        if (strstr(tok, "abstol") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%lf", &user->params.abstol);
        }
        if (strstr(tok, "outputfreq") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%d", &user->params.outputfreq);
        }    
        if (strstr(tok, "outfile") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%s", user->params.outfile);
        }    
        //advance the token
        tok = strtok_r(NULL, "\n", &savetok);
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
                currentinterface->potential = malloc(user->nc*sizeof(PetscReal));
                currentinterface->mobilityc = malloc(user->nc*sizeof(PetscReal));
                memset(currentinterface->potential,0,user->nc*sizeof(PetscReal));
                memset(currentinterface->mobilityc,0,user->nc*sizeof(PetscReal));
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
                assert(c == user->nc);
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
                assert(c == user->nc);
            }
            //advance the token
            tok = strtok_r(NULL, "\n", &savetok);
        }
        free(substr);
    }
    user->interfacelist = malloc(user->np*user->np*sizeof(uint16_t));
    
    substr = Extract(buffer, "<interface_mapping>", "</interface_mapping>");
    tok = strtok_r(substr, "\n", &savetok);
    PetscInt ctr=0;
    while (tok != NULL && ctr < user->np*user->np) {
        // process the line
        if (strstr(tok, "of") != NULL) {
            char stra[128], strb[128];
            PetscInt sof, interface;
            sscanf(tok, "%s of %s", stra, strb);
            sof = atoi(stra); interface = atoi(strb);
            for (PetscInt j=ctr;j<ctr+sof;j++) {
                PetscInt row = j/user->np, col = j%user->np;
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
                    PetscInt row = ctr/user->np, col = ctr%user->np;
                    if (col > row) user->interfacelist[ctr] = (unsigned char) interface-1;
                    ctr++;
                }    
            } else {
                for (interface=interfacefrom;interface>=interfaceto;interface--) {
                    PetscInt row = ctr/user->np, col = ctr%user->np;
                    if (col > row) user->interfacelist[ctr] = (unsigned char) interface-1;
                    ctr++;
                }    
            }   
        }
        else {
            PetscInt interface = atoi(tok);
            PetscInt row = ctr/user->np, col = ctr%user->np;
            if (col > row) user->interfacelist[ctr] = (unsigned char) interface-1;
            ctr++;
        }
        //advance the token
        tok = strtok_r(NULL, "\n", &savetok);
    }
    assert(ctr == user->np*user->np);
    free(substr);
    free(buffer);    
    PetscFunctionReturn(0);
}

/*
 SetUpProblem - initializes solution
 */
PetscErrorCode SetUpProblem(DM da_solution, Vec solution, AppCtx *user)
{
    PetscInt       xs, ys, zs, xm, ym, zm;
    PetscErrorCode ierr;
    Vec            gactivephaseset, gactivephasesuperset;
    FIELD          ***fdof;
    STATE          ***matstate;
    F2I            ***slist, ***alist, ***gslist, ***galist;
    MATERIAL       *currentmaterial;
    uint16_t       superset[MAXAP], injectionA[MAXAP], injectionB[MAXAP];
    
    PetscFunctionBeginUser;  
    /*
     Import initial microstructure from geom file
     */
    ierr = DMDAGetCorners(da_solution,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
    zm+=zs; ym+=ys; xm+=xs; 

    /* initially only self phase is active... */
    ierr = DMGetGlobalVector(user->da_phaseID,&gactivephaseset); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da_phaseID,gactivephaseset,&galist); CHKERRQ(ierr);
    for (PetscInt k=zs; k<zm; k++) {
        for (PetscInt j=ys; j<ym; j++) {
            for (PetscInt i=xs; i<xm; i++) {
                galist[k][j][i].i[0]=1;
                galist[k][j][i].i[1]=user->phasevoxelmapping[i+user->resolution[0]*(j+user->resolution[1]*k)];
            }
        }
    }        
    ierr = DMDAVecRestoreArray(user->da_phaseID,gactivephaseset,&galist); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(user->da_phaseID,gactivephaseset,INSERT_VALUES,user->activephaseset); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->da_phaseID,gactivephaseset,INSERT_VALUES,user->activephaseset); CHKERRQ(ierr);
    ierr = DMRestoreGlobalVector(user->da_phaseID,&gactivephaseset); CHKERRQ(ierr);  

    /* activate self phase on neighbours... */
    ierr = DMDAVecGetArrayRead(user->da_phaseID,user->activephaseset,&alist); CHKERRQ(ierr);
    ierr = DMGetGlobalVector(user->da_phaseID,&gactivephasesuperset); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da_phaseID,gactivephasesuperset,&gslist); CHKERRQ(ierr);
    for (PetscInt k=zs; k<zm; k++) {
        for (PetscInt j=ys; j<ym; j++) {
            for (PetscInt i=xs; i<xm; i++) {
                PetscInt ie=i+1,iw=i-1,jn=j+1,js=j-1,kt=k+1,kb=k-1;
                gslist[k][j][i].i[0] = 0;
                for (PetscInt kk=kb; kk<=kt; kk++) {
                    for (PetscInt jj=js; jj<=jn; jj++) {
                        for (PetscInt ii=iw; ii<=ie; ii++) {
                            SetUnion(superset,injectionA,injectionB,alist[kk][jj][ii].i,gslist[k][j][i].i);
                            memcpy(gslist[k][j][i].i,superset,2*MAXAP);          
                        }
                    }
                }
            }
        }
    }       
    ierr = DMDAVecRestoreArray(user->da_phaseID,gactivephasesuperset,&gslist); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(user->da_phaseID,gactivephasesuperset,INSERT_VALUES,user->activephasesuperset); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->da_phaseID,gactivephasesuperset,INSERT_VALUES,user->activephasesuperset); CHKERRQ(ierr);
    ierr = DMRestoreGlobalVector(user->da_phaseID,&gactivephasesuperset); CHKERRQ(ierr);  
    ierr = DMDAVecRestoreArrayRead(user->da_phaseID,user->activephaseset,&alist); CHKERRQ(ierr);

    /* Set initial phase field values */
    ierr = DMDAVecGetArray(da_solution,solution,&fdof); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da_matstate,user->matstate,&matstate); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da_phaseID,user->activephasesuperset,&slist); CHKERRQ(ierr);
    for (PetscInt k=zs; k<zm; k++) {
        for (PetscInt j=ys; j<ym; j++) {
            for (PetscInt i=xs; i<xm; i++) {
                for (PetscInt g=0; g<slist[k][j][i].i[0]; g++) {
                    currentmaterial = &user->material[user->phasematerialmapping[slist[k][j][i].i[g+1]]];
                    if (slist[k][j][i].i[g+1] == user->phasevoxelmapping[i+user->resolution[0]*(j+user->resolution[1]*k)]) {
                        fdof[k][j][i].p[g] = 1.0;
                    } else {
                        fdof[k][j][i].p[g] = 0.0;
                    }
                    memcpy(&matstate[k][j][i].c[g*user->nc],currentmaterial->c0,user->nc*sizeof(PetscReal)); 
                }       
                Chemicalpotential(fdof[k][j][i].m,matstate[k][j][i].c,fdof[k][j][i].p,slist[k][j][i].i,user);
            }
        }
    }        
    ierr = DMDAVecRestoreArrayRead(user->da_phaseID,user->activephasesuperset,&slist); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da_matstate,user->matstate,&matstate); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_solution,solution,&fdof); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

