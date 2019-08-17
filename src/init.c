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
        }
        if (strstr(tok, "n_phases") != NULL) {
            char *mtok, *savemtok;
            mtok = strtok_r(tok, " ", &savemtok);
            sscanf(strtok_r(NULL, " ", &savemtok), "%d", &user->np);
            user->phasematerialmapping = malloc(user->np*sizeof(uint16_t));
            user->phasevoxelmapping = malloc(user->np*sizeof(struct roaring_bitmap_t *));
            for (PetscInt phase=0; phase<user->np; phase++) user->phasevoxelmapping[phase] = roaring_bitmap_create();
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
        //advance the token
        tok = strtok_r(NULL, "\n", &savetok);
    }
    currentmaterial = &user->material[0];
    for (PetscInt m=0; m<user->nmat; m++,currentmaterial++) {
        currentmaterial->c0 = malloc(user->nc*sizeof(PetscScalar));
        memset(currentmaterial->c0,0,user->nc*sizeof(PetscScalar));
        QUAD *currentquad = &currentmaterial->energy.quad;
        currentquad->ceq = malloc(user->nc*sizeof(PetscScalar));
        currentquad->unary = malloc(user->nc*sizeof(PetscScalar));
        currentquad->binary = malloc(user->nc*sizeof(PetscScalar));
        currentquad->mobilityc = malloc(user->nc*sizeof(PetscScalar));
        CALPHAD *currentcalphad = &currentmaterial->energy.calphad;
        currentcalphad->unary = malloc(user->nc*sizeof(PetscScalar));
        currentcalphad->binary = (RK *) malloc(user->nc*(user->nc-1)*sizeof(struct RK)/2);
        currentcalphad->ternary = (RK *) malloc(user->nc*(user->nc-1)*(user->nc-2)*sizeof(struct RK)/6);
        currentcalphad->mobilityc = malloc(user->nc*sizeof(PetscScalar));
    }
    assert(user->resolution[0] > 0  );
    assert(user->resolution[1] > 0  );
    assert(user->resolution[2] > 0  );
    assert(user->nc         >  0    );
    assert(user->nc         <= MAXCP);
    assert(user->np         >  0    );
    assert(user->nmat       >  0    );
    free(substr);
    
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
                    currentmaterial->model = QUADRATIC_FE;
                    currentquad->ref = 0.0;
                    memset(currentquad->ceq,0,user->nc*sizeof(PetscScalar));
                    memset(currentquad->unary,0,user->nc*sizeof(PetscScalar));
                    memset(currentquad->binary,0,user->nc*sizeof(PetscScalar));
                    memset(currentquad->mobilityc,0,user->nc*sizeof(PetscScalar));
                } else if (strstr(mtok, "calphaddis") != NULL) {
                    currentmaterial->model = CALPHAD_FE;
                    currentcalphad->ref = 0.0;
                    currentcalphad->RT = 2436.002; // default is room temperature
                    memset(currentcalphad->unary,0,user->nc*sizeof(PetscScalar));
                    memset(currentcalphad->mobilityc,0,user->nc*sizeof(PetscScalar));
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
                    currentmaterial->model = NONE_FE;
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
                    PetscScalar mobilityc;
                    sscanf(mtok, "%lf", &mobilityc);
                    if (currentmaterial->model == QUADRATIC_FE) {
                        currentquad->mobilityc[c] = mobilityc;
                    } else if (currentmaterial->model == CALPHAD_FE) {
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
                    PetscScalar c0;
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
                    PetscScalar ceq;
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
                    PetscScalar refenthalpy;
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
                    PetscScalar unaryenthalpy;
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
                    PetscScalar binaryenthalpy;
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
                    PetscScalar refenthalpy;
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
                    PetscScalar unaryenthalpy;
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
                currentbinary->enthalpy = malloc(currentbinary->n*sizeof(PetscScalar));
                memset(currentbinary->enthalpy,0,currentbinary->n*sizeof(PetscScalar));
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
                    PetscScalar binaryenthalpy;
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
                currentternary->enthalpy = malloc(currentternary->n*sizeof(PetscScalar));
                memset(currentternary->enthalpy,0,currentternary->n*sizeof(PetscScalar));
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
                    PetscScalar ternaryenthalpy;
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
        assert(currentmaterial->molarvolume > 0.0);
        for (PetscInt c=0;c<user->nc;c++) {
            assert(currentmaterial->c0[c] > 0.0);
        }
        if        (currentmaterial->model == QUADRATIC_FE) {
            QUAD *currentquad = &currentmaterial->energy.quad;
            for (PetscInt c=0;c<user->nc;c++) {
                assert(currentquad->ceq[c] > 0.0);
                assert(currentquad->mobilityc[c] > 0.0);
            }
        } else if (currentmaterial->model ==   CALPHAD_FE) {
            CALPHAD *currentcalphad = &currentmaterial->energy.calphad;
            assert(currentcalphad->RT > 0.0);
            for (PetscInt c=0;c<user->nc;c++) {
                assert(currentcalphad->mobilityc[c] > 0.0);
            }
        } else if (currentmaterial->model ==      NONE_FE) {
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
                roaring_bitmap_add(user->phasevoxelmapping[p-1], j);
            ctrv+=sof;
        } else if (strstr(tok, "to") != NULL) {
            char stra[128], strb[128];
            PetscInt p, pfrom, pto;
            sscanf(tok, "%s to %s", stra, strb);
            pfrom = atoi(stra); pto = atoi(strb);
            if (pfrom < pto) {
                for (p=pfrom;p<=pto;p++) 
                    roaring_bitmap_add(user->phasevoxelmapping[p-1], ctrv++);
            } else {
                for (p=pfrom;p>=pto;p--)
                    roaring_bitmap_add(user->phasevoxelmapping[p-1], ctrv++);
            }   
        }
        else {
            roaring_bitmap_add(user->phasevoxelmapping[atoi(tok)-1], ctrv++);
        }
        //advance the token
        tok = strtok_r(NULL, "\n", &savetok);
    }
    assert(ctrv == user->resolution[0]*user->resolution[1]*user->resolution[2] && tok == NULL);
    for (PetscInt phase=0; phase<user->np; phase++) 
        roaring_bitmap_run_optimize(user->phasevoxelmapping[phase]);
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
                currentinterface->potential = malloc(user->nc*sizeof(PetscScalar));
                currentinterface->mobilityc = malloc(user->nc*sizeof(PetscScalar));
                memset(currentinterface->potential,0,user->nc*sizeof(PetscScalar));
                memset(currentinterface->mobilityc,0,user->nc*sizeof(PetscScalar));
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
                    PetscScalar potential;
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
                    PetscScalar mobilityc;
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
PetscErrorCode SetUpProblem(DM da,AppCtx *user,Vec X)
{
    PetscInt       xs, ys, zs, xm, ym, zm;
    PetscErrorCode ierr;
    Vec            gactivephases;
    F_DOFS         ***fdof;
    C_DOFS         ***cdof;
    F2I            ***nalist, ***alist, ***glist;
    MATERIAL       *currentmaterial;
    uint16_t       superset[MAXAP];
    
    PetscFunctionBeginUser;  
    /*
     Import initial microstructure from geom file
     */
    ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
    zm+=zs; ym+=ys; xm+=xs; 

    /* initially only self phase is active... */
    ierr = DMGetGlobalVector(user->daIdx,&gactivephases); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->daIdx,gactivephases,&glist); CHKERRQ(ierr);
    for (PetscInt phase=0;phase<user->np;phase++) {
        roaring_bitmap_t *tempbm = roaring_bitmap_copy(user->phasevoxelmapping[phase]);
        roaring_uint32_iterator_t *g = roaring_create_iterator(user->phasevoxelmapping[phase]);
        while (g->has_value) {
            PetscInt i = (g->current_value%(user->resolution[0]));
            PetscInt j = (g->current_value%(user->resolution[0]*user->resolution[1]))/user->resolution[0];
            PetscInt k = (g->current_value/(user->resolution[1]*user->resolution[0]));
            if ((xs <= i && i < xm) && (ys <= j && j < ym) && (zs <= k && k < zm)) {
                glist[k][j][i].i[0]=1;
                glist[k][j][i].i[glist[k][j][i].i[0]] = phase;
            } else {       
                /* remove off processor nodes... */
                roaring_bitmap_remove(tempbm, g->current_value);
            }
            roaring_advance_uint32_iterator(g);
        }
        roaring_free_uint32_iterator(g);
        user->phasevoxelmapping[phase] = roaring_bitmap_copy(tempbm);
        roaring_bitmap_free(tempbm);
    }
    ierr = DMDAVecRestoreArray(user->daIdx,gactivephases,&glist); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(user->daIdx,gactivephases,INSERT_VALUES,user->activephases); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->daIdx,gactivephases,INSERT_VALUES,user->activephases); CHKERRQ(ierr);
    ierr = DMRestoreGlobalVector(user->daIdx,&gactivephases); CHKERRQ(ierr);  

    /* activate self phase on neighbours... */
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

    /* Set initial phase field values */
    ierr = DMDAVecGetArray(da,X,&fdof); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->daCmp,user->cvec,&cdof); CHKERRQ(ierr);
    for (PetscInt phase=0;phase<user->np;phase++) {
        currentmaterial = &user->material[user->phasematerialmapping[phase]];
        roaring_uint32_iterator_t *g = roaring_create_iterator(user->phasevoxelmapping[phase]);
        while (g->has_value) {
            PetscInt i = (g->current_value%(user->resolution[0]));
            PetscInt j = (g->current_value%(user->resolution[0]*user->resolution[1]))/user->resolution[0];
            PetscInt k = (g->current_value/(user->resolution[1]*user->resolution[0]));
            fdof[k][j][i].p[0] = 1.0;
            uint16_t plist[2] = {1,phase};
            memcpy(cdof[k][j][i].c,currentmaterial->c0,user->nc*sizeof(PetscScalar)); 
            Chemicalpotential(fdof[k][j][i].m,cdof[k][j][i].c,fdof[k][j][i].p,plist,user);
            roaring_advance_uint32_iterator(g);
        }
        roaring_free_uint32_iterator(g);
    }    
    ierr = DMDAVecGetArray(user->daCmp,user->cvec,&cdof); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,X,&fdof); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

