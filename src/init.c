/**
 * Initialisation module file
 */
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <petscdm.h>
#include <petscdmplex.h>
#include <petscdmforest.h>
#include <petscviewerhdf5.h>
#include "material.h"
#include "utility.h"
#include "typedef.h"
#include <yaml.h>

/* Including my own header for checking by compiler */
#define INIT_IMPORT
#include "init.h"

/*
 GetProperty - Get property from YAML file buffer
 */
PetscErrorCode GetProperty(char **propval, PetscInt *propsize, 
                           const char *mappingname, const char *propname, 
                           const unsigned char *buffer, const PetscInt filesize)
{
    yaml_parser_t parser;
    yaml_event_t  event;
    PetscInt doclevel = 0, proplevel = 0;
    *propsize = 0;

    /* Initialize parser */
    if(!yaml_parser_initialize(&parser)) {
        PetscPrintf(PETSC_COMM_WORLD,"Failed to initialize parser\n");
        return 1;
    }    

    /* Set input buffer */
    yaml_parser_set_input_string(&parser, buffer, filesize);

    /* Parsing config file */
    char mappingname_tok[PETSC_MAX_PATH_LEN], *tok, *savetok;

    strcpy(mappingname_tok,mappingname);
    tok = strtok_r(mappingname_tok, "/", &savetok); proplevel++;
    while (tok != NULL) {
        while (yaml_parser_parse(&parser, &event) && event.type != YAML_DOCUMENT_END_EVENT) {
            if (event.type == YAML_STREAM_START_EVENT) continue;
            if (event.type == YAML_DOCUMENT_START_EVENT) continue;
            if (event.type == YAML_MAPPING_START_EVENT) doclevel++; 
            if (event.type == YAML_MAPPING_END_EVENT) doclevel--;
            if (event.type == YAML_SCALAR_EVENT && !strcmp((char *) event.data.scalar.value, tok) && doclevel==proplevel) break;
            if (doclevel < proplevel) {
                PetscPrintf(PETSC_COMM_WORLD,"Warning: mapping %s not found \n", mappingname);
                return 1;
            }
        } 
        //advance the token
        tok = strtok_r(NULL, "/", &savetok); proplevel++;
    }
    
    if (event.type == YAML_DOCUMENT_END_EVENT) {
        PetscPrintf(PETSC_COMM_WORLD,"Warning: mapping %s not found \n", mappingname);
        return 1;
    }

    /* Parsing config file */
    PetscInt mappinglevel = 0;
    while (yaml_parser_parse(&parser, &event) && !(event.type == YAML_MAPPING_END_EVENT && !mappinglevel)) {
        if (event.type == YAML_MAPPING_START_EVENT) mappinglevel++;
        if (event.type == YAML_MAPPING_END_EVENT) mappinglevel--;
        if (event.type == YAML_SCALAR_EVENT && !strcmp((char *) event.data.scalar.value, propname)) {
            yaml_parser_parse(&parser, &event);
            if (event.type == YAML_SEQUENCE_START_EVENT) {
                for (*propsize = 0; yaml_parser_parse(&parser, &event) && event.type != YAML_SEQUENCE_END_EVENT; ) {
                    if (event.type == YAML_SCALAR_EVENT) {
                        strcpy(propval[*propsize], (char *) event.data.scalar.value);
                        (*propsize)++;
                    }
                }
            } else if (event.type == YAML_SCALAR_EVENT) {
                *propsize = 1;
                strcpy(propval[0], (char *) event.data.scalar.value);
            }
            break;
        }
    }    

    /* Cleanup */
    if (event.type == YAML_MAPPING_END_EVENT) {
        PetscPrintf(PETSC_COMM_WORLD,"Warning: property %s not found in mapping %s \n", propname, mappingname);
        return 1;
    }
    yaml_event_delete(&event);
    yaml_parser_delete(&parser);
    return 0;
}

/*
 SetUpConfig - Parse config file
 */
PetscErrorCode SetUpConfig(AppCtx *user)
{
    char           configfile[PETSC_MAX_PATH_LEN] = "";
    unsigned char  *buffer = 0;
    PetscInt       total_cells;
    
    PetscFunctionBeginUser;

    /* Read config file to buffer */
    PetscOptionsGetString(NULL,NULL,"--config",configfile,PETSC_MAX_PATH_LEN,NULL);
    FILE *infile=NULL;
    if (user->worldrank == 0) {
        infile = fopen (configfile, "r");
        if (infile==NULL) {
            printf("Error: config file can't be opened. \n");
            printf("Error: Please make sure the file %s is in the CWD.\n",configfile);
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
    
    char **propval;
    PetscInt propsize, ierr, maxprops = 100;

    propval = malloc(maxprops*sizeof(char *));
    for (PetscInt propctr = 0; propctr < maxprops; propctr++) propval[propctr] = malloc(PETSC_MAX_PATH_LEN);

    /* Parsing config file header */
    {
     /* resolution */
     {
      ierr = GetProperty(propval, &propsize, "header", "grid", buffer, filesize); CHKERRQ(ierr);
      assert(propsize >= 1 && propsize <= 3);
      user->dim = propsize;
      user->resolution = malloc(user->dim*sizeof(PetscInt));
      user->size = malloc(user->dim*sizeof(PetscReal));
      for (PetscInt propctr = 0; propctr < propsize; propctr++) user->resolution[propctr] = atoi(propval[propctr]);
      total_cells = user->dim == 2 ? user->resolution[0]*user->resolution[1] 
                                   : user->resolution[0]*user->resolution[1]*user->resolution[2];
      user->voxelphasemapping = malloc(total_cells*sizeof(PetscInt));
      user->voxelsitemapping = malloc(total_cells*sizeof(PetscInt));
     }
     /* size */
     {
      ierr = GetProperty(propval, &propsize, "header", "size", buffer, filesize); CHKERRQ(ierr);
      assert(propsize == user->dim);
      for (PetscInt propctr = 0; propctr < propsize; propctr++) user->size[propctr] = atof(propval[propctr]);
     }
     /* number of phases */
     {
      ierr = GetProperty(propval, &propsize, "header", "n_phases", buffer, filesize); CHKERRQ(ierr);
      assert(propsize == 1 && atoi(propval[0]) > 0);
      user->npf = atoi(propval[0]);
      user->phasematerialmapping = malloc(user->npf*sizeof(PetscInt));
     }
     /* number of nucleation sites */
     {
      ierr = GetProperty(propval, &propsize, "header", "n_sites", buffer, filesize);
      if (propsize) {
          assert(propsize == 1 && atoi(propval[0]) > 0);
          user->nsites = atoi(propval[0]);
      } else {
          user->nsites = 0;
      }
      user->siteoffset = user->worldrank*(1 + ((user->nsites - 1)/user->worldsize));
      user->nsites_local = (  user->siteoffset + (1 + ((user->nsites - 1)/user->worldsize)) < user->nsites
                            ? user->siteoffset + (1 + ((user->nsites - 1)/user->worldsize)) : user->nsites)
                         - user->siteoffset;
      user->siteactivity_global = malloc(user->nsites*sizeof(char));
      user->siteactivity_local = malloc(user->nsites_local*sizeof(char));
      for (PetscInt site=0;site<user->nsites_local;site++) user->siteactivity_local[site] = 1;
      user->sitenucleusmapping = malloc(user->nsites*sizeof(PetscInt));
      user->sitephasemapping = malloc(user->nsites*sizeof(PetscInt));
     }
     /* material names */
     {
      ierr = GetProperty(propval, &propsize, "header", "materials", buffer, filesize); CHKERRQ(ierr);
      assert(propsize > 0);
      user->nmat = propsize;
      user->materialname = malloc(user->nmat*sizeof(char *));
      for (PetscInt m=0; m<user->nmat; m++) {
          user->materialname[m] = malloc(PETSC_MAX_PATH_LEN);
          strcpy(user->materialname[m],propval[m]);
      }    
      user->material = (MATERIAL *) malloc(user->nmat*sizeof(struct MATERIAL));
     }
     /* nuclei names */
     {
      ierr = GetProperty(propval, &propsize, "header", "nuclei", buffer, filesize);
      assert(propsize >= 0);
      user->nnuclei = propsize;
      user->nucleusname = malloc(user->nnuclei*sizeof(char *));
      for (PetscInt n=0; n<user->nnuclei; n++) {
          user->nucleusname[n] = malloc(PETSC_MAX_PATH_LEN);
          strcpy(user->nucleusname[n],propval[n]);
      }    
      user->nucleus = (NUCLEUS *) malloc(user->nnuclei*sizeof(struct NUCLEUS));
     }
     /* interface names */
     {
      ierr = GetProperty(propval, &propsize, "header", "interfaces", buffer, filesize); CHKERRQ(ierr);
      assert(propsize > 0);
      user->nf = propsize;
      user->interfacename = malloc(user->nf*sizeof(char *));
      for (PetscInt m=0; m<user->nf; m++) {
          user->interfacename[m] = malloc(PETSC_MAX_PATH_LEN);
          strcpy(user->interfacename[m],propval[m]);
      }    
      user->interface = (INTERFACE *) malloc(user->nf*sizeof(struct INTERFACE));
     }
     /* component names */
     {
      ierr = GetProperty(propval, &propsize, "header", "components", buffer, filesize); CHKERRQ(ierr);
      assert(propsize > 0 && propsize <= MAXCP);
      user->ncp = propsize;
      user->componentname = malloc(user->ncp*sizeof(char *));
      for (PetscInt c=0; c<user->ncp; c++) {
          user->componentname[c] = malloc(PETSC_MAX_PATH_LEN);
          strcpy(user->componentname[c],propval[c]);
      }
      user->ndp = user->ncp-1;
     }
     /* output names */
     {
      ierr = GetProperty(propval, &propsize, "header", "outputs", buffer, filesize); CHKERRQ(ierr);
      assert(propsize >= 0);
      user->noutputs = propsize;
      user->outputname = malloc(user->noutputs*sizeof(char *));
      for (PetscInt o=0; o<user->noutputs; o++) {
          user->outputname[o] = malloc(PETSC_MAX_PATH_LEN);
          strcpy(user->outputname[o],propval[o]);
      }
     }
    } 
    
    /* field offsets */
    AS_OFFSET = 0;
    AS_SIZE   = (MAXAP < user->npf ? MAXAP : user->npf) + 1;
    PF_OFFSET = AS_OFFSET+AS_SIZE;
    PF_SIZE   = (MAXAP < user->npf ? MAXAP : user->npf);
    DP_OFFSET = PF_OFFSET+PF_SIZE;
    DP_SIZE   = user->ndp;
    CP_OFFSET = DP_OFFSET+DP_SIZE;
    CP_SIZE   = PF_SIZE*user->ndp;

    /* Parsing config file material */
    {
     MATERIAL *currentmaterial = &user->material[0];
     for (PetscInt m=0; m<user->nmat; m++,currentmaterial++) {
         char materialmapping[PETSC_MAX_PATH_LEN];
         sprintf(materialmapping,"material/%s",user->materialname[m]);
         /* molar volume */
         ierr = GetProperty(propval, &propsize, materialmapping, "molarvolume", buffer, filesize); CHKERRQ(ierr);
         assert(propsize == 1);
         currentmaterial->molarvolume = atof(propval[0]);
         /* initial composition */
         ierr = GetProperty(propval, &propsize, materialmapping, "c0", buffer, filesize); CHKERRQ(ierr);
         assert(propsize == user->ncp);
         currentmaterial->c0 = malloc(user->ncp*sizeof(PetscReal));
         for (PetscInt propctr = 0; propctr < user->ncp; propctr++) currentmaterial->c0[propctr] = atof(propval[propctr]);
         /* chempot explicit kinetic coefficient */
         ierr = GetProperty(propval, &propsize, materialmapping, "chempot_ex_kineticcoeff", buffer, filesize);
         if (propsize) {assert(propsize == 1); currentmaterial->chempot_ex_kineticcoeff = atof(propval[0]);} else {currentmaterial->chempot_ex_kineticcoeff = 0.0;}
         /* chemical energy type */
         {
          ierr = GetProperty(propval, &propsize, materialmapping, "chemicalenergy", buffer, filesize); CHKERRQ(ierr);
          assert(propsize == 1);
          if        (!strcmp(propval[0], "quadratic" )) {
              currentmaterial->chemfe_model = QUADRATIC_CHEMENERGY;
              QUAD *currentquad = &currentmaterial->energy.quad;
              /* solute mobility */
              {
               currentquad->mobilityc = malloc(user->ncp*sizeof(PetscReal));
               ierr = GetProperty(propval, &propsize, materialmapping, "mobilityc", buffer, filesize); CHKERRQ(ierr);
               assert(propsize == user->ncp);
               for (PetscInt propctr = 0; propctr < user->ncp; propctr++) currentquad->mobilityc[propctr] = atof(propval[propctr]);
              } 
              /* equilibrium composition */
              {
               char ceqmapping[PETSC_MAX_PATH_LEN];
               currentquad->ceq = (TSeries *) malloc(user->ncp*sizeof(TSeries));
               TSeries *currentceq = &currentquad->ceq[0];
               for (PetscInt c=0; c<user->ncp; c++, currentceq++) {
                   sprintf(ceqmapping, "%s/ceq/%s",materialmapping,user->componentname[c]);
                   ierr = GetProperty(propval, &propsize, ceqmapping, "t_coefficient", buffer, filesize);
                   currentceq->nTser = propsize;
                   for (PetscInt propctr = 0; propctr < propsize; propctr++) currentceq->coeff[propctr] = atof(propval[propctr]);
                   ierr = GetProperty(propval, &propsize, ceqmapping, "t_exponent", buffer, filesize);
                   assert(propsize == currentceq->nTser);
                   for (PetscInt propctr = 0; propctr < propsize; propctr++) currentceq->exp[propctr] = atoi(propval[propctr]);
                   currentceq->logCoeff = 0.0;
               }
              } 
              /* reference enthalpy */
              {
               char refmapping[PETSC_MAX_PATH_LEN];
               sprintf(refmapping, "%s/ref_enthalpy",materialmapping);
               ierr = GetProperty(propval, &propsize, refmapping, "t_coefficient", buffer, filesize);
               currentquad->ref.nTser = propsize;
               for (PetscInt propctr = 0; propctr < propsize; propctr++) currentquad->ref.coeff[propctr] = atof(propval[propctr]);
               ierr = GetProperty(propval, &propsize, refmapping, "t_exponent", buffer, filesize);
               assert(propsize == currentquad->ref.nTser);
               for (PetscInt propctr = 0; propctr < propsize; propctr++) currentquad->ref.exp[propctr] = atoi(propval[propctr]);
               currentquad->ref.logCoeff = 0.0;
              } 
              /* unary enthalpy */
              {
               char unarymapping[PETSC_MAX_PATH_LEN];
               currentquad->unary = (TSeries *) malloc(user->ncp*sizeof(TSeries));
               TSeries *currentunary = &currentquad->unary[0];
               for (PetscInt c=0; c<user->ncp; c++, currentunary++) {
                   sprintf(unarymapping, "%s/unary_enthalpy/%s",materialmapping,user->componentname[c]);
                   ierr = GetProperty(propval, &propsize, unarymapping, "t_coefficient", buffer, filesize);
                   currentunary->nTser = propsize;
                   for (PetscInt propctr = 0; propctr < propsize; propctr++) currentunary->coeff[propctr] = atof(propval[propctr]);
                   ierr = GetProperty(propval, &propsize, unarymapping, "t_exponent", buffer, filesize);
                   assert(propsize == currentunary->nTser);
                   for (PetscInt propctr = 0; propctr < propsize; propctr++) currentunary->exp[propctr] = atoi(propval[propctr]);
                   currentunary->logCoeff = 0.0;
               }
              } 
              /* binary enthalpy */
              {
               char binarymapping[PETSC_MAX_PATH_LEN];
               currentquad->binary = (TSeries *) malloc(user->ncp*sizeof(TSeries));
               TSeries *currentbinary = &currentquad->binary[0];
               for (PetscInt c=0; c<user->ncp; c++, currentbinary++) {
                   sprintf(binarymapping, "%s/binary_enthalpy/%s",materialmapping,user->componentname[c]);
                   ierr = GetProperty(propval, &propsize, binarymapping, "t_coefficient", buffer, filesize);
                   currentbinary->nTser = propsize;
                   for (PetscInt propctr = 0; propctr < propsize; propctr++) currentbinary->coeff[propctr] = atof(propval[propctr]);
                   ierr = GetProperty(propval, &propsize, binarymapping, "t_exponent", buffer, filesize);
                   assert(propsize == currentbinary->nTser);
                   for (PetscInt propctr = 0; propctr < propsize; propctr++) currentbinary->exp[propctr] = atoi(propval[propctr]);
                   currentbinary->logCoeff = 0.0;
               }
              } 
          } else if (!strcmp(propval[0], "calphaddis")) {
              currentmaterial->chemfe_model = CALPHAD_CHEMENERGY;
              CALPHAD *currentcalphad = &currentmaterial->energy.calphad;
              /* mobility */
              {
               char mobmapping[PETSC_MAX_PATH_LEN];
               currentcalphad->mobilityc = (MOBILITY *) malloc(user->ncp*sizeof(MOBILITY));
               MOBILITY *currentmobility = &currentcalphad->mobilityc[0];
               for (PetscInt c=0; c<user->ncp; c++, currentmobility++) {
                   sprintf(mobmapping, "%s/mobilityc/%s",materialmapping,user->componentname[c]);
                   ierr = GetProperty(propval, &propsize, mobmapping, "mobility0", buffer, filesize);
                   assert(propsize == 1);
                   currentmobility->m0 = atof(propval[0]);
                   currentmobility->unary = (TSeries *) malloc(user->ncp*sizeof(TSeries));
                   currentmobility->binary = (RK *) malloc(user->ncp*(user->ncp-1)/2*sizeof(struct RK));
                   TSeries *currentunary = &currentmobility->unary[0];
                   RK *currentbinary = &currentmobility->binary[0];
                   for (PetscInt cj=0; cj<user->ncp; cj++, currentunary++) {
                       sprintf(mobmapping, "%s/mobilityc/%s/unary_migration/%s",materialmapping,user->componentname[c],user->componentname[cj]);
                       ierr = GetProperty(propval, &propsize, mobmapping, "t_coefficient", buffer, filesize);
                       currentunary->nTser = propsize;
                       for (PetscInt propctr = 0; propctr < propsize; propctr++) currentunary->coeff[propctr] = atof(propval[propctr]);
                       ierr = GetProperty(propval, &propsize, mobmapping, "t_exponent", buffer, filesize);
                       assert(propsize == currentunary->nTser);
                       for (PetscInt propctr = 0; propctr < propsize; propctr++) currentunary->exp[propctr] = atoi(propval[propctr]);
                       currentunary->logCoeff = 0.0;
                       for (PetscInt ci=cj+1; ci<user->ncp; ci++, currentbinary++) {
                           sprintf(mobmapping, "%s/mobilityc/%s/binary_migration/%s/%s",materialmapping,user->componentname[c],user->componentname[cj],user->componentname[ci]);
                           ierr = GetProperty(propval, &propsize, mobmapping, "nrk", buffer, filesize);
                           if (propsize) {currentbinary->n = atoi(propval[0]);} else {currentbinary->n = 0;}
                           currentbinary->enthalpy = (TSeries *) malloc(currentbinary->n*sizeof(TSeries));
                           TSeries *currententhalpy = &currentbinary->enthalpy[0];
                           for (PetscInt nrk=0; nrk < currentbinary->n; nrk++, currententhalpy++) {
                               sprintf(mobmapping, "%s/mobilityc/%s/binary_migration/%s/%s/rk_%d",materialmapping,user->componentname[c],user->componentname[cj],user->componentname[ci],nrk);
                               ierr = GetProperty(propval, &propsize, mobmapping, "t_coefficient", buffer, filesize);
                               currententhalpy->nTser = propsize;
                               for (PetscInt propctr = 0; propctr < propsize; propctr++) currententhalpy->coeff[propctr] = atof(propval[propctr]);
                               ierr = GetProperty(propval, &propsize, mobmapping, "t_exponent", buffer, filesize);
                               assert(propsize == currententhalpy->nTser);
                               for (PetscInt propctr = 0; propctr < propsize; propctr++) currententhalpy->exp[propctr] = atoi(propval[propctr]);
                               currententhalpy->logCoeff = 0.0;
                           }
                       }
                   }
               }    
              }
              /* reference enthalpy */
              {
               char refmapping[PETSC_MAX_PATH_LEN];
               sprintf(refmapping, "%s/ref_enthalpy",materialmapping);
               ierr = GetProperty(propval, &propsize, refmapping, "t_coefficient", buffer, filesize);
               currentcalphad->ref.nTser = propsize;
               for (PetscInt propctr = 0; propctr < propsize; propctr++) currentcalphad->ref.coeff[propctr] = atof(propval[propctr]);
               ierr = GetProperty(propval, &propsize, refmapping, "t_exponent", buffer, filesize);
               assert(currentcalphad->ref.nTser == propsize);
               for (PetscInt propctr = 0; propctr < propsize; propctr++) currentcalphad->ref.exp[propctr] = atoi(propval[propctr]);
               currentcalphad->ref.logCoeff = 0.0;
              }
              /* unary enthalpy */
              {
               char unarymapping[PETSC_MAX_PATH_LEN];
               currentcalphad->unary = (TSeries *) malloc(user->ncp*sizeof(TSeries));
               TSeries *currentunary = &currentcalphad->unary[0];
               for (PetscInt c=0; c<user->ncp; c++, currentunary++) {
                   sprintf(unarymapping, "%s/unary_enthalpy/%s",materialmapping,user->componentname[c]);
                   ierr = GetProperty(propval, &propsize, unarymapping, "t_coefficient", buffer, filesize);
                   currentunary->nTser = propsize;
                   for (PetscInt propctr = 0; propctr < propsize; propctr++) currentunary->coeff[propctr] = atof(propval[propctr]);
                   ierr = GetProperty(propval, &propsize, unarymapping, "t_exponent", buffer, filesize);
                   assert(propsize == currentunary->nTser);
                   for (PetscInt propctr = 0; propctr < propsize; propctr++) currentunary->exp[propctr] = atoi(propval[propctr]);
                   ierr = GetProperty(propval, &propsize, unarymapping, "tlnt_coefficient", buffer, filesize);
                   if (propsize) {currentunary->logCoeff = atof(propval[0]);} else {currentunary->logCoeff = 0.0;}
               }
              }
              /* binary enthalpy */
              {
               char binarymapping[PETSC_MAX_PATH_LEN];
               currentcalphad->binary = (RK *) malloc(user->ncp*(user->ncp-1)*sizeof(struct RK)/2);
               RK *currentbinary = &currentcalphad->binary[0];
               for (PetscInt ck=0; ck<user->ncp; ck++) {
                   for (PetscInt cj=ck+1; cj<user->ncp; cj++,currentbinary++) {
                       sprintf(binarymapping, "%s/binary_enthalpy/%s/%s",materialmapping,user->componentname[ck],user->componentname[cj]);
                       ierr = GetProperty(propval, &propsize, binarymapping, "nrk", buffer, filesize);
                       if (propsize) {currentbinary->n = atoi(propval[0]);} else {currentbinary->n = 0;}
                       currentbinary->enthalpy = (TSeries *) malloc(currentbinary->n*sizeof(TSeries));
                       TSeries *currententhalpy = &currentbinary->enthalpy[0];
                       for (PetscInt nrk=0; nrk < currentbinary->n; nrk++, currententhalpy++) {
                           sprintf(binarymapping, "%s/binary_enthalpy/%s/%s/rk_%d",materialmapping,user->componentname[ck],user->componentname[cj],nrk);
                           ierr = GetProperty(propval, &propsize, binarymapping, "t_coefficient", buffer, filesize);
                           currententhalpy->nTser = propsize;
                           for (PetscInt propctr = 0; propctr < propsize; propctr++) currententhalpy->coeff[propctr] = atof(propval[propctr]);
                           ierr = GetProperty(propval, &propsize, binarymapping, "t_exponent", buffer, filesize);
                           assert(propsize == currententhalpy->nTser);
                           for (PetscInt propctr = 0; propctr < propsize; propctr++) currententhalpy->exp[propctr] = atoi(propval[propctr]);
                           currententhalpy->logCoeff = 0.0;
                       }
                   }
               }
              }  
              /* ternary enthalpy */
              {
               char ternarymapping[PETSC_MAX_PATH_LEN];
               currentcalphad->ternary = (RK *) malloc(user->ncp*(user->ncp-1)*(user->ncp-2)*sizeof(struct RK)/6);
               RK *currentternary = &currentcalphad->ternary[0];
               for (PetscInt ck=0;ck<user->ncp;ck++) {
                   for (PetscInt cj=ck+1;cj<user->ncp;cj++) {
                       for (PetscInt ci=cj+1;ci<user->ncp;ci++, currentternary++) {
                           sprintf(ternarymapping, "%s/ternary_enthalpy/%s/%s/%s",materialmapping,user->componentname[ck],user->componentname[cj],user->componentname[ci]);
                           ierr = GetProperty(propval, &propsize, ternarymapping, "nrk", buffer, filesize);
                           if (propsize) {currentternary->n = atoi(propval[0]);} else {currentternary->n = 0;}
                           currentternary->enthalpy = (TSeries *) malloc(currentternary->n*sizeof(TSeries));
                           currentternary->i = malloc(currentternary->n*sizeof(PetscInt));
                           TSeries *currententhalpy = &currentternary->enthalpy[0];
                           for (PetscInt nrk=0; nrk < currentternary->n; nrk++, currententhalpy++) {
                               sprintf(ternarymapping, "%s/ternary_enthalpy/%s/%s/%s/rk_%d",materialmapping,user->componentname[ck],user->componentname[cj],user->componentname[ci],nrk);
                               ierr = GetProperty(propval, &propsize, ternarymapping, "t_coefficient", buffer, filesize);
                               currententhalpy->nTser = propsize;
                               for (PetscInt propctr = 0; propctr < propsize; propctr++) currententhalpy->coeff[propctr] = atof(propval[propctr]);
                               ierr = GetProperty(propval, &propsize, ternarymapping, "t_exponent", buffer, filesize);
                               assert(propsize == currententhalpy->nTser);
                               for (PetscInt propctr = 0; propctr < propsize; propctr++) currententhalpy->exp[propctr] = atoi(propval[propctr]);
                               ierr = GetProperty(propval, &propsize, ternarymapping, "index", buffer, filesize);
                               assert(propsize == 1);
                               currentternary->i[nrk] = atoi(propval[0]);
                               currententhalpy->logCoeff = 0.0;
                           }
                       }
                   }
               }        
              }  
          } else {
              currentmaterial->chemfe_model = NONE_CHEMENERGY;
          }
         }
     }
    }

    /* Parsing config file nucleation */
    {
     NUCLEUS *currentnucleus = &user->nucleus[0];
     for (PetscInt n=0; n<user->nnuclei; n++,currentnucleus++) {
         char nucleusmapping[PETSC_MAX_PATH_LEN];
         sprintf(nucleusmapping,"nucleus/%s",user->nucleusname[n]);
         /* matrix phases */
         ierr = GetProperty(propval, &propsize, nucleusmapping, "matrix", buffer, filesize);
         assert(propsize > 0); currentnucleus->matrixlist[0] = propsize;
         for (PetscInt propctr = 0; propctr < propsize; propctr++) currentnucleus->matrixlist[propctr+1] = atoi(propval[propctr]);
         ierr = GetProperty(propval, &propsize, nucleusmapping, "nucleation_model", buffer, filesize);
         assert(propsize == 1);
         if (!strcmp(propval[0], "cnt" )) {
             currentnucleus->nuc_model = CNT_NUCLEATION;
             CNT_NUC *currentcntnuc = &currentnucleus->nucleation.cnt;
             /* surface energy */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "gamma", buffer, filesize);
             assert(propsize == 1); currentcntnuc->gamma = atof(propval[0]);
             /* heterogeneous nucleation shape factor */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "shape_factor", buffer, filesize);
             if (propsize) {assert(propsize == 1); currentcntnuc->shapefactor = atof(propval[0]);} 
             else {currentcntnuc->shapefactor = 1.0;}
             /* normalized diffusion coefficient */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "D0", buffer, filesize);
             assert(propsize == 1); currentcntnuc->D0 = atof(propval[0]);
             /* normalized migration energy */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "migration", buffer, filesize);
             assert(propsize == 1); currentcntnuc->migration = atof(propval[0]);
             /* lattice parameter */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "minsize", buffer, filesize);
             if (propsize) {assert(propsize == 1); currentcntnuc->minsize = atof(propval[0]);} 
             else {currentcntnuc->minsize = 0.0;}
             /* atomic volume */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "atomic_volume", buffer, filesize);
             assert(propsize == 1); currentcntnuc->atomicvolume = atof(propval[0]);
             /* normalization length scale */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "length_scale", buffer, filesize);
             assert(propsize == 1); currentcntnuc->lengthscale = atof(propval[0]);
         } else if (!strcmp(propval[0], "constant" )) {
             currentnucleus->nuc_model = CONST_NUCLEATION;
             CONST_NUC *currentconstnuc = &currentnucleus->nucleation.constant;
             /* nucleation rate */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "nucleation_rate", buffer, filesize);
             assert(propsize == 1); currentconstnuc->nucleation_rate = atof(propval[0]);
         } else if (!strcmp(propval[0], "thermal" )) {
             currentnucleus->nuc_model = THERMAL_NUCLEATION;
             THERMAL_NUC *currentthermalnuc = &currentnucleus->nucleation.thermal;
             /* solvus temperature */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "solvus_temperature", buffer, filesize);
             assert(propsize == 1); currentthermalnuc->solvus_temperature = atof(propval[0]);
             /* enthalpy of fusion */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "enthalpy_fusion", buffer, filesize);
             assert(propsize == 1); currentthermalnuc->enthalpy_fusion = atof(propval[0]);
             /* surface energy */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "gamma", buffer, filesize);
             assert(propsize == 1); currentthermalnuc->gamma = atof(propval[0]);
             /* heterogeneous nucleation shape factor */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "shape_factor", buffer, filesize);
             if (propsize) {assert(propsize == 1); currentthermalnuc->shapefactor = atof(propval[0]);} 
             else {currentthermalnuc->shapefactor = 1.0;}
             /* normalized diffusion coefficient */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "D0", buffer, filesize);
             assert(propsize == 1); currentthermalnuc->D0 = atof(propval[0]);
             /* normalized migration energy */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "migration", buffer, filesize);
             assert(propsize == 1); currentthermalnuc->migration = atof(propval[0]);
             /* lattice parameter */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "minsize", buffer, filesize);
             if (propsize) {assert(propsize == 1); currentthermalnuc->minsize = atof(propval[0]);} 
             else {currentthermalnuc->minsize = 0.0;}
             /* atomic volume */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "atomic_volume", buffer, filesize);
             assert(propsize == 1); currentthermalnuc->atomicvolume = atof(propval[0]);
             /* normalization length scale */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "length_scale", buffer, filesize);
             assert(propsize == 1); currentthermalnuc->lengthscale = atof(propval[0]);
         } else {
             currentnucleus->nuc_model = NONE_NUCLEATION;
         }
     }
    }

    /* Parsing config file interface */
    {
     INTERFACE *currentinterface = &user->interface[0];
     for (PetscInt interface=0; interface<user->nf; interface++,currentinterface++) {
         char interfacemapping[PETSC_MAX_PATH_LEN];
         sprintf(interfacemapping,"interface/%s",user->interfacename[interface]);
         /* interface energy */
         ierr = GetProperty(propval, &propsize, interfacemapping, "energy", buffer, filesize); CHKERRQ(ierr);
         assert(propsize == 1);
         currentinterface->energy = atof(propval[0]);
         /* interface mobility */
         {
          char mobmapping[PETSC_MAX_PATH_LEN];
          currentinterface->mobility = (MOBILITY *) malloc(sizeof(MOBILITY));
          sprintf(mobmapping, "%s/mobility",interfacemapping);
          ierr = GetProperty(propval, &propsize, mobmapping, "m0", buffer, filesize); CHKERRQ(ierr);
          assert(propsize == 1); currentinterface->mobility->m0 = atof(propval[0]);
          currentinterface->mobility->unary = (TSeries *) malloc(sizeof(TSeries));
          sprintf(mobmapping, "%s/mobility/activation_energy",interfacemapping);
          ierr = GetProperty(propval, &propsize, mobmapping, "t_coefficient", buffer, filesize);
          currentinterface->mobility->unary->nTser = propsize;
          for (PetscInt propctr = 0; propctr < propsize; propctr++) currentinterface->mobility->unary->coeff[propctr] = atof(propval[propctr]);
          ierr = GetProperty(propval, &propsize, mobmapping, "t_exponent", buffer, filesize);
          assert(propsize == currentinterface->mobility->unary->nTser);
          for (PetscInt propctr = 0; propctr < propsize; propctr++) currentinterface->mobility->unary->exp[propctr] = atoi(propval[propctr]);
          currentinterface->mobility->unary->logCoeff = 0.0;
         }
         /* segregation potential */
         ierr = GetProperty(propval, &propsize, interfacemapping, "potential", buffer, filesize);
         if (propsize) {
             assert(propsize == user->ndp); currentinterface->potential = malloc(user->ndp*sizeof(PetscReal));
             for (PetscInt propctr = 0; propctr < propsize; propctr++) currentinterface->potential[propctr] = atof(propval[propctr]);
         } else {currentinterface->potential = NULL;}
         /* interface solute mobility */
         ierr = GetProperty(propval, &propsize, interfacemapping, "mobilityc", buffer, filesize);
         if (propsize) {
             assert(propsize == user->ncp); currentinterface->mobilityc = malloc(user->ncp*sizeof(PetscReal));
             for (PetscInt propctr = 0; propctr < propsize; propctr++) currentinterface->mobilityc[propctr] = atof(propval[propctr]);
         } else {currentinterface->mobilityc = NULL;}
     }    
    }

    /* Parsing config file solution parameters */
    {
     ierr = GetProperty(propval, &propsize, "solution_parameters", "interpolation", buffer, filesize); CHKERRQ(ierr);
     assert(propsize == 1);
     if        (!strcmp(propval[0], "linear"   )) {
         user->solparams.interpolation = LINEAR_INTERPOLATION;
     } else if (!strcmp(propval[0], "quadratic")) {
         user->solparams.interpolation = QUADRATIC_INTERPOLATION;
     } else if (!strcmp(propval[0], "cubic"    )) {
         user->solparams.interpolation = CUBIC_INTERPOLATION;
     } else {
         user->solparams.interpolation = NONE_INTERPOLATION;
     }
     assert(user->solparams.interpolation != NONE_INTERPOLATION);
     ierr = GetProperty(propval, &propsize, "solution_parameters", "time", buffer, filesize); CHKERRQ(ierr);
     assert(propsize > 1); user->solparams.nloadcases = propsize; user->solparams.currentloadcase = 0; 
     user->solparams.time = malloc(user->solparams.nloadcases*sizeof(PetscReal));
     for (PetscInt propctr = 0; propctr < propsize; propctr++) user->solparams.time[propctr] = atof(propval[propctr]);
     ierr = GetProperty(propval, &propsize, "solution_parameters", "temperature", buffer, filesize); CHKERRQ(ierr);
     assert(propsize == user->solparams.nloadcases);
     user->solparams.temperature = malloc(user->solparams.nloadcases*sizeof(PetscReal));
     for (PetscInt propctr = 0; propctr < propsize; propctr++) user->solparams.temperature[propctr] = atof(propval[propctr]);
     ierr = GetProperty(propval, &propsize, "solution_parameters", "timestep0", buffer, filesize); 
     if (propsize) {assert(propsize == 1); user->solparams.timestep = atof(propval[0]);} 
     else {user->solparams.timestep = 1.0e-18;}
     ierr = GetProperty(propval, &propsize, "solution_parameters", "timestepmin", buffer, filesize);
     if (propsize) {assert(propsize == 1); user->solparams.mintimestep = atof(propval[0]);} 
     else {user->solparams.mintimestep = 1.0e-18;}
     ierr = GetProperty(propval, &propsize, "solution_parameters", "timestepmax", buffer, filesize);
     if (propsize) {assert(propsize == 1); user->solparams.maxtimestep = atof(propval[0]);} 
     else {user->solparams.maxtimestep = user->solparams.time[user->solparams.nloadcases-1];}
     ierr = GetProperty(propval, &propsize, "solution_parameters", "interfacewidth", buffer, filesize); CHKERRQ(ierr);
     assert(propsize == 1); user->solparams.interfacewidth = atof(propval[0]);
     ierr = GetProperty(propval, &propsize, "solution_parameters", "initblocksize", buffer, filesize); CHKERRQ(ierr);
     assert(propsize == user->dim); 
     user->amrparams.initblocksize = malloc(user->dim*sizeof(PetscInt));
     for (PetscInt propctr = 0; propctr < propsize; propctr++) user->amrparams.initblocksize[propctr] = atoi(propval[propctr]);
     ierr = GetProperty(propval, &propsize, "solution_parameters", "initrefine", buffer, filesize); CHKERRQ(ierr);
     assert(propsize == 1); user->amrparams.initrefine = atoi(propval[0]);
     ierr = GetProperty(propval, &propsize, "solution_parameters", "initcoarsen", buffer, filesize); CHKERRQ(ierr);
     assert(propsize == 1); user->amrparams.initcoarsen = atoi(propval[0]);
     ierr = GetProperty(propval, &propsize, "solution_parameters", "maxnrefine", buffer, filesize); CHKERRQ(ierr);
     assert(propsize == 1); user->amrparams.maxnrefine = atoi(propval[0]);
     ierr = GetProperty(propval, &propsize, "solution_parameters", "minnrefine", buffer, filesize); CHKERRQ(ierr);
     assert(propsize == 1); user->amrparams.minnrefine = atoi(propval[0]);
     ierr = GetProperty(propval, &propsize, "solution_parameters", "amrinterval", buffer, filesize); CHKERRQ(ierr);
     assert(propsize == 1); user->amrparams.amrinterval = atoi(propval[0]);
     ierr = GetProperty(propval, &propsize, "solution_parameters", "reltol", buffer, filesize);
     if (propsize) {assert(propsize == 1); user->solparams.reltol = atof(propval[0]);} 
     else {user->solparams.reltol = 1.0e-6;}
     ierr = GetProperty(propval, &propsize, "solution_parameters", "abstol", buffer, filesize);
     if (propsize) {assert(propsize == 1); user->solparams.abstol = atof(propval[0]);} 
     else {user->solparams.abstol = 1.0e-6;}
     ierr = GetProperty(propval, &propsize, "solution_parameters", "random_seed", buffer, filesize);
     if (propsize) {assert(propsize == 1); user->solparams.randomseed = atoi(propval[0]) + user->worldrank;} 
     else {user->solparams.randomseed = time(0) + user->worldrank;}
     PetscPrintf(PETSC_COMM_WORLD,"Using random seed: %d\n",user->solparams.randomseed);
     srand(user->solparams.randomseed);
     ierr = GetProperty(propval, &propsize, "solution_parameters", "outputfreq", buffer, filesize);
     if (propsize) {assert(propsize == 1); user->solparams.outputfreq = atoi(propval[0]);} 
     else {user->solparams.outputfreq = 1;}
     ierr = GetProperty(propval, &propsize, "solution_parameters", "outfile", buffer, filesize);
     if (propsize) {assert(propsize == 1); strcpy(user->solparams.outfile,propval[0]);} 
     else {strcpy(user->solparams.outfile,"output");}
     ierr = GetProperty(propval, &propsize, "solution_parameters", "petscoptions", buffer, filesize);
     if (propsize) {assert(propsize == 1); strcpy(user->solparams.petscoptions,propval[0]);} 
     else {strcpy(user->solparams.petscoptions,"");}
     strcat(user->solparams.petscoptions," -dm_p4est_brick_bounds ");
     for (PetscInt dim=0; dim<user->dim; ++dim) {
         char strval[128];
         sprintf(strval,"0.0,%e",user->size[dim]);
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
    }

    /* Parsing config file mappings */
    {
     char *tok, *savetok;
     ierr = GetProperty(propval, &propsize, "mappings", "phase_material_mapping", buffer, filesize); CHKERRQ(ierr);
     assert(propsize);
     tok = strtok_r(propval[0], "\n", &savetok);
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

     ierr = GetProperty(propval, &propsize, "mappings", "voxel_phase_mapping", buffer, filesize); CHKERRQ(ierr);
     assert(propsize);
     tok = strtok_r(propval[0], "\n", &savetok);
     ctrm=0;
     while (tok != NULL && ctrm < total_cells) {
         // process the line
         if (strstr(tok, "of") != NULL) {
             char stra[128], strb[128];
             PetscInt sof,p;
             sscanf(tok, "%s of %s", stra, strb);
             sof = atoi(stra); p = atoi(strb);
             for (PetscInt j=ctrm;j<ctrm+sof;j++) 
                 user->voxelphasemapping[j] = p-1;;
             ctrm+=sof;
         } else if (strstr(tok, "to") != NULL) {
             char stra[128], strb[128];
             PetscInt p, pfrom, pto;
             sscanf(tok, "%s to %s", stra, strb);
             pfrom = atoi(stra); pto = atoi(strb);
             if (pfrom < pto) {
                 for (p=pfrom;p<=pto;p++) 
                     user->voxelphasemapping[ctrm++] = p-1;
             } else {
                 for (p=pfrom;p>=pto;p--)
                     user->voxelphasemapping[ctrm++] = p-1;
             }   
         }
         else if (atoi(tok)) {
             user->voxelphasemapping[ctrm++] = atoi(tok)-1;
         }
         //advance the token
         tok = strtok_r(NULL, "\n", &savetok);
     }
     assert(ctrm == total_cells && tok == NULL);

     ierr = GetProperty(propval, &propsize, "mappings", "interface_mapping", buffer, filesize); CHKERRQ(ierr);
     assert(propsize);
     user->interfacelist = malloc(user->npf*user->npf*sizeof(PetscInt));
     tok = strtok_r(propval[0], "\n", &savetok);
     ctrm=0;
     while (tok != NULL && ctrm < user->npf*user->npf) {
         // process the line
         if (strstr(tok, "of") != NULL) {
             char stra[128], strb[128];
             PetscInt sof, interface;
             sscanf(tok, "%s of %s", stra, strb);
             sof = atoi(stra); interface = atoi(strb);
             for (PetscInt j=ctrm;j<ctrm+sof;j++) {
                 PetscInt row = j/user->npf, col = j%user->npf;
                 if (col > row) user->interfacelist[j] = (unsigned char) interface-1;
             }
             ctrm+=sof;
         } else if (strstr(tok, "to") != NULL) {
             char stra[128], strb[128];
             PetscInt interface, interfacefrom, interfaceto;
             sscanf(tok, "%s to %s", stra, strb);
             interfacefrom = atoi(stra); interfaceto = atoi(strb);
             if (interfacefrom < interfaceto) {
                 for (interface=interfacefrom;interface<=interfaceto;interface++) {
                     PetscInt row = ctrm/user->npf, col = ctrm%user->npf;
                     if (col > row) user->interfacelist[ctrm] = (unsigned char) interface-1;
                     ctrm++;
                 }    
             } else {
                 for (interface=interfacefrom;interface>=interfaceto;interface--) {
                     PetscInt row = ctrm/user->npf, col = ctrm%user->npf;
                     if (col > row) user->interfacelist[ctrm] = (unsigned char) interface-1;
                     ctrm++;
                 }    
             }   
         } else if (atoi(tok)) {
             PetscInt interface = atoi(tok);
             PetscInt row = ctrm/user->npf, col = ctrm%user->npf;
             if (col > row) user->interfacelist[ctrm] = (unsigned char) interface-1;
             ctrm++;
         }
         //advance the token
         tok = strtok_r(NULL, "\n", &savetok);
     }
     assert(ctrm == user->npf*user->npf);

     if (user->nsites) {
         ierr = GetProperty(propval, &propsize, "mappings", "site_nucleus_mapping", buffer, filesize);
         assert(propsize);
         tok = strtok_r(propval[0], "\n", &savetok);
         ctrm=0;
         while (tok != NULL && ctrm < user->nsites) {
             // process the line
             if (strstr(tok, "of") != NULL) {
                 char stra[128], strb[128];
                 PetscInt sof,m;
                 sscanf(tok, "%s of %s", stra, strb);
                 sof = atoi(stra); m = atoi(strb);
                 assert(m <= user->nnuclei);
                 for (PetscInt j=ctrm;j<ctrm+sof;j++) 
                     user->sitenucleusmapping[j] = m-1;
                 ctrm += sof;
             } else if (strstr(tok, "to") != NULL) {
                 char stra[128], strb[128];
                 PetscInt m, mfrom, mto;
                 sscanf(tok, "%s to %s", stra, strb);
                 mfrom = atoi(stra); mto = atoi(strb);
                 assert(mfrom <= user->nnuclei);
                 assert(mto   <= user->nnuclei);
                 if (mfrom < mto) {
                     for (m=mfrom;m<=mto;m++) 
                         user->sitenucleusmapping[ctrm++] = m-1;
                 } else {
                     for (m=mfrom;m>=mto;m--)
                         user->sitenucleusmapping[ctrm++] = m-1;
                 }   
             } else if (atoi(tok)) {
                 assert(atoi(tok) <= user->nnuclei);
                 user->sitenucleusmapping[ctrm++] = atoi(tok)-1;
             }
             //advance the token
             tok = strtok_r(NULL, "\n", &savetok);
         }
         assert(ctrm == user->nsites && tok == NULL);

         ierr = GetProperty(propval, &propsize, "mappings", "site_phase_mapping", buffer, filesize);
         assert(propsize);
         tok = strtok_r(propval[0], "\n", &savetok);
         ctrm=0;
         while (tok != NULL && ctrm < user->nsites) {
             // process the line
             if (strstr(tok, "of") != NULL) {
                 char stra[128], strb[128];
                 PetscInt sof,m;
                 sscanf(tok, "%s of %s", stra, strb);
                 sof = atoi(stra); m = atoi(strb);
                 assert(m <= user->npf);
                 for (PetscInt j=ctrm;j<ctrm+sof;j++) 
                     user->sitephasemapping[j] = m-1;
                 ctrm += sof;
             } else if (strstr(tok, "to") != NULL) {
                 char stra[128], strb[128];
                 PetscInt m, mfrom, mto;
                 sscanf(tok, "%s to %s", stra, strb);
                 mfrom = atoi(stra); mto = atoi(strb);
                 assert(mfrom <= user->npf);
                 assert(mto   <= user->npf);
                 if (mfrom < mto) {
                     for (m=mfrom;m<=mto;m++) 
                         user->sitephasemapping[ctrm++] = m-1;
                 } else {
                     for (m=mfrom;m>=mto;m--)
                         user->sitephasemapping[ctrm++] = m-1;
                 }   
             } else if (atoi(tok)) {
                 assert(atoi(tok) <= user->npf);
                 user->sitephasemapping[ctrm++] = atoi(tok)-1;
             }
             //advance the token
             tok = strtok_r(NULL, "\n", &savetok);
         }
         assert(ctrm == user->nsites && tok == NULL);

         ierr = GetProperty(propval, &propsize, "mappings", "voxel_site_mapping", buffer, filesize);
         assert(propsize);
         tok = strtok_r(propval[0], "\n", &savetok);
         ctrm=0;
         while (tok != NULL && ctrm < total_cells) {
             // process the line
             if (strstr(tok, "of") != NULL) {
                 char stra[128], strb[128];
                 PetscInt sof,p;
                 sscanf(tok, "%s of %s", stra, strb);
                 sof = atoi(stra); p = atoi(strb);
                 for (PetscInt j=ctrm;j<ctrm+sof;j++) 
                     user->voxelsitemapping[j] = p-1;;
                 ctrm+=sof;
             } else if (strstr(tok, "to") != NULL) {
                 char stra[128], strb[128];
                 PetscInt p, pfrom, pto;
                 sscanf(tok, "%s to %s", stra, strb);
                 pfrom = atoi(stra); pto = atoi(strb);
                 if (pfrom < pto) {
                     for (p=pfrom;p<=pto;p++) 
                         user->voxelsitemapping[ctrm++] = p-1;
                 } else {
                     for (p=pfrom;p>=pto;p--)
                         user->voxelsitemapping[ctrm++] = p-1;
                 }   
             }
             else if (atoi(tok)) {
                 user->voxelsitemapping[ctrm++] = atoi(tok)-1;
             }
             //advance the token
             tok = strtok_r(NULL, "\n", &savetok);
         }
         assert(ctrm == total_cells && tok == NULL);
     }
    } 

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
    PetscScalar       *pcell, *dcell, *ccell;
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
    PetscReal temperature = Interpolate(0.0,user->solparams.temperature,user->solparams.time,user->solparams.currentloadcase);
    for (localcell = 0; localcell < user->nlocalcells; ++localcell) {
        cell = user->localcells[localcell];

        /* get cell state */
        offset = NULL;
        ierr = DMPlexPointGlobalRef(user->da_solution, cell, fdof, &offset); CHKERRQ(ierr);
        F2IFUNC(gslist,&offset[AS_OFFSET]);
        pcell  = &offset[PF_OFFSET];
        dcell  = &offset[DP_OFFSET];
        ccell  = &offset[CP_OFFSET];
        PetscInt phase;
        ierr = DMLabelGetValue(plabel, cell, &phase);
    
        /* set initial conditions */
        for (g=0; g<gslist[0]; g++) {
            currentmaterial = &user->material[user->phasematerialmapping[gslist[g+1]]];
            ChemicalpotentialExplicit(&ccell[g*user->ndp],currentmaterial->c0,temperature,gslist[g+1],user);
            if (gslist[g+1] == phase) {
                pcell[g] = 1.0;
                Chemicalpotential(dcell,currentmaterial->c0,temperature,gslist[g+1],user);
            } else {
                pcell[g] = 0.0;
            }
        }
    }        
    ierr = VecRestoreArray(solution, &fdof);
    PetscFunctionReturn(0);
}

