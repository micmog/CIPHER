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
    while (yaml_parser_parse(&parser, &event) && event.type != YAML_MAPPING_END_EVENT) {
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
    PetscMPIInt    rank;
    
    PetscFunctionBeginUser;

    /* Read config file to buffer */
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    PetscOptionsGetString(NULL,NULL,"--config",configfile,PETSC_MAX_PATH_LEN,NULL);
    FILE *infile=NULL;
    if (rank == 0) {
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
      user->phasevoxelmapping = malloc(user->resolution[2]*user->resolution[1]*user->resolution[0]*sizeof(uint16_t));
     }
     /* size */
     {
      ierr = GetProperty(propval, &propsize, "header", "size", buffer, filesize); CHKERRQ(ierr);
      assert(propsize == user->dim);
      for (PetscInt propctr = 0; propctr < propsize; propctr++) user->size[propctr] = atof(propval[propctr]);
     }
     /* n_phases */
     {
      ierr = GetProperty(propval, &propsize, "header", "n_phases", buffer, filesize); CHKERRQ(ierr);
      assert(propsize == 1 && propval[0] > 0);
      user->npf = atoi(propval[0]);
      user->phasematerialmapping = malloc(user->npf*sizeof(uint16_t));
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
    } 
    
    /* field offsets */
    AS_OFFSET = 0;
    AS_SIZE   = (MAXAP < user->npf ? MAXAP : user->npf) + 1;
    PF_OFFSET = AS_OFFSET+AS_SIZE;
    PF_SIZE   = (MAXAP < user->npf ? MAXAP : user->npf);
    DP_OFFSET = PF_OFFSET+PF_SIZE;
    DP_SIZE   = user->ndp;
    CP_OFFSET = DP_OFFSET+DP_SIZE;
    CP_SIZE   = PF_SIZE*user->ncp;

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
         /* chemical energy type */
         {
          ierr = GetProperty(propval, &propsize, materialmapping, "chemicalenergy", buffer, filesize); CHKERRQ(ierr);
          assert(propsize == 1);
          if        (!strcmp(propval[0], "quadratic" )) {
              currentmaterial->model = QUADRATIC_CHEMENERGY;
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
              currentmaterial->model = CALPHAD_CHEMENERGY;
              CALPHAD *currentcalphad = &currentmaterial->energy.calphad;
              /* stabilisation constant */
              ierr = GetProperty(propval, &propsize, materialmapping, "stabilisation_const", buffer, filesize);
              if (propsize) {assert(propsize == 1); currentcalphad->stabilisation_const = atof(propval[0]);} else {currentcalphad->stabilisation_const = 0.0;}
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
              currentmaterial->model = NONE_CHEMENERGY;
          }
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
         assert(ierr || propsize == user->ncp);
         currentinterface->potential = malloc(user->ncp*sizeof(PetscReal));
         for (PetscInt propctr = 0; propctr < propsize; propctr++) currentinterface->potential[propctr] = atof(propval[propctr]);
         /* interface solute mobility */
         ierr = GetProperty(propval, &propsize, interfacemapping, "mobilityc", buffer, filesize);
         assert(ierr || propsize == user->ncp);
         currentinterface->mobilityc = malloc(user->ncp*sizeof(PetscReal));
         for (PetscInt propctr = 0; propctr < propsize; propctr++) currentinterface->mobilityc[propctr] = atof(propval[propctr]);
     }    
    }

    /* Parsing config file solution parameters */
    {
     ierr = GetProperty(propval, &propsize, "solution_parameters", "interpolation", buffer, filesize); CHKERRQ(ierr);
     assert(propsize == 1);
     if        (!strcmp(propval[0], "linear"   )) {
         user->interpolation = LINEAR_INTERPOLATION;
     } else if (!strcmp(propval[0], "quadratic")) {
         user->interpolation = QUADRATIC_INTERPOLATION;
     } else if (!strcmp(propval[0], "cubic"    )) {
         user->interpolation = CUBIC_INTERPOLATION;
     } else {
         user->interpolation = NONE_INTERPOLATION;
     }
     assert(user->interpolation != NONE_INTERPOLATION);
     ierr = GetProperty(propval, &propsize, "solution_parameters", "finaltime", buffer, filesize); CHKERRQ(ierr);
     assert(propsize == 1); user->solparams.finaltime = atof(propval[0]);
     ierr = GetProperty(propval, &propsize, "solution_parameters", "timestep0", buffer, filesize); CHKERRQ(ierr);
     assert(propsize == 1); user->solparams.timestep = atof(propval[0]);
     ierr = GetProperty(propval, &propsize, "solution_parameters", "timestepmin", buffer, filesize); CHKERRQ(ierr);
     assert(propsize == 1); user->solparams.mintimestep = atof(propval[0]);
     ierr = GetProperty(propval, &propsize, "solution_parameters", "timestepmax", buffer, filesize); CHKERRQ(ierr);
     assert(propsize == 1); user->solparams.maxtimestep = atof(propval[0]);
     ierr = GetProperty(propval, &propsize, "solution_parameters", "interfacewidth", buffer, filesize); CHKERRQ(ierr);
     assert(propsize == 1); user->solparams.interfacewidth = atof(propval[0]);
     ierr = GetProperty(propval, &propsize, "solution_parameters", "temperature", buffer, filesize); CHKERRQ(ierr);
     assert(propsize > 0); user->solparams.n_temperature = propsize;
     user->solparams.temperature_T = malloc(user->solparams.n_temperature*sizeof(PetscReal));
     for (PetscInt propctr = 0; propctr < propsize; propctr++) user->solparams.temperature_T[propctr] = atof(propval[propctr]);
     ierr = GetProperty(propval, &propsize, "solution_parameters", "time", buffer, filesize); CHKERRQ(ierr);
     assert(propsize == user->solparams.n_temperature);
     user->solparams.temperature_t = malloc(user->solparams.n_temperature*sizeof(PetscReal));
     for (PetscInt propctr = 0; propctr < propsize; propctr++) user->solparams.temperature_t[propctr] = atof(propval[propctr]);
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
     ierr = GetProperty(propval, &propsize, "solution_parameters", "reltol", buffer, filesize); CHKERRQ(ierr);
     assert(propsize == 1); user->solparams.reltol = atof(propval[0]);
     ierr = GetProperty(propval, &propsize, "solution_parameters", "abstol", buffer, filesize); CHKERRQ(ierr);
     assert(propsize == 1); user->solparams.abstol = atof(propval[0]);
     ierr = GetProperty(propval, &propsize, "solution_parameters", "outputfreq", buffer, filesize); CHKERRQ(ierr);
     assert(propsize == 1); user->solparams.outputfreq = atoi(propval[0]);
     ierr = GetProperty(propval, &propsize, "solution_parameters", "outfile", buffer, filesize); CHKERRQ(ierr);
     assert(propsize == 1); strcpy(user->solparams.outfile,propval[0]);
     ierr = GetProperty(propval, &propsize, "solution_parameters", "petscoptions", buffer, filesize); CHKERRQ(ierr);
     assert(propsize == 1); strcpy(user->solparams.petscoptions,propval[0]);
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
     tok = strtok_r(propval[0], "\n", &savetok);
     PetscInt ctrv=0, total_cells = user->dim == 2 ? user->resolution[0]*user->resolution[1] 
                                                   : user->resolution[0]*user->resolution[1]*user->resolution[2];
     while (tok != NULL && ctrv < total_cells) {
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
     assert(ctrv == total_cells && tok == NULL);

     ierr = GetProperty(propval, &propsize, "mappings", "interface_mapping", buffer, filesize); CHKERRQ(ierr);
     user->interfacelist = malloc(user->npf*user->npf*sizeof(uint16_t));
     tok = strtok_r(propval[0], "\n", &savetok);
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
    PetscInt          localcell, cell, face, phase, g, c;
    PetscInt          conesize, supp, nsupp;
    const PetscInt    *cone, *scells;
    DMLabel           plabel = NULL;
    uint16_t          gslist[AS_SIZE], nalist[AS_SIZE];
    PetscScalar       *pcell, *dcell, *ccell, chempot_im[DP_SIZE], chempot_ex[DP_SIZE];
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
    PetscReal temperature = Interpolate(0.0,user->solparams.temperature_T,user->solparams.temperature_t,user->solparams.n_temperature);
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
            memcpy(&ccell[g*user->ncp],currentmaterial->c0,user->ncp*sizeof(PetscReal));
            if (gslist[g+1] == phase) {
                pcell[g] = 1.0;
                Chemicalpotential(dcell,&ccell[g*user->ncp],temperature,gslist[g+1],user);
            } else {
                pcell[g] = 0.0;
            }
        }
    }        
    ierr = VecRestoreArray(solution, &fdof);
    PetscFunctionReturn(0);
}

