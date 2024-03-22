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
    while (yaml_parser_parse(&parser, &event)) {
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
        if (event.type == YAML_MAPPING_END_EVENT && !mappinglevel) break;
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
    unsigned char  *buffer = 0, *mappingbuffer = 0;
    PetscInt       total_cells, max_phase_neighbours = MAXAP;
    
    PetscFunctionBeginUser;

    /* Read config file to buffer */
    PetscOptionsGetString(NULL,NULL,"--config",configfile,PETSC_MAX_PATH_LEN,NULL);
    FILE *infile=NULL;
	infile = fopen (configfile, "r");
	if (infile==NULL) {
		printf("Error: config file can't be opened. \n");
		printf("Error: Please make sure the file %s is in the CWD.\n",configfile);
		return 1;
	}
    size_t filesize;
    if (infile) {
        fseek(infile, 0, SEEK_END);
        filesize = ftell(infile);
        fseek(infile, 0, SEEK_SET);
    }    
    buffer = malloc(filesize);
    if (infile) {
        fread(buffer, 1, filesize, infile);
        fclose(infile);
    }
    if (user->worldrank == 0) {
        printf("Config file size %lu \n", filesize);
    }
    
    char **propval;
    PetscInt propsize, ierr, maxprops = 100;

    propval = malloc(maxprops*sizeof(char *));
    for (PetscInt propctr = 0; propctr < maxprops; propctr++) propval[propctr] = malloc(filesize);

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
      user->interfacelist  = malloc(user->npf * user->npf * sizeof(PetscInt));
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
      PetscReal sitesperproc = ((PetscReal) user->nsites)/((PetscReal) user->worldsize);
      PetscInt sitesperprocctr[user->worldsize];
      memset(sitesperprocctr,0,user->worldsize*sizeof(PetscInt));
      for (PetscInt site = 0; site<user->nsites; site++) {
          sitesperprocctr[(PetscInt) (((PetscReal) site)/sitesperproc)]++;
      }
      user->nsites_local = sitesperprocctr[user->worldrank];
      user->siteoffset = malloc(user->worldsize*sizeof(PetscInt));
      memset(user->siteoffset,0,user->worldsize*sizeof(PetscInt));
      for (PetscInt rank = 0; rank<user->worldsize-1; rank++) {
          user->siteoffset[rank+1] = user->siteoffset[rank] + sitesperprocctr[rank];
      }
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
     /* boundary names */
     {
      ierr = GetProperty(propval, &propsize, "header", "boundaries", buffer, filesize);
      assert(propsize <= 6);
      user->nbcs = propsize;
      user->bcname = malloc(user->nbcs*sizeof(char *));
      BOUNDARYCONDITIONS *currentbc = &user->bcs[0];
      for (PetscInt bc=0; bc<user->nbcs; bc++, currentbc++) {
          user->bcname[bc] = malloc(PETSC_MAX_PATH_LEN);
          strcpy(user->bcname[bc],propval[bc]);
      }
      user->bcs = (BOUNDARYCONDITIONS *) malloc(user->nbcs*sizeof(struct BOUNDARYCONDITIONS));
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
    
    /* Parsing config file material */
    MAXSITES = 0;
    {
     MATERIAL *currentmaterial = &user->material[0];
     for (PetscInt m=0; m<user->nmat; m++,currentmaterial++) {
         char materialmapping[PETSC_MAX_PATH_LEN];
         sprintf(materialmapping,"material/%s",user->materialname[m]);
         /* molar volume */
         ierr = GetProperty(propval, &propsize, materialmapping, "molarvolume", buffer, filesize); CHKERRQ(ierr);
         assert(propsize == 1);
         currentmaterial->molarvolume = atof(propval[0]);
         /* chempot explicit kinetic coefficient */
         ierr = GetProperty(propval, &propsize, materialmapping, "chempot_ex_kineticcoeff", buffer, filesize);
         if (propsize) {assert(propsize == 1); currentmaterial->chempot_ex_kineticcoeff = atof(propval[0]);} else {currentmaterial->chempot_ex_kineticcoeff = 0.0;}
         /* sources */
         {
          char sourcemapping[PETSC_MAX_PATH_LEN];
          sprintf(sourcemapping,"%s/sources",materialmapping);
          ierr = GetProperty(propval, &propsize, sourcemapping, "activesources", buffer, filesize);
          currentmaterial->sources.nsources = propsize;
          currentmaterial->sources.sourcename = malloc(currentmaterial->sources.nsources*sizeof(char *));
          for (PetscInt s=0; s<currentmaterial->sources.nsources; s++) {
              currentmaterial->sources.sourcename[s] = malloc(PETSC_MAX_PATH_LEN);
              strcpy(currentmaterial->sources.sourcename[s],propval[s]);
          }    
          currentmaterial->sources.source = (SOURCE *) malloc(currentmaterial->sources.nsources*sizeof(struct SOURCE));
          SOURCE *currentsource = &currentmaterial->sources.source[0];
          for (PetscInt s=0; s<currentmaterial->sources.nsources; s++,currentsource++) {
              char activesourcemapping[PETSC_MAX_PATH_LEN];
              sprintf(activesourcemapping,"%s/%s",sourcemapping,currentmaterial->sources.sourcename[s]);
              ierr = GetProperty(propval, &propsize, activesourcemapping, "source", buffer, filesize); CHKERRQ(ierr);
              assert(propsize == 1);
              /* equilibration source model */
              if        (!strcmp(propval[0], "sink"   )) {
                  currentsource->source_model = SINK_SOURCE;
                  SOURCE_SINK *currentsinksource = &currentsource->source.sink;
                  /* equilibration rate */
                  ierr = GetProperty(propval, &propsize, activesourcemapping, "rate", buffer, filesize);
                  currentsinksource->rate = malloc(user->ncp*sizeof(PetscReal));
                  if (propsize) {
                      assert(propsize == user->ncp);
                      for (PetscInt propctr = 0; propctr < propsize; propctr++) {
                          currentsinksource->rate[propctr] = atof(propval[propctr]);
                      }    
                  } else {
                      memset(currentsinksource->rate,0,user->ncp*sizeof(PetscReal));
                  }
                  /* equilibrium composition */
                  ierr = GetProperty(propval, &propsize, activesourcemapping, "ceq", buffer, filesize);
                  currentsinksource->ceq = malloc(user->ncp*sizeof(PetscReal));
                  if (propsize) {
                      assert(propsize == user->ncp);
                      for (PetscInt propctr = 0; propctr < propsize; propctr++) {
                          currentsinksource->ceq[propctr] = atof(propval[propctr]);
                      }    
                  } else {
                      memset(currentsinksource->ceq,0,user->ncp*sizeof(PetscReal));
                  }
              /* random fluctuation source model */
              } else if (!strcmp(propval[0], "random" )) {
                  currentsource->source_model = RND_SOURCE;
                  SOURCE_RND *currentrndsource = &currentsource->source.rnd;
                  /* fluctuation amplitude */
                  ierr = GetProperty(propval, &propsize, activesourcemapping, "fluctuation_amplitute", buffer, filesize);
                  currentrndsource->fluctuation_amplitute = malloc(user->ncp*sizeof(PetscReal));
                  if (propsize) {
                      assert(propsize == user->ncp);
                      for (PetscInt propctr = 0; propctr < propsize; propctr++) {
                          currentrndsource->fluctuation_amplitute[propctr] = atof(propval[propctr]);
                      }    
                  } else {
                      memset(currentrndsource->fluctuation_amplitute,0,user->ncp*sizeof(PetscReal));
                  } 
              }
          }
         } 
         /* chemical energy type */
         {
          ierr = GetProperty(propval, &propsize, materialmapping, "chemicalenergy", buffer, filesize); CHKERRQ(ierr);
          assert(propsize == 1);
          if        (!strcmp(propval[0], "none" )) {
              currentmaterial->chemfe_model = NONE_CHEMENERGY;
              currentmaterial->nsites = 0;
              CHEMNONE *currentnone = &currentmaterial->energy.none;
              /* reference enthalpy */
              {
               char refmapping[PETSC_MAX_PATH_LEN];
               sprintf(refmapping, "%s/ref_enthalpy",materialmapping);
               ierr = GetProperty(propval, &propsize, refmapping, "t_coefficient", buffer, filesize);
               currentnone->ref.nTser = propsize;
               for (PetscInt propctr = 0; propctr < propsize; propctr++) currentnone->ref.coeff[propctr] = atof(propval[propctr]);
               ierr = GetProperty(propval, &propsize, refmapping, "t_exponent", buffer, filesize);
               assert(propsize == currentnone->ref.nTser);
               for (PetscInt propctr = 0; propctr < propsize; propctr++) currentnone->ref.exp[propctr] = atoi(propval[propctr]);
               currentnone->ref.logCoeff = 0.0;
              } 
          } else if (!strcmp(propval[0], "quadratic")) {
              currentmaterial->chemfe_model = QUADRATIC_CHEMENERGY;
              currentmaterial->nsites = 1;
              QUAD *currentquad = &currentmaterial->energy.quad;
              /* initial composition */
              {
               ierr = GetProperty(propval, &propsize, materialmapping, "c0", buffer, filesize); CHKERRQ(ierr);
               assert(propsize == currentmaterial->nsites*user->ncp);
               currentmaterial->c0 = malloc(user->ncp*sizeof(PetscReal));
               for (PetscInt propctr = 0; propctr < user->ncp; propctr++) currentmaterial->c0[propctr] = atof(propval[propctr]);
              }
              /* mobility */
              {
               char mobmapping[PETSC_MAX_PATH_LEN];
               currentquad->mobilityc = (TACTIVATIONPROP *) malloc(user->ncp*sizeof(TACTIVATIONPROP));
               TACTIVATIONPROP *currentmobility = &currentquad->mobilityc[0];
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
              currentmaterial->stochiometry = malloc(currentmaterial->nsites*sizeof(PetscReal));
              currentmaterial->stochiometry[0] = 1.0;
          } else if (!strcmp(propval[0], "calphaddis")) {
              currentmaterial->chemfe_model = CALPHADDIS_CHEMENERGY;
              currentmaterial->nsites = 1;
              CALPHADDIS *currentcalphaddis = &currentmaterial->energy.calphaddis;
              /* initial composition */
              {
               ierr = GetProperty(propval, &propsize, materialmapping, "c0", buffer, filesize); CHKERRQ(ierr);
               assert(propsize == currentmaterial->nsites*user->ncp);
               currentmaterial->c0 = malloc(user->ncp*sizeof(PetscReal));
               for (PetscInt propctr = 0; propctr < user->ncp; propctr++) currentmaterial->c0[propctr] = atof(propval[propctr]);
              }
              /* mobility */
              {
               char mobmapping[PETSC_MAX_PATH_LEN];
               currentcalphaddis->mobilityc = (TACTIVATIONPROP *) malloc(user->ncp*sizeof(TACTIVATIONPROP));
               TACTIVATIONPROP *currentmobility = &currentcalphaddis->mobilityc[0];
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
               currentcalphaddis->ref.nTser = propsize;
               for (PetscInt propctr = 0; propctr < propsize; propctr++) currentcalphaddis->ref.coeff[propctr] = atof(propval[propctr]);
               ierr = GetProperty(propval, &propsize, refmapping, "t_exponent", buffer, filesize);
               assert(currentcalphaddis->ref.nTser == propsize);
               for (PetscInt propctr = 0; propctr < propsize; propctr++) currentcalphaddis->ref.exp[propctr] = atoi(propval[propctr]);
               currentcalphaddis->ref.logCoeff = 0.0;
              }
              /* unary enthalpy */
              {
               char unarymapping[PETSC_MAX_PATH_LEN];
               currentcalphaddis->unary = (TSeries *) malloc(user->ncp*sizeof(TSeries));
               TSeries *currentunary = &currentcalphaddis->unary[0];
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
               currentcalphaddis->binary = (RK *) malloc(user->ncp*(user->ncp-1)*sizeof(struct RK)/2);
               RK *currentbinary = &currentcalphaddis->binary[0];
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
               currentcalphaddis->ternary = (RK *) malloc(user->ncp*(user->ncp-1)*(user->ncp-2)*sizeof(struct RK)/6);
               RK *currentternary = &currentcalphaddis->ternary[0];
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
              currentmaterial->stochiometry = malloc(currentmaterial->nsites*sizeof(PetscReal));
              currentmaterial->stochiometry[0] = 1.0;
          } else if (!strcmp(propval[0], "calphad2sl")) {
              currentmaterial->chemfe_model = CALPHAD2SL_CHEMENERGY;
              currentmaterial->nsites = 2;
              CALPHAD2SL *currentcalphad2sl = &currentmaterial->energy.calphad2sl;
              /* initial composition */
              {
               currentmaterial->c0 = malloc(currentmaterial->nsites*user->ncp*sizeof(PetscReal));
               ierr = GetProperty(propval, &propsize, materialmapping, "c0_p", buffer, filesize); CHKERRQ(ierr);
               assert(propsize == user->ncp);
               for (PetscInt propctr = 0; propctr < user->ncp; propctr++) currentmaterial->c0[propctr] = atof(propval[propctr]);
               ierr = GetProperty(propval, &propsize, materialmapping, "c0_q", buffer, filesize); CHKERRQ(ierr);
               assert(propsize == user->ncp);
               for (PetscInt propctr = 0; propctr < user->ncp; propctr++) currentmaterial->c0[user->ncp+propctr] = atof(propval[propctr]);
              }
              /* stochiometry */
              {
               ierr = GetProperty(propval, &propsize, materialmapping, "stochiometry", buffer, filesize); CHKERRQ(ierr);
               assert(propsize == 2);
               currentcalphad2sl->p = atof(propval[0]); currentcalphad2sl->q = atof(propval[1]);
               assert(fabs(currentcalphad2sl->p + currentcalphad2sl->q - 1.0) < TOL);
              } 
              /* mobility */
              {
               char mobmapping[PETSC_MAX_PATH_LEN];
               currentcalphad2sl->mobilityc = (TACTIVATIONPROP *) malloc(user->ncp*sizeof(TACTIVATIONPROP));
               TACTIVATIONPROP *currentmobility = &currentcalphad2sl->mobilityc[0];
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
              /* unary enthalpy */
              {
               char unarymapping[PETSC_MAX_PATH_LEN];
               currentcalphad2sl->unary = (TSeries *) malloc(user->ncp*user->ncp*sizeof(TSeries));
               TSeries *currentunary = &currentcalphad2sl->unary[0];
               for (PetscInt cp=0; cp<user->ncp; cp++) {
                   for (PetscInt cq=0; cq<user->ncp; cq++, currentunary++) {
                       sprintf(unarymapping, "%s/unary_enthalpy/(%s:%s)",materialmapping,user->componentname[cp],user->componentname[cq]);
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
              }
              /* binary enthalpy */
              {
               char binarymapping[PETSC_MAX_PATH_LEN];
               RK *currentbinary;
               currentcalphad2sl->binaryp = (RK *) malloc(user->ncp*user->ncp*(user->ncp-1)*sizeof(struct RK)/2);
               currentbinary = &currentcalphad2sl->binaryp[0];
               for (PetscInt cp_i=0; cp_i<user->ncp; cp_i++) {
                   for (PetscInt cp_j=cp_i+1; cp_j<user->ncp; cp_j++) {
                       for (PetscInt cq=0; cq<user->ncp; cq++, currentbinary++) {
                           sprintf(binarymapping, "%s/binary_enthalpy/(%s,%s:%s)",materialmapping,user->componentname[cp_i],
                                                                                                  user->componentname[cp_j],user->componentname[cq]);
                           ierr = GetProperty(propval, &propsize, binarymapping, "nrk", buffer, filesize);
                           if (propsize) {currentbinary->n = atoi(propval[0]);} else {currentbinary->n = 0;}
                           currentbinary->enthalpy = (TSeries *) malloc(currentbinary->n*sizeof(TSeries));
                           TSeries *currententhalpy = &currentbinary->enthalpy[0];
                           for (PetscInt nrk=0; nrk < currentbinary->n; nrk++, currententhalpy++) {
                               sprintf(binarymapping, "%s/binary_enthalpy/(%s,%s:%s)/rk_%d",materialmapping,user->componentname[cp_i],
                                                                                                            user->componentname[cp_j],user->componentname[cq],nrk);
                               ierr = GetProperty(propval, &propsize, binarymapping, "t_coefficient", buffer, filesize);
                               currententhalpy->nTser = propsize;
                               for (PetscInt propctr = 0; propctr < propsize; propctr++) currententhalpy->coeff[propctr] = atof(propval[propctr]);
                               ierr = GetProperty(propval, &propsize, binarymapping, "t_exponent", buffer, filesize);
                               assert(propsize == currententhalpy->nTser);
                               for (PetscInt propctr = 0; propctr < propsize; propctr++) currententhalpy->exp[propctr] = atoi(propval[propctr]);
                           }
                       }
                   }
               }        
               currentcalphad2sl->binaryq = (RK *) malloc(user->ncp*user->ncp*(user->ncp-1)*sizeof(struct RK)/2);
               currentbinary = &currentcalphad2sl->binaryq[0];
               for (PetscInt cp=0; cp<user->ncp; cp++) {
                   for (PetscInt cq_i=0; cq_i<user->ncp; cq_i++) {
                       for (PetscInt cq_j=cq_i+1; cq_j<user->ncp; cq_j++, currentbinary++) {
                           sprintf(binarymapping, "%s/binary_enthalpy/(%s:%s,%s)",materialmapping,user->componentname[cp],user->componentname[cq_i],
                                                                                                                          user->componentname[cq_j]);
                           ierr = GetProperty(propval, &propsize, binarymapping, "nrk", buffer, filesize);
                           if (propsize) {currentbinary->n = atoi(propval[0]);} else {currentbinary->n = 0;}
                           currentbinary->enthalpy = (TSeries *) malloc(currentbinary->n*sizeof(TSeries));
                           TSeries *currententhalpy = &currentbinary->enthalpy[0];
                           for (PetscInt nrk=0; nrk < currentbinary->n; nrk++, currententhalpy++) {
                               sprintf(binarymapping, "%s/binary_enthalpy/(%s:%s,%s)/rk_%d",materialmapping,user->componentname[cp],user->componentname[cq_i],
                                                                                                                                    user->componentname[cq_j],nrk);
                               ierr = GetProperty(propval, &propsize, binarymapping, "t_coefficient", buffer, filesize);
                               currententhalpy->nTser = propsize;
                               for (PetscInt propctr = 0; propctr < propsize; propctr++) currententhalpy->coeff[propctr] = atof(propval[propctr]);
                               ierr = GetProperty(propval, &propsize, binarymapping, "t_exponent", buffer, filesize);
                               assert(propsize == currententhalpy->nTser);
                               for (PetscInt propctr = 0; propctr < propsize; propctr++) currententhalpy->exp[propctr] = atoi(propval[propctr]);
                           }
                       }
                   }
               }        
              }  
              /* ternary enthalpy */
              {
               char ternarymapping[PETSC_MAX_PATH_LEN];
               RK *currentternary;
               currentcalphad2sl->ternaryp = (RK *) malloc(user->ncp*user->ncp*(user->ncp-1)*(user->ncp-2)*sizeof(struct RK)/6);
               currentternary = &currentcalphad2sl->ternaryp[0];
               for (PetscInt cp_i=0; cp_i<user->ncp; cp_i++) {
                   for (PetscInt cp_j=cp_i+1; cp_j<user->ncp; cp_j++) {
                       for (PetscInt cp_k=cp_j+1; cp_k<user->ncp; cp_k++) {
                           for (PetscInt cq=0; cq<user->ncp; cq++, currentternary++) {
                               sprintf(ternarymapping, "%s/ternary_enthalpy/(%s,%s,%s:%s)",materialmapping,user->componentname[cp_i],
                                                                                                           user->componentname[cp_j],
                                                                                                           user->componentname[cp_k],user->componentname[cq]);
                               ierr = GetProperty(propval, &propsize, ternarymapping, "nrk", buffer, filesize);
                               if (propsize) {currentternary->n = atoi(propval[0]);} else {currentternary->n = 0;}
                               currentternary->enthalpy = (TSeries *) malloc(currentternary->n*sizeof(TSeries));
                               TSeries *currententhalpy = &currentternary->enthalpy[0];
                               for (PetscInt nrk=0; nrk < currentternary->n; nrk++, currententhalpy++) {
                                   sprintf(ternarymapping, "%s/ternary_enthalpy/(%s,%s,%s:%s)/rk_%d",materialmapping,user->componentname[cp_i],
                                                                                                                     user->componentname[cp_j],
                                                                                                                     user->componentname[cp_k],user->componentname[cq],nrk);
                                   ierr = GetProperty(propval, &propsize, ternarymapping, "t_coefficient", buffer, filesize);
                                   currententhalpy->nTser = propsize;
                                   for (PetscInt propctr = 0; propctr < propsize; propctr++) currententhalpy->coeff[propctr] = atof(propval[propctr]);
                                   ierr = GetProperty(propval, &propsize, ternarymapping, "t_exponent", buffer, filesize);
                                   assert(propsize == currententhalpy->nTser);
                                   for (PetscInt propctr = 0; propctr < propsize; propctr++) currententhalpy->exp[propctr] = atoi(propval[propctr]);
                               }
                           }
                       }
                   }
               }        
               currentcalphad2sl->ternaryq = (RK *) malloc(user->ncp*user->ncp*(user->ncp-1)*(user->ncp-2)*sizeof(struct RK)/6);
               currentternary = &currentcalphad2sl->ternaryq[0];
               for (PetscInt cp=0; cp<user->ncp; cp++) {
                   for (PetscInt cq_i=0; cq_i<user->ncp; cq_i++) {
                       for (PetscInt cq_j=cq_i+1; cq_j<user->ncp; cq_j++) {
                           for (PetscInt cq_k=cq_j+1; cq_k<user->ncp; cq_k++, currentternary++) {
                               sprintf(ternarymapping, "%s/ternary_enthalpy/(%s:%s,%s,%s)",materialmapping,user->componentname[cp],user->componentname[cq_i],
                                                                                                                                   user->componentname[cq_j],
                                                                                                                                   user->componentname[cq_k]);
                               ierr = GetProperty(propval, &propsize, ternarymapping, "nrk", buffer, filesize);
                               if (propsize) {currentternary->n = atoi(propval[0]);} else {currentternary->n = 0;}
                               currentternary->enthalpy = (TSeries *) malloc(currentternary->n*sizeof(TSeries));
                               TSeries *currententhalpy = &currentternary->enthalpy[0];
                               for (PetscInt nrk=0; nrk < currentternary->n; nrk++, currententhalpy++) {
                                   sprintf(ternarymapping, "%s/ternary_enthalpy/(%s:%s,%s,%s)/rk_%d",materialmapping,user->componentname[cp],user->componentname[cq_i],
                                                                                                                                             user->componentname[cq_j],
                                                                                                                                             user->componentname[cq_k],nrk);
                                   ierr = GetProperty(propval, &propsize, ternarymapping, "t_coefficient", buffer, filesize);
                                   currententhalpy->nTser = propsize;
                                   for (PetscInt propctr = 0; propctr < propsize; propctr++) currententhalpy->coeff[propctr] = atof(propval[propctr]);
                                   ierr = GetProperty(propval, &propsize, ternarymapping, "t_exponent", buffer, filesize);
                                   assert(propsize == currententhalpy->nTser);
                                   for (PetscInt propctr = 0; propctr < propsize; propctr++) currententhalpy->exp[propctr] = atoi(propval[propctr]);
                               }
                           }
                       }
                   }
               }        
              }  
              currentmaterial->stochiometry = malloc(currentmaterial->nsites*sizeof(PetscReal));
              currentmaterial->stochiometry[0] = currentcalphad2sl->p;
              currentmaterial->stochiometry[1] = currentcalphad2sl->q;
          } else {
              currentmaterial->chemfe_model = NONE_CHEMENERGY;
              currentmaterial->nsites = 0;
          }
         }
         MAXSITES = MAXSITES > currentmaterial->nsites
                  ? MAXSITES : currentmaterial->nsites;

         /* thermal transport */
         {
          char thermalmapping[PETSC_MAX_PATH_LEN];
          THERMAL *currentthermal = &currentmaterial->thermal;
          /* initial temperature */
          {
           ierr = GetProperty(propval, &propsize, materialmapping, "temperature0", buffer, filesize);
           assert(propsize == 1);
           currentthermal->temperature0 = atof(propval[0]);
          }
          /* thermal specific heat */
          {
           sprintf(thermalmapping, "%s/specific_heat",materialmapping);
           ierr = GetProperty(propval, &propsize, thermalmapping, "t_coefficient", buffer, filesize);
           currentthermal->specific_heat.nTser = propsize;
           for (PetscInt propctr = 0; propctr < propsize; propctr++) currentthermal->specific_heat.coeff[propctr] = atof(propval[propctr]);
           ierr = GetProperty(propval, &propsize, thermalmapping, "t_exponent", buffer, filesize);
           assert(propsize == currentthermal->specific_heat.nTser);
           for (PetscInt propctr = 0; propctr < propsize; propctr++) currentthermal->specific_heat.exp[propctr] = atoi(propval[propctr]);
           currentthermal->specific_heat.logCoeff = 0.0;
          }
          /* thermal latent heat */
          {
           sprintf(thermalmapping, "%s/latent_heat",materialmapping);
           ierr = GetProperty(propval, &propsize, thermalmapping, "t_coefficient", buffer, filesize);
           currentthermal->latent_heat.nTser = propsize;
           for (PetscInt propctr = 0; propctr < propsize; propctr++) currentthermal->latent_heat.coeff[propctr] = atof(propval[propctr]);
           ierr = GetProperty(propval, &propsize, thermalmapping, "t_exponent", buffer, filesize);
           assert(propsize == currentthermal->latent_heat.nTser);
           for (PetscInt propctr = 0; propctr < propsize; propctr++) currentthermal->latent_heat.exp[propctr] = atoi(propval[propctr]);
           currentthermal->latent_heat.logCoeff = 0.0;
          }
          /* thermal conductivity */
          {
           sprintf(thermalmapping, "%s/tconductivity",materialmapping);
           ierr = GetProperty(propval, &propsize, thermalmapping, "t_coefficient", buffer, filesize);
           currentthermal->tconductivity.nTser = propsize;
           for (PetscInt propctr = 0; propctr < propsize; propctr++) currentthermal->tconductivity.coeff[propctr] = atof(propval[propctr]);
           ierr = GetProperty(propval, &propsize, thermalmapping, "t_exponent", buffer, filesize);
           assert(propsize == currentthermal->tconductivity.nTser);
           for (PetscInt propctr = 0; propctr < propsize; propctr++) currentthermal->tconductivity.exp[propctr] = atoi(propval[propctr]);
           currentthermal->tconductivity.logCoeff = 0.0;
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
         assert(propsize > 0); 
         currentnucleus->matrixlist = malloc((propsize+1)*sizeof(uint16_t));
         currentnucleus->matrixlist[0] = propsize;
         for (PetscInt propctr = 0; propctr < propsize; propctr++) currentnucleus->matrixlist[propctr+1] = atoi(propval[propctr]);
         /* active solute components */
         ierr = GetProperty(propval, &propsize, nucleusmapping, "active_solutes", buffer, filesize);
         if (propsize) {
             currentnucleus->activesolutes = malloc(user->ncp*sizeof(char));
             memset(currentnucleus->activesolutes,0,user->ncp*sizeof(char));
             for (PetscInt propctr = 0; propctr < propsize; propctr++) {
                 for (PetscInt cp=0; cp<user->ncp; cp++) {
                     if (!strcmp(propval[propctr], user->componentname[cp])) currentnucleus->activesolutes[cp] = 1;
                 }
             }
         } else {
             currentnucleus->activesolutes = malloc(user->ncp*sizeof(char));
             for (PetscInt cp=0; cp<user->ncp; cp++) currentnucleus->activesolutes[cp] = 1;
         }
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
             /* min incubation time */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "incubation_time", buffer, filesize);
             if (propsize) {assert(propsize == 1); currentcntnuc->incubationtimemin = atof(propval[0]);} 
             else {currentcntnuc->incubationtimemin = 0.0;}
             /* liquidus temperature */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "liquidus_temperature", buffer, filesize);
             if (propsize) {assert(propsize == 1); currentcntnuc->liquidus = atof(propval[0]);} 
             else {currentcntnuc->liquidus = LARGE;}                
         } else if (!strcmp(propval[0], "constant" )) {
             currentnucleus->nuc_model = CONST_NUCLEATION;
             CONST_NUC *currentconstnuc = &currentnucleus->nucleation.constant;
             /* nucleation rate */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "nucleation_rate", buffer, filesize);
             assert(propsize == 1); currentconstnuc->nucleation_rate = atof(propval[0]);
             /* min incubation time */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "incubation_time", buffer, filesize);
             if (propsize) {assert(propsize == 1); currentconstnuc->incubationtimemin = atof(propval[0]);} 
             else {currentconstnuc->incubationtimemin = 0.0;}
         } else if (!strcmp(propval[0], "thermal" )) {
             currentnucleus->nuc_model = THERMAL_NUCLEATION;
             THERMAL_NUC *currentthermalnuc = &currentnucleus->nucleation.thermal;
             /* solvus temperature */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "solvus_temperature_0", buffer, filesize);
             assert(propsize == 1); currentthermalnuc->solvus_temperature_0 = atof(propval[0]);
             ierr = GetProperty(propval, &propsize, nucleusmapping, "solvus_temperature_c", buffer, filesize);
             if (propsize) {
                 assert(propsize == user->ncp); currentthermalnuc->solvus_temperature_c = malloc(user->ncp*sizeof(PetscReal));
                 for (PetscInt propctr = 0; propctr < propsize; propctr++) currentthermalnuc->solvus_temperature_c[propctr] = atof(propval[propctr]);
             } else {
                 currentthermalnuc->solvus_temperature_c = NULL;
             }    
             /* enthalpy of fusion */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "enthalpy_fusion_0", buffer, filesize);
             assert(propsize == 1); currentthermalnuc->enthalpy_fusion_0 = atof(propval[0]);
             ierr = GetProperty(propval, &propsize, nucleusmapping, "enthalpy_fusion_c", buffer, filesize);
             if (propsize) {
                 assert(propsize == user->ncp); currentthermalnuc->enthalpy_fusion_c = malloc(user->ncp*sizeof(PetscReal));
                 for (PetscInt propctr = 0; propctr < propsize; propctr++) currentthermalnuc->enthalpy_fusion_c[propctr] = atof(propval[propctr]);
             } else {
                 currentthermalnuc->enthalpy_fusion_c = NULL;
             }    
             /* surface energy */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "gamma", buffer, filesize);
             assert(propsize == 1); currentthermalnuc->gamma = atof(propval[0]);
             /* heterogeneous nucleation shape factor */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "shape_factor", buffer, filesize);
             if (propsize) {assert(propsize == 1); currentthermalnuc->shapefactor = atof(propval[0]);} 
             else {currentthermalnuc->shapefactor = 1.0;}
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
             /* min incubation time */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "incubation_time", buffer, filesize);
             if (propsize) {assert(propsize == 1); currentthermalnuc->incubationtimemin = atof(propval[0]);} 
             else {currentthermalnuc->incubationtimemin = 0.0;}
             /* liquidus temperature */
             ierr = GetProperty(propval, &propsize, nucleusmapping, "liquidus_temperature", buffer, filesize);
             if (propsize) {assert(propsize == 1); currentthermalnuc->liquidus = atof(propval[0]);} 
             else {currentthermalnuc->liquidus = LARGE;}                 
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
         /* sources */
         {
          char sourcemapping[PETSC_MAX_PATH_LEN];
          sprintf(sourcemapping,"%s/sources",interfacemapping);
          ierr = GetProperty(propval, &propsize, sourcemapping, "activesources", buffer, filesize);
          currentinterface->isources.nsources = propsize;
          currentinterface->isources.sourcename = malloc(currentinterface->isources.nsources*sizeof(char *));
          for (PetscInt s=0; s<currentinterface->isources.nsources; s++) {
              currentinterface->isources.sourcename[s] = malloc(PETSC_MAX_PATH_LEN);
              strcpy(currentinterface->isources.sourcename[s],propval[s]);
          }    
          currentinterface->isources.source = (SOURCE *) malloc(currentinterface->isources.nsources*sizeof(struct SOURCE));
          SOURCE *currentsource = &currentinterface->isources.source[0];
          for (PetscInt s=0; s<currentinterface->isources.nsources; s++,currentsource++) {
              char activesourcemapping[PETSC_MAX_PATH_LEN];
              sprintf(activesourcemapping,"%s/%s",sourcemapping,currentinterface->isources.sourcename[s]);
              ierr = GetProperty(propval, &propsize, activesourcemapping, "source", buffer, filesize); CHKERRQ(ierr);
              assert(propsize == 1);
              /* equilibration source model */
              if        (!strcmp(propval[0], "sink"   )) {
                  currentsource->source_model = SINK_SOURCE;
                  SOURCE_SINK *currentsinksource = &currentsource->source.sink;
                  /* equilibration rate */
                  ierr = GetProperty(propval, &propsize, activesourcemapping, "rate", buffer, filesize);
                  currentsinksource->rate = malloc(user->ncp*sizeof(PetscReal));
                  if (propsize) {
                      assert(propsize == user->ncp);
                      for (PetscInt propctr = 0; propctr < propsize; propctr++) {
                          currentsinksource->rate[propctr] = atof(propval[propctr]);
                      }    
                  } else {
                      memset(currentsinksource->rate,0,user->ncp*sizeof(PetscReal));
                  }
                  /* equilibrium composition */
                  ierr = GetProperty(propval, &propsize, activesourcemapping, "ceq", buffer, filesize);
                  currentsinksource->ceq = malloc(user->ncp*sizeof(PetscReal));
                  if (propsize) {
                      assert(propsize == user->ncp);
                      for (PetscInt propctr = 0; propctr < propsize; propctr++) {
                          currentsinksource->ceq[propctr] = atof(propval[propctr]);
                      }    
                  } else {
                      memset(currentsinksource->ceq,0,user->ncp*sizeof(PetscReal));
                  }
              /* random fluctuation source model */
              } else if (!strcmp(propval[0], "random" )) {
                  currentsource->source_model = RND_SOURCE;
                  SOURCE_RND *currentrndsource = &currentsource->source.rnd;
                  /* fluctuation amplitude */
                  ierr = GetProperty(propval, &propsize, activesourcemapping, "fluctuation_amplitute", buffer, filesize);
                  currentrndsource->fluctuation_amplitute = malloc(user->ncp*sizeof(PetscReal));
                  if (propsize) {
                      for (PetscInt propctr = 0; propctr < propsize; propctr++) {
                          currentrndsource->fluctuation_amplitute[propctr] = atof(propval[propctr]);
                      }    
                  } else {
                      memset(currentrndsource->fluctuation_amplitute,0,user->ncp*sizeof(PetscReal));
                  } 
              }
          }
         } 
         /* interface mobility */
         {
          char mobmapping[PETSC_MAX_PATH_LEN];
          currentinterface->imobility = (IMOBILITY *) malloc(sizeof(IMOBILITY));
          currentinterface->imobility->m = (TACTIVATIONPROP *) malloc(sizeof(TACTIVATIONPROP));
          sprintf(mobmapping, "%s/mobility",interfacemapping);
          ierr = GetProperty(propval, &propsize, mobmapping, "m0", buffer, filesize); CHKERRQ(ierr);
          assert(propsize == 1); currentinterface->imobility->m->m0 = atof(propval[0]);
          currentinterface->imobility->m->unary = (TSeries *) malloc(sizeof(TSeries));
          sprintf(mobmapping, "%s/mobility/activation_energy",interfacemapping);
          ierr = GetProperty(propval, &propsize, mobmapping, "t_coefficient", buffer, filesize);
          currentinterface->imobility->m->unary->nTser = propsize;
          for (PetscInt propctr = 0; propctr < propsize; propctr++) currentinterface->imobility->m->unary->coeff[propctr] = atof(propval[propctr]);
          ierr = GetProperty(propval, &propsize, mobmapping, "t_exponent", buffer, filesize);
          assert(propsize == currentinterface->imobility->m->unary->nTser);
          for (PetscInt propctr = 0; propctr < propsize; propctr++) currentinterface->imobility->m->unary->exp[propctr] = atoi(propval[propctr]);
          currentinterface->imobility->m->unary->logCoeff = 0.0;
          sprintf(mobmapping, "%s/mobility",interfacemapping);
          ierr = GetProperty(propval, &propsize, mobmapping, "anisotropy_values", buffer, filesize);
          if (propsize) {
              assert(propsize > 0); 
              currentinterface->imobility->n = propsize;
              currentinterface->imobility->val = malloc(currentinterface->imobility->n*sizeof(PetscReal));
              for (PetscInt propctr = 0; propctr < propsize; propctr++) currentinterface->imobility->val[propctr] = atof(propval[propctr]);
          } else {
              currentinterface->imobility->n = 0;
          }    
          ierr = GetProperty(propval, &propsize, mobmapping, "anisotropy_directions", buffer, filesize);
          if (propsize) {
              assert(propsize == user->dim*currentinterface->imobility->n); 
              currentinterface->imobility->dir = malloc(user->dim*currentinterface->imobility->n*sizeof(PetscReal));
              for (PetscInt propctr = 0; propctr < propsize; propctr += user->dim) {
                  PetscReal vecnorm = 0.0;
                  for (PetscInt dim = 0; dim < user->dim; dim++) {
                      currentinterface->imobility->dir[propctr+dim] = atof(propval[propctr+dim]);
                      vecnorm += currentinterface->imobility->dir[propctr+dim]*currentinterface->imobility->dir[propctr+dim];
                  }
                  vecnorm = sqrt(vecnorm);
                  for (PetscInt dim = 0; dim < user->dim; dim++) {
                      assert(vecnorm > 0.0);
                      currentinterface->imobility->dir[propctr+dim] /= vecnorm;
                  }
              }    
          }    
         }
         /* interface energy */
         {
          char energymapping[PETSC_MAX_PATH_LEN];
          currentinterface->ienergy = (IENERGY *) malloc(sizeof(IENERGY));
          currentinterface->ienergy->e = (TACTIVATIONPROP *) malloc(sizeof(TACTIVATIONPROP));
          sprintf(energymapping, "%s/energy",interfacemapping);
          ierr = GetProperty(propval, &propsize, energymapping, "e0", buffer, filesize); CHKERRQ(ierr);
          assert(propsize == 1); currentinterface->ienergy->e->m0 = atof(propval[0]);
          currentinterface->ienergy->e->unary = (TSeries *) malloc(sizeof(TSeries));
          sprintf(energymapping, "%s/energy/activation_energy",interfacemapping);
          ierr = GetProperty(propval, &propsize, energymapping, "t_coefficient", buffer, filesize);
          currentinterface->ienergy->e->unary->nTser = propsize;
          for (PetscInt propctr = 0; propctr < propsize; propctr++) currentinterface->ienergy->e->unary->coeff[propctr] = atof(propval[propctr]);
          ierr = GetProperty(propval, &propsize, energymapping, "t_exponent", buffer, filesize);
          assert(propsize == currentinterface->ienergy->e->unary->nTser);
          for (PetscInt propctr = 0; propctr < propsize; propctr++) currentinterface->ienergy->e->unary->exp[propctr] = atoi(propval[propctr]);
          currentinterface->ienergy->e->unary->logCoeff = 0.0;
          sprintf(energymapping, "%s/energy",interfacemapping);
          ierr = GetProperty(propval, &propsize, energymapping, "anisotropy_values", buffer, filesize);
          if (propsize) {
              assert(propsize > 0); 
              currentinterface->ienergy->n = propsize;
              currentinterface->ienergy->val = malloc(currentinterface->ienergy->n*sizeof(PetscReal));
              for (PetscInt propctr = 0; propctr < propsize; propctr++) currentinterface->ienergy->val[propctr] = atof(propval[propctr]);
          } else {
              currentinterface->ienergy->n = 0;
          }    
          ierr = GetProperty(propval, &propsize, energymapping, "anisotropy_directions", buffer, filesize);
          if (propsize) {
              assert(propsize == user->dim*currentinterface->ienergy->n); 
              currentinterface->ienergy->dir = malloc(user->dim*currentinterface->ienergy->n*sizeof(PetscReal));
              for (PetscInt propctr = 0; propctr < propsize; propctr += user->dim) {
                  PetscReal vecnorm = 0.0;
                  for (PetscInt dim = 0; dim < user->dim; dim++) {
                      currentinterface->ienergy->dir[propctr+dim] = atof(propval[propctr+dim]);
                      vecnorm += currentinterface->ienergy->dir[propctr+dim]*currentinterface->ienergy->dir[propctr+dim];
                  }
                  vecnorm = sqrt(vecnorm);
                  for (PetscInt dim = 0; dim < user->dim; dim++) {
                      assert(vecnorm > 0.0);
                      currentinterface->ienergy->dir[propctr+dim] /= vecnorm;
                  }
              }    
          }    
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
         /* interface width */
         ierr = GetProperty(propval, &propsize, interfacemapping, "width", buffer, filesize);
         assert(propsize == 1); 
         currentinterface->width = atof(propval[0]);
     }    
    }

    /* Parsing config file boundary conditions */
    {
     BOUNDARYCONDITIONS *currentboundary = &user->bcs[0];
     for (PetscInt bc=0; bc<user->nbcs; bc++,currentboundary++) {
         char boundarymapping[PETSC_MAX_PATH_LEN];
         sprintf(boundarymapping,"boundary/%s",user->bcname[bc]);
         /* boundary ID */
         ierr = GetProperty(propval, &propsize, boundarymapping, "boundary_id", buffer, filesize);
         assert(propsize == 1); currentboundary->boundaryid = atoi(propval[0]);
         /* chem boundary conditions */
         {
          char chemboundarymapping[PETSC_MAX_PATH_LEN];
          sprintf(chemboundarymapping,"%s/chem",boundarymapping);
          /* boundary type */
          currentboundary->chem_bctype = NONE_BC;
          ierr = GetProperty(propval, &propsize, chemboundarymapping, "type", buffer, filesize);
          if (propsize) {
              assert(propsize == 1); 
              if      (!strcmp(propval[0], "neumann"  )) {currentboundary->chem_bctype = NEUMANN_BC;  }
              else if (!strcmp(propval[0], "dirichlet")) {currentboundary->chem_bctype = DIRICHLET_BC;}
              else                                       {currentboundary->chem_bctype = NONE_BC;     }
          }
          /* boundary val */
          ierr = GetProperty(propval, &propsize, chemboundarymapping, "value", buffer, filesize);
          if (!(currentboundary->chem_bctype == NONE_BC)) {
              assert(propsize == user->ncp); 
              currentboundary->chem_bcval = malloc(user->ncp*sizeof(PetscReal));
              currentboundary->chem_bcbool = malloc(user->ncp*sizeof(char));
              for (PetscInt propctr = 0; propctr < propsize; propctr++) {
                  if (!strcmp(propval[propctr], "*" )) {
                      currentboundary->chem_bcbool[propctr] = 0;
                      currentboundary->chem_bcval[propctr] = 0.0;
                  } else {
                      currentboundary->chem_bcbool[propctr] = 1;
                      currentboundary->chem_bcval[propctr] = atof(propval[propctr]);
                  }
              }    
          }    
         }
         /* thermal boundary conditions */
         {
          char thermalboundarymapping[PETSC_MAX_PATH_LEN];
          sprintf(thermalboundarymapping,"%s/thermal",boundarymapping);
          /* boundary type */
          currentboundary->thermal_bctype = NONE_BC;
          ierr = GetProperty(propval, &propsize, thermalboundarymapping, "type", buffer, filesize);
          if (propsize) {
              assert(propsize == 1); 
              if      (!strcmp(propval[0], "neumann"  )) {currentboundary->thermal_bctype = NEUMANN_BC;  }
              else if (!strcmp(propval[0], "dirichlet")) {currentboundary->thermal_bctype = DIRICHLET_BC;}
              else                                       {currentboundary->thermal_bctype = NONE_BC;     }
          }
          /* boundary val */
          ierr = GetProperty(propval, &propsize, thermalboundarymapping, "value", buffer, filesize);
          if (!(currentboundary->thermal_bctype == NONE_BC)) {
              assert(propsize == 1); currentboundary->thermal_bcval = atof(propval[0]);
          }    
         } 
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
     assert(propsize > 0); user->solparams.nloadcases = propsize; 
     user->solparams.time = malloc(user->solparams.nloadcases*sizeof(PetscReal));
     for (PetscInt propctr = 0; propctr < propsize; propctr++) user->solparams.time[propctr] = atof(propval[propctr]);
     ierr = GetProperty(propval, &propsize, "solution_parameters", "temperature_rate", buffer, filesize);
     if (propsize) {
         assert(propsize == user->solparams.nloadcases);
         user->solparams.temperature_rate = malloc(user->solparams.nloadcases*sizeof(PetscReal));
         for (PetscInt propctr = 0; propctr < propsize; propctr++) user->solparams.temperature_rate[propctr] = atof(propval[propctr]);
     } else {
         user->solparams.temperature_rate = NULL;
     }
     ierr = GetProperty(propval, &propsize, "solution_parameters", "timestep0", buffer, filesize); 
     if (propsize) {assert(propsize == 1); user->solparams.timestep = atof(propval[0]);} 
     else {user->solparams.timestep = 1.0e-18;}
     ierr = GetProperty(propval, &propsize, "solution_parameters", "timestepmin", buffer, filesize);
     if (propsize) {assert(propsize == 1); user->solparams.mintimestep = atof(propval[0]);} 
     else {user->solparams.mintimestep = 1.0e-18;}
     ierr = GetProperty(propval, &propsize, "solution_parameters", "timestepmax", buffer, filesize);
     if (propsize) {assert(propsize == 1); user->solparams.maxtimestep = atof(propval[0]);} 
     else {user->solparams.maxtimestep = user->solparams.time[0];}
     ierr = GetProperty(propval, &propsize, "solution_parameters", "junctionpenalty", buffer, filesize);
     if (propsize) {assert(propsize == 1); user->solparams.junctionpenalty = atof(propval[0]);} 
     else {user->solparams.junctionpenalty = 1.0;}
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
     ierr = GetProperty(propval, &propsize, "solution_parameters", "gradient_calculation", buffer, filesize);
     if (propsize) {assert(propsize == 1 && atoi(propval[0]) >= 0 && atoi(propval[0]) <= 1); user->gradient_calculation = atoi(propval[0]);}
     else {user->gradient_calculation = 0;}
     ierr = GetProperty(propval, &propsize, "solution_parameters", "max_phase_neighbours", buffer, filesize);
     if (propsize) {assert(propsize == 1); max_phase_neighbours = atoi(propval[0]);}
    }

    /* field offsets */
    SF_SIZE = MAXSITES*user->ncp;
    SP_SIZE = MAXSITES*user->ndp;

    AS_OFFSET = 0;
    AS_SIZE   = (max_phase_neighbours < user->npf ? max_phase_neighbours : user->npf) + 1;
    PF_OFFSET = AS_OFFSET+AS_SIZE;
    PF_SIZE   = (max_phase_neighbours < user->npf ? max_phase_neighbours : user->npf);
    DP_OFFSET = PF_OFFSET+PF_SIZE;
    DP_SIZE   = user->ndp;
    EX_OFFSET = DP_OFFSET+DP_SIZE;
    EX_SIZE   = PF_SIZE*SP_SIZE;
    TM_OFFSET = EX_OFFSET+EX_SIZE;
    TM_SIZE   = 1;
    
    /* Parsing config file mappings */
    {
     if (user->worldrank == 0) {
        printf("Parsing config file mappings...\n");
     }     
     char *tok, *savetok;
     PetscInt ctrm = 0;
     if (user->worldrank == 0) {             
        ierr = GetProperty(propval, &propsize, "mappings", "phase_material_mapping", buffer, filesize); CHKERRQ(ierr);
        assert(propsize);
        tok = strtok_r(propval[0], "\n", &savetok);
        if (strstr(tok, ".txt") != NULL) {
            printf("Reading in file '%s' to populate the phase-material mapping.\n", tok);
            FILE *mappingfile=NULL;
            mappingfile = fopen (propval[0], "r");
            if (mappingfile==NULL) {
                printf("Error: Voxel phase mapping file can't be opened. \n");
                printf("Error: Please make sure the file %s is in the CWD.\n",propval[0]);
                return 1;
            }
            long int mappingfilesize;
            if (mappingfile) {
                fseek(mappingfile, 0, SEEK_END);
                mappingfilesize = ftell(infile);
                fseek(mappingfile, 0, SEEK_SET);
            }
            mappingbuffer = malloc(mappingfilesize);
            if (mappingfile) {
                fread(mappingbuffer, 1, mappingfilesize, mappingfile);
                fclose(mappingfile);
            }
            tok = strtok_r(mappingbuffer, "\n", &savetok);
        }
        printf("Processing the phase-material mapping.\n");
        ctrm=0;
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
        if (mappingbuffer) {free(mappingbuffer);mappingbuffer=0;}
        printf("Finished processing the phase-material mapping.\n");
     }
     ierr = MPI_Bcast(user->phasematerialmapping, user->npf, MPIU_INT, 0, PETSC_COMM_WORLD);

     if (user->worldrank == 0) {        
        ierr = GetProperty(propval, &propsize, "mappings", "voxel_phase_mapping", buffer, filesize); CHKERRQ(ierr);
        assert(propsize);
        tok = strtok_r(propval[0], "\n", &savetok);
        if (strstr(tok, ".txt") != NULL) {
            printf("Reading in file '%s' to populate the voxel-phase mapping.\n", tok);
            FILE *mappingfile=NULL;
            mappingfile = fopen (propval[0], "r");
            if (mappingfile==NULL) {
                printf("Error: Voxel phase mapping file can't be opened. \n");
                printf("Error: Please make sure the file %s is in the CWD.\n",propval[0]);
                return 1;
            }
            long int mappingfilesize;
            if (mappingfile) {
                fseek(mappingfile, 0, SEEK_END);
                mappingfilesize = ftell(infile);
                fseek(mappingfile, 0, SEEK_SET);
            }
            mappingbuffer = malloc(mappingfilesize);
            if (mappingfile) {
                fread(mappingbuffer, 1, mappingfilesize, mappingfile);
                fclose(mappingfile);
            }
            tok = strtok_r(mappingbuffer, "\n", &savetok);
        }
        printf("Processing the voxel-phase mapping.\n");
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
        if (mappingbuffer) {free(mappingbuffer);mappingbuffer=0;}
        printf("Finished processing voxel-phase mapping.\n");
     }
     ierr = MPI_Bcast(user->voxelphasemapping, total_cells, MPIU_INT, 0, PETSC_COMM_WORLD);

     if (user->worldrank == 0) {
        ierr = GetProperty(propval, &propsize, "mappings", "interface_mapping", buffer, filesize); CHKERRQ(ierr);
        assert(propsize);
        tok = strtok_r(propval[0], "\n", &savetok);
        if (strstr(tok, ".txt") != NULL) {
            printf("Reading in file '%s' to populate the interface mapping.\n", tok);
            FILE *mappingfile=NULL;
            mappingfile = fopen (propval[0], "r");
            if (mappingfile==NULL) {
                printf("Error: Interface mapping file can't be opened. \n");
                printf("Error: Please make sure the file %s is in the CWD.\n",propval[0]);
                return 1;
            }
            long int mappingfilesize;
            if (mappingfile) {
                fseek(mappingfile, 0, SEEK_END);
                mappingfilesize = ftell(infile);
                fseek(mappingfile, 0, SEEK_SET);
            }
            mappingbuffer = malloc(mappingfilesize);
            if (mappingfile) {
                fread(mappingbuffer, 1, mappingfilesize, mappingfile);
                fclose(mappingfile);
            }
            tok = strtok_r(mappingbuffer, "\n", &savetok);
        }
        printf("Processing the interface mapping.\n");
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
        if (mappingbuffer) {free(mappingbuffer);mappingbuffer=0;}
        printf("Finished processing the interface mapping.\n");
     }
     ierr = MPI_Bcast(user->interfacelist, user->npf * user->npf, MPIU_INT, 0, PETSC_COMM_WORLD);

     if (user->nsites) {
         ierr = GetProperty(propval, &propsize, "mappings", "site_nucleus_mapping", buffer, filesize);
         assert(propsize);
         tok = strtok_r(propval[0], "\n", &savetok);
		 if (strstr(tok, ".txt") != NULL) {
			 FILE *mappingfile=NULL;
			 mappingfile = fopen (propval[0], "r");
			 if (mappingfile==NULL) {
				 printf("Error: Interface mapping file can't be opened. \n");
				 printf("Error: Please make sure the file %s is in the CWD.\n",propval[0]);
				 return 1;
			 }
			 long int mappingfilesize;
			 if (mappingfile) {
				 fseek(mappingfile, 0, SEEK_END);
				 mappingfilesize = ftell(infile);
				 fseek(mappingfile, 0, SEEK_SET);
			 }
			 mappingbuffer = malloc(mappingfilesize);
			 if (mappingfile) {
				 fread(mappingbuffer, 1, mappingfilesize, mappingfile);
				 fclose(mappingfile);
			 }
			 tok = strtok_r(mappingbuffer, "\n", &savetok);
		 }
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
	     if (mappingbuffer) {free(mappingbuffer);mappingbuffer=0;}

         ierr = GetProperty(propval, &propsize, "mappings", "site_phase_mapping", buffer, filesize);
         assert(propsize);
         tok = strtok_r(propval[0], "\n", &savetok);
		 if (strstr(tok, ".txt") != NULL) {
			 FILE *mappingfile=NULL;
			 mappingfile = fopen (propval[0], "r");
			 if (mappingfile==NULL) {
				 printf("Error: Interface mapping file can't be opened. \n");
				 printf("Error: Please make sure the file %s is in the CWD.\n",propval[0]);
				 return 1;
			 }
			 long int mappingfilesize;
			 if (mappingfile) {
				 fseek(mappingfile, 0, SEEK_END);
				 mappingfilesize = ftell(infile);
				 fseek(mappingfile, 0, SEEK_SET);
			 }
			 mappingbuffer = malloc(mappingfilesize);
			 if (mappingfile) {
				 fread(mappingbuffer, 1, mappingfilesize, mappingfile);
				 fclose(mappingfile);
			 }
			 tok = strtok_r(mappingbuffer, "\n", &savetok);
		 }
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
	     if (mappingbuffer) {free(mappingbuffer);mappingbuffer=0;}

         ierr = GetProperty(propval, &propsize, "mappings", "voxel_site_mapping", buffer, filesize);
         assert(propsize);
         tok = strtok_r(propval[0], "\n", &savetok);
		 if (strstr(tok, ".txt") != NULL) {
			 FILE *mappingfile=NULL;
			 mappingfile = fopen (propval[0], "r");
			 if (mappingfile==NULL) {
				 printf("Error: Interface mapping file can't be opened. \n");
				 printf("Error: Please make sure the file %s is in the CWD.\n",propval[0]);
				 return 1;
			 }
			 long int mappingfilesize;
			 if (mappingfile) {
				 fseek(mappingfile, 0, SEEK_END);
				 mappingfilesize = ftell(infile);
				 fseek(mappingfile, 0, SEEK_SET);
			 }
			 mappingbuffer = malloc(mappingfilesize);
			 if (mappingfile) {
				 fread(mappingbuffer, 1, mappingfilesize, mappingfile);
				 fclose(mappingfile);
			 }
			 tok = strtok_r(mappingbuffer, "\n", &savetok);
		 }
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
	     if (mappingbuffer) {free(mappingbuffer);mappingbuffer=0;}
     }
    }     
    free(buffer); 
    if (user->worldrank == 0) {
        printf("Finished parsing config file.\n");
    }    
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
    PetscInt          localcell, cell, face, phase, g, c, s;
    PetscInt          conesize, supp, nsupp;
    const PetscInt    *cone, *scells;
    DMLabel           plabel = NULL;
    uint16_t          gslist[AS_SIZE], nalist[AS_SIZE];
    PetscScalar       *pcell, *dcell, *ccell, *tcell;
    uint16_t          setunion[AS_SIZE], injectionL[AS_SIZE], injectionR[AS_SIZE];
    MATERIAL          *currentmaterial;
    PetscReal         sitepot[SP_SIZE];

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
        ccell  = &offset[EX_OFFSET];
        tcell  = &offset[TM_OFFSET];
        PetscInt phase;
        ierr = DMLabelGetValue(plabel, cell, &phase);
    
        /* set initial conditions */
        for (g=0; g<gslist[0]; g++) {
            if (gslist[g+1] == phase) {
                currentmaterial = &user->material[user->phasematerialmapping[gslist[g+1]]];
                *tcell = currentmaterial->thermal.temperature0;
            }
        }    
        for (g=0; g<gslist[0]; g++) {
            currentmaterial = &user->material[user->phasematerialmapping[gslist[g+1]]];
            SitepotentialExplicit(&ccell[g*SP_SIZE],currentmaterial->c0,(*tcell),gslist[g+1],user);
            if (gslist[g+1] == phase) {
                pcell[g] = 1.0;
                Sitepotential(sitepot,currentmaterial->c0,(*tcell),gslist[g+1],user);
                memset(dcell,0,user->ndp*sizeof(PetscReal));
                for (s=0; s<currentmaterial->nsites; s++) {
                    for (c=0; c<user->ndp; c++) {
                        dcell[c] += sitepot[s*user->ndp+c];
                    }
                }
            } else {
                pcell[g] = 0.0;
            }
        }
    }        
    ierr = VecRestoreArray(solution, &fdof);
    PetscFunctionReturn(0);
}

