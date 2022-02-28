#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "edefinition.h"
#include "eutility.h"




int ReadParFile(int ARGC, char *ARGV[], simulationparameters * simpa_p)
{
  /* Declaration */
  FILE * parfile_p;
  int count = 0;
  int i;
  struct tm *simtime;
  char dump[101];
   /* Code */
  printf("\n");
  printf("***************************\n");
  DBG(CALL OF ReadParFile);
  printf("***************************\n");
  if(ARGC>1){
    strcpy(simpa_p->ParFileName,ARGV[1]);
  }
  else {
    strcpy(simpa_p->ParFileName,"parfile2.txt");
  }
  
  if((parfile_p = fopen(simpa_p->ParFileName,"r"))==NULL){
    printf("%s:%d Error opening parameter file %s.\n",__FILE__,__LINE__,simpa_p->ParFileName);
    return -1;
  }
  fgets(dump,100,parfile_p);
  count = fscanf(parfile_p,"Simulation Tag: %s\n",simpa_p->tag);
  count += fscanf(parfile_p,"Number of Particles: %ld\n",&simpa_p->npart);
  count += fscanf(parfile_p,"Number of Intervals in probability table: %lf\n", &simpa_p->q);
  count += fscanf(parfile_p,"Number of Intervals in (0 PI/2) for reflection/conversion tables: %d\n",&simpa_p->Nangle);
  count += fscanf(parfile_p,"Grid size (x, y, z, theta, t): %d %d %d %d %d\n",&simpa_p->GridSize[0],&simpa_p->GridSize[1],&simpa_p->GridSize[2],&simpa_p->GridSize[3],&simpa_p->GridSize[4]);
  count += fscanf(parfile_p,"Grid spacing (x, y, z, theta, t) in m, degree or ms: %lf %lf %lf %lf %lf\n",&simpa_p->GridSpacing[0], &simpa_p->GridSpacing[1], &simpa_p->GridSpacing[2], &simpa_p->GridSpacing[3], &simpa_p->GridSpacing[4]);
  count += fscanf(parfile_p,"Simulation time interval in ms: %lf\n",&simpa_p->TimeStep);
  count += fscanf(parfile_p,"Simulation length in ms: %f\n",&simpa_p->SimLen);
  count += fscanf(parfile_p,"Reference 1/Q (internal amplification): %f\n",&simpa_p->Qref);
  count += fscanf(parfile_p,"Output type (0 for 4D-data, 1 for time-distance section): %d\n",&simpa_p->outputtype);

  fgets(dump,100,parfile_p);
  count += fscanf(parfile_p,"Source position (x, y, z) in m from grid origin: %lf %lf %lf\n",&simpa_p->SourcePos[0], &simpa_p->SourcePos[1], &simpa_p->SourcePos[2]);
  count += fscanf(parfile_p,"Source type (0 for P, 1 for S-Source): %d\n",&simpa_p->SourceType);
  count += fscanf(parfile_p,"Source orientation (strike, dip, rake) in degrees: %lf %lf %lf\n",&simpa_p->SourceOri[0], &simpa_p->SourceOri[1], &simpa_p->SourceOri[2]); 
  count += fscanf(parfile_p,"Angular frequency *1000: %lf\n", &simpa_p->frc);
  count += fscanf(parfile_p,"Source body: %d\n", &simpa_p->sourcebody);

  fgets(dump,100,parfile_p);
  count += fscanf(parfile_p,"File for Fractional fluctuation e: %s\n", simpa_p->Espsilon_file);
  count += fscanf(parfile_p,"File for Intrinsic Attenuation 1/Q_p: %s\n", simpa_p->Qp_file);
  count += fscanf(parfile_p,"File for Intrinsic Attenuation 1/Q_s: %s\n", simpa_p->Qs_file);
  count += fscanf(parfile_p,"Number of bodies: %d\n",&simpa_p->Nbody);
  simpa_p->body = (bodyparm *) calloc((size_t) simpa_p->Nbody, sizeof(bodyparm));
  for(i=0;i<simpa_p->Nbody;i++){
    fgets(dump,100,parfile_p);
    count += fscanf(parfile_p,"X coordinate: %lf %lf\n",&simpa_p->body[i].x[0], &simpa_p->body[i].x[1]);
    count += fscanf(parfile_p,"Y coordinate: %lf %lf\n",&simpa_p->body[i].y[0], &simpa_p->body[i].y[1]);
    count += fscanf(parfile_p,"Z coordinate: %lf %lf\n",&simpa_p->body[i].z[0], &simpa_p->body[i].z[1]);
    count += fscanf(parfile_p,"S Velocity and gamma: %lf %lf\n",&simpa_p->body[i].v[1], &simpa_p->body[i].gam);
    count += fscanf(parfile_p,"Density: %lf\n",&simpa_p->body[i].rho);
    count += fscanf(parfile_p,"Type of ACF (gauss, expo, karman): %s\n",simpa_p->body[i].type);
    count += fscanf(parfile_p,"Kappa (only evaluated for karnam ACF): %lf\n",&simpa_p->body[i].kap);
    count += fscanf(parfile_p,"Correlation distance a: %lf\n",&simpa_p->body[i].a);
    count += fscanf(parfile_p,"Fractional fluctuation e: %lf\n",&simpa_p->body[i].e);
    count += fscanf(parfile_p,"Birch ny: %lf\n",&simpa_p->body[i].ny);
    count += fscanf(parfile_p,"Intrinsic Attenuation (1/Q_p 1/Q_s): %lf %lf\n",&simpa_p->body[i].Qi[0], &simpa_p->body[i].Qi[1]);
    simpa_p->body[i].v[0] = simpa_p->body[i].v[1] * simpa_p->body[i].gam;
    simpa_p->body[i].l = simpa_p->frc / simpa_p->body[i].v[1];
    simpa_p->body[i].q = simpa_p->q;
    simpa_p->body[i].path[0] = simpa_p->body[i].v[0] * simpa_p->TimeStep;
    simpa_p->body[i].path[1] = simpa_p->body[i].v[1] * simpa_p->TimeStep;
    if(!strcmp(simpa_p->body[i].type,"karman")){
      sprintf(simpa_p->body[i].ProbTabFileName,"tables/%s_%1.2f_%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_%1.2e_scat",simpa_p->body[i].type,simpa_p->body[i].kap,simpa_p->body[i].gam,simpa_p->body[i].l, simpa_p->body[i].a, 1., simpa_p->body[i].ny, simpa_p->q);
    }
    else{
      sprintf(simpa_p->body[i].ProbTabFileName,"tables/%s_%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_%1.2e_scat",simpa_p->body[i].type,simpa_p->body[i].gam,simpa_p->body[i].l, simpa_p->body[i].a, 1., simpa_p->body[i].ny, simpa_p->q);
    }
  }
    if(simpa_p->Nbody==1){
     sprintf(simpa_p->InterfaceName,"tables/interf1_%1.4f_%1.4f_0_0_%1.4f_0_%1.2e.int",simpa_p->body[0].v[0],simpa_p->body[0].v[1],simpa_p->body[0].rho,(float) simpa_p->Nangle);}
  else{
  sprintf(simpa_p->InterfaceName,"tables/interf1_%1.4f_%1.4f_%1.4f_%1.4f_%1.4f_%1.4f_%1.2e.int",simpa_p->body[0].v[0],simpa_p->body[0].v[1],simpa_p->body[1].v[0],simpa_p->body[1].v[1],simpa_p->body[0].rho,simpa_p->body[1].rho,(float) simpa_p->Nangle);}
  sprintf(simpa_p->SurfaceName,"tables/surf1_%1.4f_%1.4f_%1.4f_%1.2e.int",simpa_p->body[0].v[0],simpa_p->body[0].v[1],simpa_p->body[0].rho,(float) simpa_p->Nangle);

printf("\n"); 
printf("***************************\n");
DBG(MARK2);
printf("***************************\n");
printf("\n"); 


  simpa_p->SimulationTime = time(0);
  simtime = localtime(&simpa_p->SimulationTime);
  sprintf(simpa_p->dstring,"%s_%d_%d_%d-%d_%d_%d",simpa_p->tag,simtime->tm_year+1900,simtime->tm_mon+1,simtime->tm_mday,simtime->tm_hour,simtime->tm_min,simtime->tm_sec);
  return 0;
}

int ReadModelFiles(simulationparameters * simpa)
{
    int i,j;
    printf("\n");
    printf("***************************\n");
    DBG(CALL OF ReadModelFile);
    printf("***************************\n");
    /* Allocate space and read the model files for the epsilon*/
    /* x - direction */
    if((simpa->body[0].ane = (double ***) malloc((size_t) simpa->GridSize[0]*sizeof(double **)))==NULL){
        fprintf(stderr,"%s:%d Error allocating memory\n",__FILE__,__LINE__);
        return 0;
    }
    /* y - direction */
    for(i=0;i<simpa->GridSize[0];i++){
        if((simpa->body[0].ane[i] = (double **) malloc((size_t) simpa->GridSize[1]*sizeof(double *)))==NULL){
            fprintf(stderr,"%s:%d Error allocating memory\n",__FILE__,__LINE__);
            return 0;
        }
        /* z - direction */
        for(j=0;j<simpa->GridSize[1];j++){
            if((simpa->body[0].ane[i][j] = (double*) malloc((size_t) simpa->GridSize[2]*sizeof(double)))==NULL){
                fprintf(stderr,"%s:%d Error allocating memory\n",__FILE__,__LINE__);
                return 0;

            }
        }
    }
    
    /* Allocate space and read the model files for the Q*/
    /* x - direction */
    if((simpa->body[0].anQp = (double ***) malloc((size_t) simpa->GridSize[0]*sizeof(double **)))==NULL & (simpa->body[0].anQs = (double ***) malloc((size_t) simpa->GridSize[0]*sizeof(double **)))==NULL){
        fprintf(stderr,"%s:%d Error allocating memory\n",__FILE__,__LINE__);
        return 0;
    }
    /* y - direction */
    for(i=0;i<simpa->GridSize[0];i++){
        if((simpa->body[0].anQp[i] = (double **) malloc((size_t) simpa->GridSize[1]*sizeof(double *)))==NULL & (simpa->body[0].anQs[i] = (double **) malloc((size_t) simpa->GridSize[1]*sizeof(double *))) ==NULL){
            fprintf(stderr,"%s:%d Error allocating memory\n",__FILE__,__LINE__);
            return 0;
        }
        /* z - direction */
        for(j=0;j<simpa->GridSize[1];j++){
            if((simpa->body[0].anQp[i][j] = (double*) malloc((size_t) simpa->GridSize[2]*sizeof(double)))==NULL & (simpa->body[0].anQs[i][j] = (double*) malloc((size_t) simpa->GridSize[2]*sizeof(double)))==NULL){
                fprintf(stderr,"%s:%d Error allocating memory\n",__FILE__,__LINE__);
                return 0;
                
            }
        }
    }
    
    
    /* Read the model file and fill the ane array */
    
    FILE * fp;
    fp = fopen(simpa->Espsilon_file,"r");
    for(i=0;i<simpa->GridSize[0];i++){
        for(j=0;j<simpa->GridSize[1];j++){
            fread(simpa->body[0].ane[i][j],sizeof(double),simpa->GridSize[2],fp);
        }
    }
    
    /* Read the model file and fill the anQp and anQs array */
    fp = fopen(simpa->Qp_file,"r");
    for(i=0;i<simpa->GridSize[0];i++){
        for(j=0;j<simpa->GridSize[1];j++){
            fread(simpa->body[0].anQp[i][j],sizeof(double),simpa->GridSize[2],fp);
        }
    }
    fclose(fp);
    
    fp = fopen(simpa->Qs_file,"r");
    for(i=0;i<simpa->GridSize[0];i++){
        for(j=0;j<simpa->GridSize[1];j++){
            fread(simpa->body[0].anQs[i][j],sizeof(double),simpa->GridSize[2],fp);
        }
    }
    fclose(fp);
    
    printf("\n");
    printf("***************************\n");
    DBG(MARK3);
    printf("***************************\n");
    printf("\n");
    
    /*for(i=0;i<simpa->GridSize[0];i++){
        for(j=0;j<simpa->GridSize[1];j++){
            printf("%lf\t", simpa->body[0].ane[i][j][3]);
        }
        printf("\n");
    }*/
    return 0;
}

gridtype ***** GridAlloc(simulationparameters simpa)
{
  /* Declaration */
  gridtype ***** grid;
  int i,j,k,l;

  /* Code */
  DBG(CALL OF GridAlloc);
  /* x - direction */
  if((grid = (gridtype *****) malloc((size_t) simpa.GridSize[0]*sizeof(gridtype ***)))==NULL){
    fprintf(stderr,"%s:%d Error allocating memory\n",__FILE__,__LINE__);
    return NULL;
  }
  /* y - direction */
  for(i=0;i<simpa.GridSize[0];i++){
    if((grid[i] = (gridtype ****) malloc((size_t) simpa.GridSize[1]*sizeof(gridtype **)))==NULL){
      fprintf(stderr,"%s:%d Error allocating memory\n",__FILE__,__LINE__);
      return NULL;
    }
    /* z - direction */
    for(j=0;j<simpa.GridSize[1];j++){
      if((grid[i][j] = (gridtype ***) malloc((size_t) simpa.GridSize[2]*sizeof(gridtype *)))==NULL){
        fprintf(stderr,"%s:%d Error allocating memory\n",__FILE__,__LINE__);
        return NULL;
      }
      /* theta - angle */
      for(k=0;k<simpa.GridSize[2];k++){
        if((grid[i][j][k] = (gridtype **) malloc((size_t) simpa.GridSize[3]*sizeof(gridtype *)))==NULL){
          fprintf(stderr,"%s:%d Error allocating memory\n",__FILE__,__LINE__);
          return NULL;
        }
        /* time axis */
        for(l=0;l<simpa.GridSize[3];l++){
          if((grid[i][j][k][l] = (gridtype *) calloc((size_t) simpa.GridSize[4], sizeof(gridtype)))==NULL){
            fprintf(stderr,"%s:%d Error allocating memory\n",__FILE__,__LINE__);
            return NULL;
          }
        }
      }
    }
  }
  return grid;
}


int GridFree(simulationparameters simpa, gridtype ***** grid)
{
   /* Declaration */
   int i,j,k,l;

   /* Code */
   DBG(CALL OF GridFree);
   for(i=0;i<simpa.GridSize[0];i++){
      for(j=0;j<simpa.GridSize[1];j++){
         for(k=0;k<simpa.GridSize[2];k++){
             for(l=0;l<simpa.GridSize[3];l++){
                 free(grid[i][j][k][l]);
             }
             free(grid[i][j][k]);
         }
         free(grid[i][j]);
      }
      free(grid[i]);
   }
   free(grid);
   return 0;
}



int WriteOutput(simulationparameters simpa, gridtype ***** P, gridtype ***** S, bodytable * scattab)
    {
  /* Declaration */
  char call[200], call2[200];

  char call_p[200];
  char call_s[200];

  char call_psec[200];
  char call_ssec[200];

  int i,j,k,l,m;
  FILE * fp1;
  FILE * fp3;

  FILE * f_p;
  FILE * f_s;

  FILE * f_psec;
  FILE * f_ssec;

  double * refDamp;
  double QrefExp;
  float tmp;



  double **** Psec;
  double **** Ssec;

  int * stknum;
  double dist,dmin;
  double *x;
  double *y;
  double az;
  float azrange[2] = {-180. , 180.}; /*{-0.5, 0.5};*/
  double GridSpacing, GridSize;
  int gind;

  /* Code */
  printf("***************************\n");
  DBG(CALL OF WriteOutput);
  printf("***************************\n");
  printf("\n");
  refDamp = (double*)calloc((size_t)simpa.GridSize[4],sizeof(double));
  QrefExp = exp(simpa.frc * simpa.Qref * simpa.TimeStep);  /* artificial amplification of particle energy per timestep (introduced with intrinsic attenuation) */
  for(l=0;l<simpa.GridSize[4];l++){
    refDamp[l] = pow(QrefExp,l)*simpa.emodulo*simpa.npart;
  }

/* ---------------------------------------------------------------------------------------------------------------------------------------------------------------- */
/* OUTPUT OF COMPELTE 5D DATA CUBE */
if(simpa.outputtype == 0){    
    sprintf(call_p, "./sims/%s_P.sim",simpa.dstring);
    sprintf(call_s, "./sims/%s_S.sim",simpa.dstring);
    if((f_p      = fopen(call_p, "w")) ==NULL)    printf("Error opening output file %s\n",call_p);
    if((f_s      = fopen(call_s, "w")) ==NULL)    printf("Error opening output file %s\n",call_s);
    /* OUTPUT OF P,S,TL,TQ,TT,RL,RQ,RT*/ 
   for(m=0;m<simpa.GridSize[4];m++){
    for(l=0;l<simpa.GridSize[3];l++){
      for(k=0;k<simpa.GridSize[2];k++){
        for(j=0;j<simpa.GridSize[1];j++){
          for(i=0;i<simpa.GridSize[0];i++){
            tmp  = (float)(((double)P[i][j][k][l][m]));
            fwrite(&tmp,(size_t) sizeof(float), (size_t) 1, f_p);
            tmp  = (float)(((double)S[i][j][k][l][m]));
            fwrite(&tmp,(size_t) sizeof(float), (size_t) 1, f_s);    
           }
         }
       }
     }
    }
    fclose(f_p);                                                                                            
    fclose(f_s); 
}


/* OUTPUT OF SECTION STORING DEPTH, RADIUS AND TIME */
if(simpa.outputtype == 1){ 
    GridSpacing = (double)(simpa.GridSpacing[0] + simpa.GridSpacing[1])/2;
    dmin = MAX(MAX(simpa.GridSize[0]*simpa.GridSpacing[0]-simpa.SourcePos[0],simpa.SourcePos[0]),MAX(simpa.GridSize[1]*simpa.GridSpacing[1]-simpa.SourcePos[1],simpa.SourcePos[1]));
    GridSize = dmin/GridSpacing;
    x = (double *)calloc((size_t)simpa.GridSize[0],(size_t) sizeof(double));
    y = (double *)calloc((size_t)simpa.GridSize[1],(size_t) sizeof(double));
    for(i=0;i<simpa.GridSize[0];i++)x[i]=(0.5+i)*simpa.GridSpacing[0]-simpa.SourcePos[0];
    for(i=0;i<simpa.GridSize[1];i++)y[i]=(0.5+i)*simpa.GridSpacing[1]-simpa.SourcePos[1];
    stknum = (int *)calloc((size_t)GridSize,(size_t)sizeof(int));
    Psec  = (double ****)malloc((size_t)GridSize*sizeof(double***));
    Ssec  = (double ****)malloc((size_t)GridSize*sizeof(double***));
    for(i=0;i<GridSize;i++){
      Psec[i] = (double ***)calloc((size_t)simpa.GridSize[4],(size_t)sizeof(double**));
      Ssec[i] = (double ***)calloc((size_t)simpa.GridSize[4],(size_t)sizeof(double**));
      for(j=0;j<simpa.GridSize[2];j++){
          Psec[i][j]  = (double **)calloc((size_t)simpa.GridSize[4],(size_t)sizeof(double*));
          Ssec[i][j]  = (double **)calloc((size_t)simpa.GridSize[4],(size_t)sizeof(double*));
          for(k=0;k<simpa.GridSize[3];k++){
              Psec[i][j][k]  = (double *)calloc((size_t)simpa.GridSize[4],(size_t)sizeof(double));
              Ssec[i][j][k]  = (double *)calloc((size_t)simpa.GridSize[4],(size_t)sizeof(double));
          }
      }
    }
    for(i=0;i<simpa.GridSize[0];i++){
      for(j=0;j<simpa.GridSize[1];j++){
        dist = sqrt(x[i]*x[i] + y[j]*y[j]);
        gind = (int)round(dist/GridSpacing);
	az = atan2(x[i],y[j]);  /* stack only in defined azimuth range */
        if((gind<GridSize) & (azrange[0]<az) & (azrange[1]>az)){
        if((gind<GridSize)){
          stknum[gind] ++;
	  for(l=0;l<simpa.GridSize[2];l++){
	    for(k=0;k<simpa.GridSize[3];k++){
            for(m=0;m<simpa.GridSize[4];m++){
	      Psec[gind][l][k][m]  += P[i][j][l][k][m];
	      Ssec[gind][l][k][m]  += S[i][j][l][k][m];
            }
           }
          }
        }
      }
    }
	}
    for(i=0;i<GridSize;i++){
      for(m=0;m<simpa.GridSize[4];m++){
        if(stknum[i]>0){
    for(l=0;l<simpa.GridSize[3];l++){
	  for(k=0;k<simpa.GridSize[2];k++){
	    Psec[i][k][l][m]    /= stknum[i];
	    Ssec[i][k][l][m]    /= stknum[i];
	  }
        }
      }
    }
    sprintf(call_psec, "./sims/%s_Psec.sec",simpa.dstring);
    sprintf(call_ssec, "./sims/%s_Ssec.sec",simpa.dstring);

    if((f_psec  = fopen(call_psec, "w"))==NULL){printf("Error opening output file %s\n",call_psec);  return -1;}
    if((f_ssec  = fopen(call_ssec, "w"))==NULL){printf("Error opening output file %s\n",call_ssec);  return -1;}
for(m=0;m<simpa.GridSize[3];m++){
    for(l=0;l<simpa.GridSize[2];l++){
      for(i=0;i<GridSize;i++){
    for(k=0;k<simpa.GridSize[4];k++){
      tmp  = (float)(Psec[i][l][m][k]/refDamp[k]);
      fwrite(&tmp,(size_t) sizeof(float), (size_t) 1, f_psec);
      tmp  = (float)(Ssec[i][l][m][k]/refDamp[k]);
      fwrite(&tmp,(size_t) sizeof(float), (size_t) 1, f_ssec);
    }
      }
      }
     }
    }
    fclose(f_psec);
    fclose(f_ssec);
}
/* ---------------------------------------------------------------------------------------------------------------------------------------------------------------- */



  /* Print header file */
  simpa.SimulationTime = time(0) - simpa.SimulationTime;
  sprintf(call,"./sims/%s.hed",simpa.dstring);
  if((fp1 = fopen(call,"w"))==NULL) printf("Error opening output file %s\n",call);
  else{

    fprintf(fp1,"#### SIMULATION ######################################\n");
    fprintf(fp1,"Simulation Tag: %s\n",simpa.tag);
    fprintf(fp1,"Number of Particles: %ld\n",simpa.npart);
    fprintf(fp1,"Number of Intervals in probability table: %ld\n", (long) simpa.q);
    fprintf(fp1,"Number of Intervals in (0 PI/2) for reflection/conversion tables: %d\n", simpa.Nangle);
    fprintf(fp1,"Grid size (x, y, z, theta, t): %d %d %d %d %d\n",simpa.GridSize[0],simpa.GridSize[1],simpa.GridSize[2],simpa.GridSize[3],simpa.GridSize[4]);
    fprintf(fp1,"Grid spacing (x, y, z, theta, t) in km, degree or s: %lf %lf %lf %lf %lf\n",simpa.GridSpacing[0], simpa.GridSpacing[1], simpa.GridSpacing[2], simpa.GridSpacing[3], simpa.GridSpacing[4]);
    fprintf(fp1,"Simulation time interval in s: %f\n",simpa.TimeStep);
    fprintf(fp1,"Simulation length in s: %f\n",simpa.SimLen);
    fprintf(fp1,"Reference 1/Q (internal amplification): %f\n",simpa.Qref);
    fprintf(fp1,"Energy modulo = %llu\n",simpa.emodulo);
    fprintf(fp1,"Execution time in seconds = %ld\n",(long)simpa.SimulationTime);
    fprintf(fp1,"Simulation time in seconds = %ld\n",(long)simpa.looptime);
    fprintf(fp1,"Output type (0 for 5D, 1 for section): %d\n",simpa.outputtype);

    fprintf(fp1,"\n");
    fprintf(fp1,"#### SOURCE PARAMETERS ###############################\n");
    fprintf(fp1,"Source position (x, y, z) in km from grid origin: %f %f %f\n",simpa.SourcePos[0], simpa.SourcePos[1], simpa.SourcePos[2]); 
    fprintf(fp1,"Source type (0 for P, 1 for S-Source): %d\n",simpa.SourceType);
    fprintf(fp1,"Source orientation (strike, dip, rake) in degrees: %f %f %f\n",simpa.SourceOri[0], simpa.SourceOri[1], simpa.SourceOri[2]); 
    fprintf(fp1,"Angular frequency: %f\n", simpa.frc);
    fprintf(fp1,"Source body: %d\n", simpa.sourcebody);

    fprintf(fp1,"\n");
    fprintf(fp1,"#### STRUCTURAL SETTINGS #############################\n");
    fprintf(fp1,"Number of bodies: %d\n",simpa.Nbody);
    for(i=0;i<simpa.Nbody;i++){
      fprintf(fp1,"---- body %d -----------------------------------------\n",i);
      fprintf(fp1,"X coordinate: %f %f\n",simpa.body[i].x[0], simpa.body[i].x[1]);
      fprintf(fp1,"Y coordinate: %f %f\n",simpa.body[i].y[0], simpa.body[i].y[1]);
      fprintf(fp1,"Z coordinate: %f %f\n",simpa.body[i].z[0], simpa.body[i].z[1]);
      fprintf(fp1,"Time Step fraction: %d\n",simpa.body[i].TimeStepFrac);
      fprintf(fp1,"S Velocity and gamma: %f %f\n",simpa.body[i].v[1], simpa.body[i].gam);
      fprintf(fp1,"Density: %f\n",simpa.body[i].rho);
      fprintf(fp1,"Type of ACF (gauss, expo, karman): %s\n",simpa.body[i].type);
      fprintf(fp1,"Kappa (only evaluated for karnam ACF): %f\n",simpa.body[i].kap);
      fprintf(fp1,"Correlation distance a: %f\n",simpa.body[i].a);
      fprintf(fp1,"Fractional fluctuation e: %f\n",simpa.body[i].e);
      fprintf(fp1,"Birch ny: %f\n",simpa.body[i].ny);
      fprintf(fp1,"Intrinsic Attenuation (1/Q_p 1/Q_s): %f %f\n",simpa.body[i].Qi[0], simpa.body[i].Qi[1]);
      fprintf(fp1,"Exponent of absobtion factor of single timestep: %f %f\n",simpa.body[i].Qex[0], simpa.body[i].Qex[1]);
      fprintf(fp1,"Inverse of free path of P- and S-waves: %f %f\n",simpa.body[i].pg[0], simpa.body[i].pg[1]);
      fprintf(fp1,"Mean free times  of P- and S-waves: %f %f\n",simpa.body[i].tau[0], simpa.body[i].tau[1]);
      fprintf(fp1,"Length of pathes in one time step for P- and S-waves: %f %f\n",simpa.body[i].path[0], simpa.body[i].path[1]);
      fprintf(fp1,"Total scattaring coefficients:\n");
      fprintf(fp1,"gpp0 = %e\ngps0 = %e\ngsp0 = %e\ngssl0 = %e\ngssr0 = %e\n",scattab[i].gpp0, scattab[i].gps0, scattab[i].gsp0, scattab[i].gssl0, scattab[i].gssr0);
      fprintf(fp1,"Forward-weighted scattaring coefficients as in Truner 1998 (BSSA):\n");
      fprintf(fp1,"<cos(theta)>_pp = %e\n<cos(theta)>_ps = %e\n<cos(theta)>_sp = %e\n<cos(theta)>_ssl = %e\n<cos(theta)>_ssr = %e\n", scattab[i].cos_pp, scattab[i].cos_ps, scattab[i].cos_sp, scattab[i].cos_ssl, scattab[i].cos_ssr);
      fprintf(fp1,"Transport mean free paths (Truner 1998, BSSA):\n");
      fprintf(fp1,"l^*_p = %e\n",scattab[i].l_p);
      fprintf(fp1,"l^*_s = %e\n",scattab[i].l_s);
    }
  }

  fclose(fp1);
  return 0;
}


double gammln(double xx)
{
  double x, y, tmp, ser;
  static double cof[6] = {76.18009172947146, -86.50532032941677, 24.01409824083091,
			  -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};
  int j;
  /*copied from numerical recipes page 214.*/
  y = x = xx;
  tmp = x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser = 1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}


int stokesRot(double phi, int mode, double IP, double I1,double I2, double U, double V, double * IP_p, double * I1_p,double * I2_p, double * U_p, double * V_p)
{
  /* Looking against the ray this function rotates counterclock wise phi according to Ishimaru[1978] p.35 */
  double tmpa[5];
  int retval = 0;
  if(mode == 1){       /* only necessarry for S-waves */
    tmpa[0] = pow(cos(phi),2);
    tmpa[1] = 1-tmpa[0];
    tmpa[2] = sin(2.*phi);
    *IP_p = IP;
    *I1_p = I1 * tmpa[0] + I2 * tmpa[1] + 0.5*tmpa[2] * U;
    *I2_p = I1 * tmpa[1] + I2 * tmpa[0] - 0.5*tmpa[2] * U;
    *U_p = -I1 * tmpa[2] + I2 * tmpa[2] + cos(2.*phi) * U;
    *V_p = V;
    if(*I1_p<0){     /* there are certain rotation angles where this might happen due to roundoff errors */
	printf("*I1_p = %e\n", *I1_p);
	*I2_p = I1 + I2;
	*I1_p = 0.;
	*U_p = 0.;
	retval = -1;
    }
    else if(*I2_p<0){
	printf("*I2_p = %e\n", *I2_p);
	*I1_p = I1 + I2;
	*I2_p = 0.;
	*U_p = 0.;
	retval = -1;
    }
  }
  else{
    *IP_p = IP;
    *I1_p = I1;
    *I2_p = I2;
    *U_p = U;
    *V_p = V;
  }
  return retval;
}



int checkSimulationParameters(simulationparameters * simpa_p)
{
  /* Declaration */
  double ls, lp;      /* fraction of simulqtion path length to mean (scattering) free path */
  int mc;             /* maximun number of counts from one particle in a cell in one time step */
  gridtype smc;       /* maximum number of counts in the source cell */
  gridtype mpc;       /* maximum possible number handled by gridtype */
  int mbit;           /* number of bytes of a variable of type gridtype */
  gridtype mweight;   /* maximum weigth of particle */
  int ntimeStep;      /* numbe rof simulation time steps */
  double MinQex = 1.; /* minimum Q exponent (or rather factor) */
  double MaxQex = 0;  /* maximum Q exponent (or rather factor) */
  double MinQ = 9999999; /* minimum Q */
  double MaxQ = 0;    /* maximum Q */
  double DiffQ;       /* differential Q */
  double amp;
  double RefQex;      /* reference Q exponent (or rather factor) */
  double TMinQex;     /* Total maximum attenuation factor */
  double TMaxQex;     /* Total maximum amplification factor */
  double pg_max[2];
  double epsilon_max;
  int i,j,k;

  /* Code */
  printf("***************************\n");
  DBG(CALL OF checkSimulationParameters);
  printf("***************************\n");
  printf("\n"); 
  /* This routine checks for reasonable parameters of the simulation */
  
  /* Find the maximum of epsilon over the whole domain (anomalious area and background). */
    epsilon_max=0;
    for (i=0; i< simpa_p->GridSize[0]; i++){
        for (j=0; j< simpa_p->GridSize[1]; j++){
            for (k=0; k< simpa_p->GridSize[2]; k++){
                if( simpa_p->body[0].ane[i][j][k]>epsilon_max){
                    epsilon_max=simpa_p->body[0].ane[i][j][k];
                }
            }
        }
    }
    if( simpa_p->body[0].e > epsilon_max){
        epsilon_max=simpa_p->body[0].e;
    }
    
  /* length of simulation time steps compared to mean free path and absorption length */
    for(i=0;i<simpa_p->Nbody;i++){
        
    /* Calculate pg[0] and pg[1] corresponding to that maximum */
    pg_max[0]=pow(epsilon_max,2) * simpa_p->body[i].pg[0];                              /* scattering coefficient from p */
    pg_max[1]=pow(epsilon_max,2) * simpa_p->body[i].pg[1];                              /* scattering coefficient from s */
      
    lp = 1./(simpa_p->body[i].path[0] * pg_max[0]);
    ls = 1./(simpa_p->body[i].path[1] * pg_max[1]);
    simpa_p->body[i].TimeStepFrac = 1;
    if(ls<lp){
      if(ls < 10){
        simpa_p->body[i].TimeStepFrac = (int)ceil(10./ls);
        simpa_p->body[i].path[0] = simpa_p->body[i].path[0]/simpa_p->body[i].TimeStepFrac;
        simpa_p->body[i].path[1] = simpa_p->body[i].path[1]/simpa_p->body[i].TimeStepFrac;
      }
    }
    else{
      if(lp < 10){
        simpa_p->body[i].TimeStepFrac = (int)ceil(10./lp);
        simpa_p->body[i].path[0] = simpa_p->body[i].path[0]/simpa_p->body[i].TimeStepFrac;
        simpa_p->body[i].path[1] = simpa_p->body[i].path[1]/simpa_p->body[i].TimeStepFrac;
      }
    }
    if(lp < 10.)printf("Warning: ratio of simulation path length and mean scattering free path of P-waves in body %d was only %4.3f.\nNow reduced %d\n", i, lp,simpa_p->body[i].TimeStepFrac);
    if(ls < 10.)printf("Warning: ratio of simulation path length and mean scattering free path of S-waves in body %d was only %4.3f.\nNow reduced %d\n", i, ls,simpa_p->body[i].TimeStepFrac);
  }

   /* printf("TimeStepFrac is : %d  and epsilon max is :%lf", simpa_p->body[0].TimeStepFrac, epsilon_max);*/
  /* maximum number of counts in one cell */
  mc = (int) ceil(simpa_p->GridSpacing[3]/simpa_p->TimeStep);

  /* maximum number of count in the source cell */
  smc = simpa_p->npart * mc;

  /* number of bits in gridtype */
  mbit = 8 * (int)sizeof(gridtype);
  mpc = (gridtype) pow(2,mbit)-1;

  /* maximum weight of particle */
  simpa_p->emodulo = mpc/smc;
  //simpa_p->emodulo = 100000;

  /* number of simulation time steps */
  ntimeStep = ceil((double)(simpa_p->GridSize[3] * simpa_p->GridSpacing[3])/simpa_p->TimeStep);

  /* minimum weight of particle at the end of the simulation */
  for(i=0;i<2;i++){  /* check only crust and mantle */
    MinQex = MIN(MinQex,exp(-simpa_p->frc * simpa_p->body[i].Qi[0] * simpa_p->SimLen));
    MinQex = MIN(MinQex,exp(-simpa_p->frc * simpa_p->body[i].Qi[1] * simpa_p->SimLen));
    MaxQex = MAX(MaxQex,exp(-simpa_p->frc * simpa_p->body[i].Qi[0] * simpa_p->SimLen));
    MaxQex = MAX(MaxQex,exp(-simpa_p->frc * simpa_p->body[i].Qi[1] * simpa_p->SimLen));
  }

  /* minimum and maximum Qi */
  for(i=0;i<2;i++){  /* check only crust and mantle */
    MinQ = MIN(MinQ,simpa_p->body[i].Qi[0]);
    MinQ = MIN(MinQ,simpa_p->body[i].Qi[1]);
    MaxQ = MAX(MaxQ,simpa_p->body[i].Qi[0]);
    MaxQ = MAX(MaxQ,simpa_p->body[i].Qi[1]);
  }
  DiffQ = MaxQ-MinQ;
  amp = exp(simpa_p->frc * DiffQ * simpa_p->SimLen);
  if(simpa_p->emodulo * exp(-simpa_p->frc * MaxQ * simpa_p->SimLen) < 1){  /* Internal amplification is nesseccary because otherwise the particle energy might be zero before end of simulation */
    //simpa_p->Qref = MaxQ - log(simpa_p->emodulo)/simpa_p->frc/simpa_p->SimLen;  /* adjust attenuation such that maximal attenuated particle has energy 1 at the end of the simulation */
    simpa_p->Qref = MinQ;         /* minimaly attenuated particle remains with constant energy */
    RefQex = exp(simpa_p->frc * (simpa_p->Qref-MinQ) * simpa_p->SimLen);
    if(RefQex>smc){printf("Warning: due to compensation for inrinsic attenuation the particles may overflow.\n");}
  }
  /*simpa_p->Qref = 0;*/
  if(amp>mpc){
    printf("Warning: the difference in attenuation is larger than the dynamic range of the recording.\n");
  }
  return 0;
}


int whichbody(particle * pat_p, simulationparameters * simpa_p)
{
  /* checks particle position and body geomety and returns 1 if the particle is in body 2 (only above the interface). */
#if BODY == 1
  /* cylinder ring around the source */
  double rsq;  /* squared distance from source */
  rsq = pow(pat_p->pos[0]-simpa_p->SourcePos[0],2) + pow(pat_p->pos[1]-simpa_p->SourcePos[1],2);
  if (rsq>simpa_p->body[2].x[0] & rsq<simpa_p->body[2].x[1]){
#elif BODY==2
  /* box body */
  if ((pat_p->pos[0]>simpa_p->body[2].x[0]) & (pat_p->pos[0]<simpa_p->body[2].x[1]) & (pat_p->pos[1]>simpa_p->body[2].y[0]) & (pat_p->pos[1]<simpa_p->body[2].y[1])){
#else
  if (0){
#endif
    return 1;
  }
  else{
    return 0;
  }
}
