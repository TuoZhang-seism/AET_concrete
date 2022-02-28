
/***************************************************************************
 *   Copyright (C) 2006 by Christoph Sens-Schönfelder                      *
 *   sens-schoenfelder@uni-leipzig.de                                      *
 *                                                                         *
 ***************************************************************************/


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "edefinition.h"



int read_sim_head(simulationparameters *, bodytable *);
int read_ref_trace(char *, double **, double **, int*);
float ** read_sec(simulationparameters *, char*);
int print_tmp_sim_parms(simulationparameters, simulationparameters, double);

int main(int ARGC, char *ARGV[])
{
  simulationparameters simpaP, simpaS;       /* stores all parameters that define the simulation */
  bodytable * scattab;
  int status;
  char ref_trace_name[] = "westPyr.asc";
  double * trc, * slow;
  int N;
  float ** PP, ** PS, ** SP, ** SS;          /* data of seismic section */
  double * reftr;                            /* sesimogram that is to be compared to the data */
  double * avtr;                             /* simulated trace avaraged over a certain distance range (slownes stack) */
  int * avcnt;                               /* counter for stacking */
  double mf = 0;                             /* misfit */
  double meanmf = 0;                         /* mean misfit */
  int cnt = 0;                               /* counter */  
  double dt = 1.;
  float distrange[2] = {600, 900};
  float slowrange[2] = {0.124, 0.47};
  double refdist;
  double gridspacing;
  int i,j;
  int index;
  double spart = 0.95;
  FILE * fp;
  int twstart, twend;
  double pmf, smf;
  double refint, avint;

  refdist = (distrange[0] + distrange[1])/2;   /* reference distance for the stack */
  if(ARGC>2){
    strcpy(simpaP.ParFileName,ARGV[1]);
    strcpy(simpaS.ParFileName,ARGV[2]);
  }
  else{
    printf("%s:%d Please give header file name of P-source as first\nand header file name of S-source as second  argument.\n",__FILE__,__LINE__);
    return -1;
  }
  /* read header */
  status = read_sim_head(&simpaP, scattab);
  status = read_sim_head(&simpaS, scattab);
  /* read reference trace */
  status = read_ref_trace(ref_trace_name,&trc,&slow,&N);
  /* read data */
  PP = read_sec(&simpaP,"P");
  PS = read_sec(&simpaP,"S");
  SP = read_sec(&simpaS,"P");
  SS = read_sec(&simpaS,"S");
  gridspacing = (simpaP.GridSpacing[0] + simpaP.GridSpacing[1])/2;
  reftr = (double*) calloc((size_t) simpaP.GridSize[3],sizeof(double));
  avtr = (double*) calloc((size_t) simpaP.GridSize[3],sizeof(double));
  avcnt = (int *) calloc((size_t) simpaP.GridSize[3],sizeof(int));
  /* stack data of different distances */
  for(i=0;i<simpaP.GridSize[0];i++){
    if((i*gridspacing > distrange[0]) & (i*gridspacing < distrange[1])){
      for(j=0;j<simpaP.GridSize[3];j++){
        index = (int)round((double)j*simpaP.GridSpacing[3]/(i*gridspacing)/(dt/refdist));
        if(index<simpaP.GridSize[3]){
          avtr[index] += (double)((PP[i][j] + PS[i][j]) * (1.-spart) + (SP[i][j] + SS[i][j]) * spart); 
          avcnt[index] ++;
        }
      }
    }
  }
  for(j=0;j<simpaP.GridSize[3];j++){
    avtr[j] = sqrt(avtr[j]/avcnt[j]);  /* take sqrt of energy to compare with amplitude */
    avcnt[j] = 0;
  }
  /* resample reference trace */
  for(j=0;j<N;j++){
    index = (int)round((slow[j]-0.001)/(dt/refdist));  /* slightly correct the slowness since the data trace is too late */
    reftr[index] += (double)MAX(trc[j],0.);
    avcnt[index] ++;
  }
  for(j=0;j<simpaP.GridSize[3];j++){
    reftr[j] /= avcnt[j];
  }
  /* estimate constant offsett factor */
  twstart = (int)ceil(slowrange[0]/(simpaP.GridSpacing[3]/refdist));
  twend = (int)floor(MIN(slowrange[1]/(simpaP.GridSpacing[3]/refdist),index));  /* MIN to be sure that the reference trace is long enough */
  for(j=twstart;j<=twend;j++){
    if(avtr[j] != 0){  /* evaluate only in this slownes window */ 
      cnt ++;
      meanmf += log(reftr[j]) - log(avtr[j]);
    }
  }
  meanmf /= cnt;

  /* estimate misfit from offset corrected traces */
  fp = fopen("ref_av_trace","w");
  for(j=twstart;j<=twend;j++){
    if(avtr[j] != 0){  /* evaluate misfit only in this slownes window */ 
      mf += fabs(log(reftr[j]) - log(avtr[j]) - meanmf);
      fprintf(fp,"%d %e %e %e\n",j,reftr[j],avtr[j],fabs(log(reftr[j]) - log(avtr[j]) - meanmf));
    }
  }
  fclose(fp);
  mf /= cnt;

  /* estimate cumulative misfit around P-arrival */
  twstart = (int)floor(0.162/(simpaP.GridSpacing[3]/refdist));  /* the 0.003s/km slownes correction is in here and the end is 5% later */
  twend = (int)floor(0.17/(simpaP.GridSpacing[3]/refdist));  /* MIN to be sure that the reference trace is long enough */
  refint = 0;
  avint = 0;
  for(j=twstart;j<=twend;j++){
    refint += reftr[j];
    avint += avtr[j];
  }
  pmf = fabs(log(refint) - log(avint) - meanmf);


  /* estimate cumulative misfit around S-arrival */
  twstart = (int)floor(0.278/(simpaP.GridSpacing[3]/refdist));  /* the 0.003s/km slownes correction is in here and the end is 15% later */
  twend = (int)floor(0.32/(simpaP.GridSpacing[3]/refdist));  /* MIN to be sure that the reference trace is long enough */
  refint = 0;
  avint = 0;
  for(j=twstart;j<=twend;j++){
    refint += reftr[j];
    avint += avtr[j];
  }
  smf = fabs(log(refint) - log(avint) - meanmf);

  /* Average the misfits */
  //mf = (mf + (pmf + smf)/2.)/2.;  /* Average of total misfit and integral misfit in P- and S-windows





  if(isnan(mf))mf = 111.;  /*avoid nans */
  if(cnt<0.95*(twend-twstart))mf = 112.;
  status = print_tmp_sim_parms(simpaP,simpaS,mf);

  return 0;
}




int print_tmp_sim_parms(simulationparameters simpaP, simulationparameters simpaS, double mf){
  FILE * fp;
  char * snameendpos;
  snameendpos = index(simpaP.tag,'q');
  sprintf(snameendpos,"\0");
  fp = fopen("tmp_sim_parms","w");
  fprintf(fp,"%s\n",simpaP.tag);
  fprintf(fp,"%s_P.sec\n",simpaP.ParFileName);
  fprintf(fp,"%s_S.sec\n",simpaP.ParFileName);
  fprintf(fp,"%s_P.sec\n",simpaS.ParFileName);
  fprintf(fp,"%s_S.sec\n",simpaS.ParFileName);
  fprintf(fp,"%12.11e\n",simpaP.body[0].e);
  fprintf(fp,"%12.11e\n",simpaP.body[1].e);
  fprintf(fp,"%12.11e\n",simpaP.body[2].e);

  fprintf(fp,"%12.11e\n",simpaP.body[0].a);
  fprintf(fp,"%12.11e\n",simpaP.body[1].a);
  fprintf(fp,"%12.11e\n",simpaP.body[2].a);

  fprintf(fp,"%12.11e\n",simpaP.body[0].Qi[0]);
  fprintf(fp,"%12.11e\n",simpaP.body[0].Qi[1]);

  fprintf(fp,"%12.11e\n",simpaP.body[1].Qi[0]);
  fprintf(fp,"%12.11e\n",simpaP.body[1].Qi[1]);

  fprintf(fp,"%12.11e\n",simpaP.body[2].Qi[0]);
  fprintf(fp,"%12.11e\n",simpaP.body[2].Qi[1]);

  fprintf(fp,"%12.11e\n",mf);
  fclose(fp);
  return 0;
}


int read_sim_head(simulationparameters * simpa_p, bodytable * scattab){
  /* Declaration */
  FILE * parfile_p;
  int count = 0;
  int i;
  char dump[101];
  char fname[201];
  int status;
  /* Code */
  /* DBG(CALL OF ReadParFile); */
  status = sprintf(fname,"%s.hed",simpa_p->ParFileName);
  if((parfile_p = fopen(fname,"r"))==NULL){
    return -1;
  }
  fgets(dump,100,parfile_p);
  count = fscanf(parfile_p,"Simulation Tag: %s\n",simpa_p->tag);
  count += fscanf(parfile_p,"Number of Particles: %ld\n",&simpa_p->npart);
  count += fscanf(parfile_p,"Number of Intervals in probability table: %lf\n", &simpa_p->q);
  count += fscanf(parfile_p,"Number of Intervals in (0 PI/2) for reflection/conversion tables: %d\n",&simpa_p->Nangle);
  count += fscanf(parfile_p,"Grid size (x, y, z, t): %d %d %d %d\n",&simpa_p->GridSize[0],&simpa_p->GridSize[1],&simpa_p->GridSize[2],&simpa_p->GridSize[3]);
  count += fscanf(parfile_p,"Grid spacing (x, y, z, t) in km or s: %f %f %f %f\n",&simpa_p->GridSpacing[0], &simpa_p->GridSpacing[1], &simpa_p->GridSpacing[2], &simpa_p->GridSpacing[3]);
  count += fscanf(parfile_p,"Simulation time interval in s: %f\n",&simpa_p->TimeStep);
  count += fscanf(parfile_p,"Simulation length in s: %d\n",&simpa_p->SimLen);
  count += fscanf(parfile_p,"Reference 1/Q (internal amplification): %f\n",&simpa_p->Qref);
  count += fscanf(parfile_p,"Energy modulo = %d\n",&simpa_p->emodulo);
  count += fscanf(parfile_p,"Execution time in seconds = %d\n",&simpa_p->SimulationTime);
  count += fscanf(parfile_p,"Simulation time in seconds = %d\n",&simpa_p->looptime);
  count += fscanf(parfile_p,"Output type (0 for 4D, 1 for section): %d\n",&simpa_p->outputtype);

  fgets(dump,100,parfile_p);
  count += fscanf(parfile_p,"Source position (x, y, z) in km from grid origin: %lf %lf %lf\n",&simpa_p->SourcePos[0], &simpa_p->SourcePos[1], &simpa_p->SourcePos[2]);
  count += fscanf(parfile_p,"Source type (0 for P, 1 for S-Source): %d\n",&simpa_p->SourceType);
  count += fscanf(parfile_p,"Source orientation (strike, dip, rake) in degrees: %lf %lf %lf\n",&simpa_p->SourceOri[0], &simpa_p->SourceOri[1], &simpa_p->SourceOri[2]); 
  count += fscanf(parfile_p,"Angular frequency: %lf\n", &simpa_p->frc);
  count += fscanf(parfile_p,"Source body: %d\n", &simpa_p->sourcebody);

  fgets(dump,100,parfile_p);
  count += fscanf(parfile_p,"Number of bodies: %d\n",&simpa_p->Nbody);
  simpa_p->body = (bodyparm *) calloc((size_t) simpa_p->Nbody, sizeof(bodyparm));
  scattab = (bodytable *) malloc((size_t) simpa_p->Nbody * sizeof(bodytable));
  for(i=0;i<simpa_p->Nbody;i++){
    fgets(dump,100,parfile_p);
    count += fscanf(parfile_p,"X coordinate: %lf %lf\n",&simpa_p->body[i].x[0], &simpa_p->body[i].x[1]);
    count += fscanf(parfile_p,"Y coordinate: %lf %lf\n",&simpa_p->body[i].y[0], &simpa_p->body[i].y[1]);
    count += fscanf(parfile_p,"Z coordinate: %lf %lf\n",&simpa_p->body[i].z[0], &simpa_p->body[i].z[1]);
    count += fscanf(parfile_p,"Time Step fraction: %d\n",&simpa_p->body[i].TimeStepFrac);

    count += fscanf(parfile_p,"S Velocity and gamma: %lf %lf\n",&simpa_p->body[i].v[1], &simpa_p->body[i].gam);
    count += fscanf(parfile_p,"Density: %lf\n",&simpa_p->body[i].rho);
    count += fscanf(parfile_p,"Type of ACF (gauss, expo, karman): %s\n",simpa_p->body[i].type);
    count += fscanf(parfile_p,"Kappa (only evaluated for karnam ACF): %lf\n",&simpa_p->body[i].kap);
    count += fscanf(parfile_p,"Correlation distance a: %lf\n",&simpa_p->body[i].a);
    count += fscanf(parfile_p,"Fractional fluctuation e: %lf\n",&simpa_p->body[i].e);
    count += fscanf(parfile_p,"Birch ny: %lf\n",&simpa_p->body[i].ny);
    count += fscanf(parfile_p,"Intrinsic Attenuation (1/Q_p 1/Q_s): %lf %lf\n",&simpa_p->body[i].Qi[0], &simpa_p->body[i].Qi[1]);
    count += fscanf(parfile_p,"Exponent of absobtion factor of single timestep: %f %f\n",&simpa_p->body[i].Qex[0], &simpa_p->body[i].Qex[1]);
    count += fscanf(parfile_p,"Inverse of free path of P- and S-waves: %f %f\n",&simpa_p->body[i].pg[0], &simpa_p->body[i].pg[1]);
    count += fscanf(parfile_p,"Mean free times  of P- and S-waves: %f %f\n",&simpa_p->body[i].tau[0], &simpa_p->body[i].tau[1]);
    count += fscanf(parfile_p,"Length of pathes in one time step for P- and S-waves: %f %f\n",&simpa_p->body[i].path[0], &simpa_p->body[i].path[1]);
    fgets(dump,100,parfile_p);
    count += fscanf(parfile_p,"gpp0 = %e\ngps0 = %e\ngsp0 = %e\ngssl0 = %e\ngssr0 = %e\n",&scattab[i].gpp0, &scattab[i].gps0, &scattab[i].gsp0, &scattab[i].gssl0, &scattab[i].gssr0);
    fgets(dump,100,parfile_p);
    count += fscanf(parfile_p,"<cos(theta)>_pp = %e\n<cos(theta)>_ps = %e\n<cos(theta)>_sp = %e\n<cos(theta)>_ssl = %e\n<cos(theta)>_ssr = %e\n", &scattab[i].cos_pp, &scattab[i].cos_ps, &scattab[i].cos_sp, &scattab[i].cos_ssl, &scattab[i].cos_ssr);
    fgets(dump,100,parfile_p);
    count += fscanf(parfile_p,"l^*_p = %e\n",&scattab[i].l_p);
    count += fscanf(parfile_p,"l^*_s = %e\n",&scattab[i].l_s);
    simpa_p->body[i].v[0] = simpa_p->body[i].v[1] * simpa_p->body[i].gam;
    simpa_p->body[i].l = simpa_p->frc / simpa_p->body[i].v[1];
    simpa_p->body[i].q = simpa_p->q;
    simpa_p->body[i].path[0] = simpa_p->body[i].v[0] * simpa_p->TimeStep;
    simpa_p->body[i].path[1] = simpa_p->body[i].v[1] * simpa_p->TimeStep;
    if(!strcmp(simpa_p->body[i].type,"karman")){
      sprintf(simpa_p->body[i].ProbTabFileName,"tables/%s_%1.2f_%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_%1.2e_scat",simpa_p->body[i].type,simpa_p->body[i].kap,simpa_p->body[i].gam,simpa_p->body[i].l, simpa_p->body[i].a, simpa_p->body[i].e, simpa_p->body[i].ny, simpa_p->q);
    }
    else{
      sprintf(simpa_p->body[i].ProbTabFileName,"tables/%s_%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_%1.2e_scat",simpa_p->body[i].type,simpa_p->body[i].gam,simpa_p->body[i].l, simpa_p->body[i].a, simpa_p->body[i].e, simpa_p->body[i].ny, simpa_p->q);
    }
  }
  fclose(parfile_p);
  return 0;
}


int read_ref_trace(char * ref_trace_name, double ** trc_p, double ** tim_p, int * N){
  FILE * dfile;
  int lcount = 0;
  char dump[101];
  int i;
  if((dfile = fopen(ref_trace_name,"r"))==NULL){
    return -1;
  }
  while(!feof(dfile)){
    fgets(dump,100,dfile);
    lcount ++;
  }
  lcount --;
  rewind(dfile);
  *trc_p = (double *) malloc((size_t)lcount * sizeof(double));
  *tim_p = (double *) malloc((size_t)lcount * sizeof(double));

  for(i=0;i<lcount;i++){
    fscanf(dfile," %lf %lf\n",*tim_p+i, *trc_p+i);
  }
  *N = lcount;
  return 0;
}


float ** read_sec(simulationparameters * simpa, char * comp){
  FILE * Pfile;
  int GridSize;
  char fname[201];
  int status, i;
  float ** P;
  status = sprintf(fname,"%s_%s.sec",simpa->ParFileName,comp);
   if((Pfile = fopen(fname,"r"))==NULL){
    return P;
  }

  /* determin filesize */
  status = fseek(Pfile, 0, SEEK_END);
  status = ftell(Pfile);
  /* using 32bit floats results in */
  simpa->GridSize[0] = status/4 / simpa->GridSize[3];
  rewind(Pfile);

  /* Allocate arrays */
  P = (float **)malloc((size_t) simpa->GridSize[0] * sizeof(float*));
  for(i=0;i<simpa->GridSize[0];i++){
    P[i] = (float *)calloc((size_t)simpa->GridSize[3],sizeof(float));
    status = fread(P[i], (size_t) sizeof(float), (size_t) simpa->GridSize[3], Pfile);
  }
  fclose(Pfile);
  return P;
}
