#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "edefinition.h"
#include "scat_coeff_2D.h"
#include "mk_source_pattern.h"
  /* #include "collect.h"*/


int MkSourcePattern(simulationparameters simpa, float ** fltpt)
{
  char pr_file_name[] = "./tables/p_rad_tab.dat";
  char st_file_name[] = "./tables/st_rad_tab.dat";
  char sf_file_name[] = "./tables/sf_rad_tab.dat";
  char cos_phi_file_name[] = "./tables/cos2phi_tab.dat";
  char sin_phi_file_name[] = "./tables/sin2phi_tab.dat";

  FILE *pr_file, *st_file, *sf_file, *cos_phi_file, *sin_phi_file;
  float *pr_rad_t, *st_rad_t, *sf_rad_t, *cos_phi_rad, *sin_phi_rad;
  int i, tmpi, pr_num=0, st_num=0, sf_num=0, phi_num=0;
  float *r, *r2;
  double *pr_b, *st_b, *sf_b, *phi_b, tmp;   /* values of pattern */
  float *pr_r2, *st_r2, *sf_r2, *phi_r2;
  double *pr_b_save, *st_b_save, *sf_b_save, *phi_b_save;   /* backups of pattern and random numbers */
  double pr_max=0, pr_sum=0, st_max=0, st_sum=0, sf_max=0, sf_sum=0, phi_max=0, phi_sum=0;
  double pr_max_save=0, st_max_save=0, sf_max_save=0, phi_max_save=0;

  pr_rad_t = (float *) calloc((size_t) simpa.q, (size_t) sizeof(float));
  st_rad_t = (float *) calloc((size_t) simpa.q, (size_t) sizeof(float));
  sf_rad_t = (float *) calloc((size_t) simpa.q, (size_t) sizeof(float));
  cos_phi_rad = (float *) calloc((size_t) simpa.q, (size_t) sizeof(float));
  sin_phi_rad = (float *) calloc((size_t) simpa.q, (size_t) sizeof(float));

  fltpt[0] = pr_rad_t;
  fltpt[1] = st_rad_t;
  fltpt[2] = sf_rad_t;
  fltpt[3] = cos_phi_rad;
  fltpt[4] = sin_phi_rad;

  pr_b = (double *) calloc((size_t) simpa.q, (size_t) sizeof(double));
  st_b = (double *) calloc((size_t) simpa.q, (size_t) sizeof(double));
  sf_b = (double *) calloc((size_t) simpa.q, (size_t) sizeof(double));
  phi_b = (double *) calloc((size_t) simpa.q, (size_t) sizeof(double));
  pr_b_save = (double *) calloc((size_t) simpa.q, (size_t) sizeof(double));
  st_b_save = (double *) calloc((size_t) simpa.q, (size_t) sizeof(double));
  sf_b_save = (double *) calloc((size_t) simpa.q, (size_t) sizeof(double));
  phi_b_save = (double *) calloc((size_t) simpa.q, (size_t) sizeof(double));
  pr_r2 = (float *) calloc((size_t) simpa.q, (size_t) sizeof(float));
  st_r2 = (float *) calloc((size_t) simpa.q, (size_t) sizeof(float));
  sf_r2 = (float *) calloc((size_t) simpa.q, (size_t) sizeof(float));
  phi_r2 = (float *) calloc((size_t) simpa.q, (size_t) sizeof(float));
  r = (float *) calloc((size_t) simpa.q, (size_t) sizeof(float));
  r2 = (float *) calloc((size_t) simpa.q, (size_t) sizeof(float));

  while((pr_num<simpa.q) | (st_num<simpa.q) | (sf_num<simpa.q) | (phi_num<simpa.q)){
    for(i=0;i<simpa.q;i++){           /* scan pattern */
      r[i] = PI*rand()/RAND_MAX;
      tmp = sin(r[i]);        /* differential of spherical coordinates */
      pr_b[i] = BPR_THETA(r[i])*tmp;
      st_b[i] = BST_THETA(r[i])*tmp;
      sf_b[i] = BSF_THETA(r[i])*tmp;
      phi_b[i] = BST_PHI(r[i] * 2.);    /* PHI part of P-wave pattern (* 2. because PHI = [0,2PI]) */
      r2[i] = ((float)rand())/RAND_MAX;
    }
    for(i=0;i<simpa.q;i++){          /* find maximum and sum */
      if(pr_b[i]>pr_max)pr_max=pr_b[i];
      if(st_b[i]>st_max)st_max=st_b[i];
      if(sf_b[i]>sf_max)sf_max=sf_b[i];
      if(phi_b[i]>phi_max)phi_max=phi_b[i];
      pr_sum += pr_b[i];
      st_sum += st_b[i];
      sf_sum += sf_b[i];
      phi_sum += phi_b[i];
    }
    if(pr_max != pr_max_save){    /* in case maximum has changed, i.e. probability distribution is better sampled */
      tmpi = pr_num;
      pr_num = 0;
      pr_max_save = pr_max;
      for(i=0;i<tmpi;i++){
	if(pr_b_save[i]/pr_max/pr_r2[i]>=1){
	  pr_rad_t[pr_num] = pr_rad_t[i];
	  pr_b_save[pr_num] = pr_b_save[i];
	  pr_r2[pr_num] = pr_r2[i];
	  pr_num ++;
	}
      }
    }
    if(st_max != st_max_save){    /* in case maximum has changed, i.e. probability distribution is better sampled */
      tmpi = st_num;
      st_num = 0;
      st_max_save = st_max;
      for(i=0;i<tmpi;i++){
	if(st_b_save[i]/st_max/st_r2[i]>=1){
	  st_rad_t[st_num] = st_rad_t[i];
	  st_b_save[st_num] = st_b_save[i];
	  st_r2[st_num] = st_r2[i];
	  st_num ++;
	}
      }
    }
    if(sf_max != sf_max_save){    /* in case maximum has changed, i.e. probability distribution is better sampled */
      tmpi = sf_num;
      sf_num = 0;
      sf_max_save = sf_max;
      for(i=0;i<tmpi;i++){
	if(sf_b_save[i]/sf_max/sf_r2[i]>=1){
	  sf_rad_t[sf_num] = sf_rad_t[i];
	  sf_b_save[sf_num] =sf_b_save[i];
	  sf_r2[sf_num] = sf_r2[i];
	  sf_num ++;
	}
      }
    }
    if(phi_max != phi_max_save){    /* in case maximum has changed, i.e. probability distribution is better sampled */
      tmpi = phi_num;
      phi_num = 0;
      phi_max_save = phi_max;
      for(i=0;i<tmpi;i++){
	if(phi_b_save[i]/phi_max/phi_r2[i]>=1){
	  sin_phi_rad[phi_num] = sin_phi_rad[i];
	  phi_b_save[phi_num] = phi_b_save[i];
	  phi_r2[phi_num] = phi_r2[i];
	  phi_num ++;
	}
      }
    }
    for(i=0;i<simpa.q;i++){          /* find new directions */
      if((pr_b[i]/pr_max/r2[i]>=1) & (pr_num<simpa.q)){
	pr_rad_t[pr_num] = r[i];
	pr_b_save[pr_num] = pr_b[i];
	pr_r2[pr_num] = r2[i];
	pr_num ++;
      }
      if((st_b[i]/st_max/r2[i]>=1) & (st_num<simpa.q)){
	st_rad_t[st_num] =  r[i];
	st_b_save[st_num] = st_b[i];
	st_r2[st_num] = r2[i];
	st_num ++;
      }
      if((sf_b[i]/sf_max/r2[i]>=1) & (sf_num<simpa.q)){
	sf_rad_t[sf_num] = r[i];
	sf_b_save[sf_num] = sf_b[i];
	sf_r2[sf_num] = r2[i];
	sf_num ++;
      }
      if((phi_b[i]/phi_max/r2[i]>=1) & (phi_num<simpa.q)){
	sin_phi_rad[phi_num] = r[i] * 2.;
	phi_b_save[pr_num] = phi_b[i];
	phi_r2[pr_num] = r2[i];
	phi_num ++;
      }
    }
  }
  for(i=0;i<simpa.q;i++){          /* calculate cos_phi_rad */
    cos_phi_rad[i] = sin_phi_rad[i]+PI4;
    if(cos_phi_rad[i]>P2I){
      cos_phi_rad[i] -= P2I;
    }
  }

  /* Write probability tables to files */
  if((pr_file = fopen(pr_file_name,"w"))==NULL){
    printf("%s:%d Error creating file %s\n",__FILE__,__LINE__,pr_file_name);
    return -1;
  }
  if((st_file = fopen(st_file_name,"w"))==NULL){
    printf("%s:%d Error creating file %s\n",__FILE__,__LINE__,st_file_name);
    return -1;
  }
  if((sf_file = fopen(sf_file_name,"w"))==NULL){
    printf("%s:%d Error creating file %s\n",__FILE__,__LINE__,sf_file_name);
    return -1;
  }
  if((cos_phi_file = fopen(cos_phi_file_name,"w"))==NULL){
    printf("%s:%d Error creating file %s\n",__FILE__,__LINE__,cos_phi_file_name);
    return -1;
  }
  if((sin_phi_file = fopen(sin_phi_file_name,"w"))==NULL){
    printf("%s:%d Error creating file %s\n",__FILE__,__LINE__,sin_phi_file_name);
    return -1;
  }
  fprintf(pr_file,"Probability table for THETA angle of P-wave radiation of a double couple source with fault normal (1,0,0) and slip in (0,1,0) direction.\n");
  fprintf(st_file,"Probability table for THETA angle of S-wave radiation polarized in THETA direction of a double couple source with fault normal (1,0,0) and slip in (0,1,0) direction.\n");
  fprintf(sf_file,"Probability table for THETA angle of S-wave radiation polarized in PHI direction of a double couple source with fault normal (1,0,0) and slip in (0,1,0) direction.\n");
  fprintf(cos_phi_file,"Probability table of (cos(2*phi))^2\n");
  fprintf(sin_phi_file,"Probability table of (sin(2*phi))^2\n");

  for(i=0;i<simpa.q;i++){
    fprintf(pr_file,"%1.12f\n",pr_rad_t[i]);
    fprintf(st_file,"%1.12f\n",st_rad_t[i]);
    fprintf(sf_file,"%1.12f\n",sf_rad_t[i]);
    fprintf(cos_phi_file,"%1.12f\n",cos_phi_rad[i]);
    fprintf(sin_phi_file,"%1.12f\n",sin_phi_rad[i]);
  }
  fclose(pr_file);
  fclose(st_file);
  fclose(sf_file);
  fclose(cos_phi_file);
  fclose(sin_phi_file);
  *(fltpt+0) = pr_rad_t;
  *(fltpt+1) = st_rad_t;
  *(fltpt+2) = sf_rad_t;
  *(fltpt+3) = cos_phi_rad;
  *(fltpt+4) = sin_phi_rad;
  return 0;
}


int ReadRadPat(simulationparameters simpa, float ** fltpt)
{
  char pr_file_name[] = "./tables/p_rad_tab.dat";
  char st_file_name[] = "./tables/st_rad_tab.dat";
  char sf_file_name[] = "./tables/sf_rad_tab.dat";
  char cos_phi_file_name[] = "./tables/cos2phi_tab.dat";
  char sin_phi_file_name[] = "./tables/sin2phi_tab.dat";
  FILE * pr_file, * st_file, * sf_file, *cos_phi_file, *sin_phi_file;
  float *pr_rad_t, *st_rad_t, *sf_rad_t, *cos_phi_rad, *sin_phi_rad;
  char dump[501];
  int count = 0, i, tmp;
  int pr_num=0, st_num=0, sf_num=0, cos_num=0, sin_num=0;

  pr_rad_t = (float *) calloc((size_t) simpa.q, (size_t) sizeof(float));
  st_rad_t = (float *) calloc((size_t) simpa.q, (size_t) sizeof(float));
  sf_rad_t = (float *) calloc((size_t) simpa.q, (size_t) sizeof(float));
  cos_phi_rad = (float *) calloc((size_t) simpa.q, (size_t) sizeof(float));
  sin_phi_rad = (float *) calloc((size_t) simpa.q, (size_t) sizeof(float));

  fltpt[0] = pr_rad_t;
  fltpt[1] = st_rad_t;
  fltpt[2] = sf_rad_t;
  fltpt[3] = cos_phi_rad;
  fltpt[4] = sin_phi_rad;

  /*###################################*/
  if((pr_file = fopen(pr_file_name,"r"))==NULL){
    printf("%s:%d Error opening file %s\n",__FILE__,__LINE__,pr_file_name);
    return -1;
  }
  while(!feof(pr_file)){
    fgets(dump,500,pr_file);
    count ++;
  }
  if(count-1!=simpa.q)return -1;
  rewind(pr_file);
  fgets(dump,500,pr_file);
  for(i=0;i<count;i++){
    tmp = fscanf(pr_file,"%f", &pr_rad_t[i]);
    if(tmp == 1)pr_num ++;
  }
  fclose(pr_file);
  if(pr_num!=simpa.q)return -1;
  /*###################################*/
  count = 0;
  if((st_file = fopen(st_file_name,"r"))==NULL){
    printf("%s:%d Error opening file %s\n",__FILE__,__LINE__,st_file_name);
    return -1;
  }
  while(!feof(st_file)){
    fgets(dump,500,st_file);
    count ++;
  }
  if(count-1!=simpa.q)return -1;
  rewind(st_file);
  fgets(dump,500,st_file);
  for(i=0;i<count-1;i++){
   tmp = fscanf(st_file,"%f",&st_rad_t[i]);
   if(tmp == 1)st_num ++;
  }
  fclose(st_file);
  if(st_num!=simpa.q)return -1;
  /*###################################*/
  count = 0;
  if((sf_file = fopen(sf_file_name,"r"))==NULL){
    printf("%s:%d Error opening file %s\n",__FILE__,__LINE__,sf_file_name);
    return -1;
  }
  while(!feof(sf_file)){
    fgets(dump,500,sf_file);
    count ++;
  }
  if(count-1!=simpa.q)return -1;
  rewind(sf_file);
  fgets(dump,500,sf_file);
  for(i=0;i<count-1;i++){
    tmp = fscanf(sf_file,"%f",&sf_rad_t[i]);
    if(tmp == 1)sf_num ++;
  }
  fclose(sf_file);
  if(sf_num!=simpa.q)return -1;
  /*###################################*/
  count = 0;
  if((cos_phi_file = fopen(cos_phi_file_name,"r"))==NULL){
    printf("%s:%d Error opening file %s\n",__FILE__,__LINE__,cos_phi_file_name);
    return -1;
  }
  while(!feof(cos_phi_file)){
    fgets(dump,500,cos_phi_file);
    count ++;
  }
  if(count-1!=simpa.q)return -1;
  rewind(cos_phi_file);
  fgets(dump,500,cos_phi_file);
  for(i=0;i<count-1;i++){
    if(fscanf(cos_phi_file,"%f",&cos_phi_rad[i])==1)cos_num++;
  }
  fclose(cos_phi_file);
  if(cos_num!=simpa.q)return -1;
  /*###################################*/
  count = 0;
  if((sin_phi_file = fopen(sin_phi_file_name,"r"))==NULL){
    printf("%s:%d Error opening file %s\n",__FILE__,__LINE__,sin_phi_file_name);
    return -1;
  }
  while(!feof(sin_phi_file)){
    fgets(dump,500,sin_phi_file);
    count ++;
  }
  if(count-1!=simpa.q)return -1;
  rewind(sin_phi_file);
  fgets(dump,500,sin_phi_file);
  for(i=0;i<count-1;i++){
    if(fscanf(sin_phi_file,"%f",&sin_phi_rad[i])==1)sin_num++;;
  }
  fclose(sin_phi_file);
  if(sin_num!=simpa.q)return -1;
  /*###################################*/
  return 0;
}


int source_rot_mat(simulationparameters simpa, double * l, double * m, double * n)
{  
  double c[4], s[4];
  int i;
/* Rotate directions according to source orientation */
  for(i=1;i<4;i++){
    c[i]=cos(simpa.SourceOri[i-1]*PI/180);
    s[i]=sin(simpa.SourceOri[i-1]*PI/180);
  }
  /* COMBINED ROTATION */
  l[0] = c[3]*c[1]+s[3]*c[2]*s[1];
  l[1] = -c[3]*s[1]+s[3]*c[2]*c[1];
  l[2] = -s[3]*s[2];
  m[0] = s[2]*s[1];
  m[1] = s[2]*c[1];
  m[2] = c[2];
  n[0] = s[3]*c[1]-c[3]*c[2]*s[1];
  n[1] = -s[3]*s[1]-c[3]*c[2]*c[1];
  n[2] = c[3]*s[2];


  l[0] = c[1]*c[3] + s[1]*s[2]*s[3];
  l[1] = s[1]*c[2];
  l[2] = -c[1]*s[3] + s[1]*s[2]*c[3];
  m[0] = -s[1]*c[3] + c[1]*s[2]*s[3];
  m[1] = c[1]*c[2];
  m[2] = s[1]*s[3] + c[1]*s[2]*c[3];
  n[0] = c[2]*s[3];
  n[1] = -s[2];
  n[2] = c[2]*c[3];

  /* Rotation is first around Y-axis for rake second around X-axis for 90°-dip and finaly around Z-axis for strike */
  /* ROTATION MATRIX TO ROTATE SOURCE WITH n=(1,0,0) AND s=(0,1,0) (S0)
     INTO SOURCE GIVEN BY simpa.SourceOri (SR)
     /SR[0]\   /l[0] l[1] l[2]\   /S0[0]\
     |SR[1]| = |m[0] m[1] m[2]| * |S0[1]|
     \SR[2]/   \n[0] n[1] n[2]/   \S0[2]/
  */
  return 0;
}

int SourceRot(double * theta, double * phi, double * l, double * m, double * n)
{
  double x[3], x_old[3];
  double tmp;
  tmp = sin(*theta);
  x_old[0] = tmp*cos(*phi);
  x_old[1] = tmp*sin(*phi);
  x_old[2] = cos(*theta);
  
  x[0] = l[0]*x_old[0] + m[0]*x_old[1] + n[0]*x_old[2];
  x[1] = l[1]*x_old[0] + m[1]*x_old[1] + n[1]*x_old[2];
  x[2] = l[2]*x_old[0] + m[2]*x_old[1] + n[2]*x_old[2];

  *theta = atan2(sqrt(pow(x[0],2)+pow(x[1],2)),x[2]);
  *phi = atan2(x[1],x[0]);
  return 0;
}
