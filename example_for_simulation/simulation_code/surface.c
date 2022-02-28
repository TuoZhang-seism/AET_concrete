#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "edefinition.h"
#include "scat_coeff_2D.h"
#include "surface.h"
#include "eutility.h"


int MkSurface(simulationparameters simpa, float ** surfS, float ** surfP)
{
  int i;
  double p, c[2];
  float ssv, ssh, sp, ps, pp;
  double alpha,beta;
  FILE * fp;
  fp = fopen(simpa.SurfaceName,"w+");
  fprintf(fp,"Reflection and conversion coefficients for incidence on the interface.\n");
  fprintf(fp,"Rows 6 to %d contain:\n",(int) 2*simpa.Nangle+6);
  fprintf(fp,"theta of reflected S, and of reflected P, SVR, SHR, PR for incident S-wave\n");
  fprintf(fp,"Rows %d to %d contain:\n",(int) 2.*simpa.Nangle+7,(int) 4*simpa.Nangle+7);
  fprintf(fp,"theta of reflected S, and of reflected P, PPR, PSR for incident P-wave\n");
  for(i=0;i<simpa.Nangle;i++){    /* S-incidence */
    alpha = PI2/simpa.Nangle*(i+0.5);
    p = sin(alpha)/simpa.body[0].v[1];         /* ray parameter */
    if(p*simpa.body[0].v[0]<=1){  /* angle of converted P-wave < PI */
      beta = asin(p*simpa.body[0].v[0]);
      c[0] = 1./pow(simpa.body[0].v[1],2)-2.*pow(p,2);
      c[1] = 4.* pow(p,2)*cos(beta)/simpa.body[0].v[0]*cos(alpha)/simpa.body[0].v[1];
      /* The following equations result from Aki[1980] p.140 and eq.5.40 on p.151
	 and correspond to eq.32-34 of margerin[2004] */
      sp = /* pow(simpa.body[0].gam,-2.) * */ cos(beta)/cos(alpha)*simpa.body[0].gam*pow((4./simpa.body[0].v[0]*p*cos(alpha)*c[0])/(pow(c[0],2) + c[1]),2);
      ssv = pow((pow(c[0],2)-c[1])/(pow(c[0],2)+c[1]),2);
    }
    else{
      beta = -1.;
      sp = 0.;
      ssv = 1.;
    }
    ssh = 1;
    surfS[0][i] = alpha;
    surfS[1][i] = PI - beta;
    surfS[2][i] = ssv;
    surfS[3][i] = ssh;
    surfS[4][i] = sp;
    fprintf(fp,"%e %e %e %e %e\n",surfS[0][i], surfS[1][i], surfS[2][i], surfS[3][i], surfS[4][i]);
  }
  for(i=0;i<simpa.Nangle;i++){    /* P-incidence */
    beta = PI2/simpa.Nangle*(i+0.5);
    p = sin(beta)/simpa.body[0].v[0];         /* ray parameter */
    alpha =  asin(p*simpa.body[0].v[1]);
    c[0] = 1./pow(simpa.body[0].v[1],2)-2.*pow(p,2);
    c[1] = 4.* pow(p,2)*cos(beta)/simpa.body[0].v[0]*cos(alpha)/simpa.body[0].v[1];
    pp = pow((pow(c[0],2)-c[1])/(pow(c[0],2)+c[1]),2);
    ps = /* pow(simpa.body[0].gam,2.)* */cos(alpha)/cos(beta)/simpa.body[0].gam*pow((4./simpa.body[0].v[1]*p*cos(beta)*c[0])/(pow(c[0],2) + c[1]),2);
    surfP[0][i] = beta;
    surfP[1][i] = PI - alpha;
    surfP[2][i]= pp;
    surfP[3][i] = ps;
    fprintf(fp,"%e %e %e %e\n",surfP[0][i], surfP[1][i], surfP[2][i], surfP[3][i]);
  }
  fclose(fp);
  return 0;
}


int ReadSurface(simulationparameters simpa, float ** surfS, float ** surfP)
{
  FILE * fp;
  char dump[301];
  int i;
  int count =0;
  for(i=0;i<5;i++){
    surfS[i] = (float *) calloc((size_t) simpa.Nangle, sizeof(float));
  }
  for(i=0;i<4;i++){
    surfP[i] = (float *) calloc((size_t) simpa.Nangle, sizeof(float));
  }
  if((fp = fopen(simpa.SurfaceName,"r"))==NULL){
    printf("%s:%d Error opening file %s\n",__FILE__,__LINE__,simpa.SurfaceName);
    return -1;
  }
  fgets(dump,300,fp);fgets(dump,300,fp);fgets(dump,300,fp);fgets(dump,300,fp);fgets(dump,300,fp);
  for(i=0;i<simpa.Nangle;i++){
    count += fscanf(fp,"%f %f %f %f %f\n",&surfS[0][i],&surfS[1][i],&surfS[2][i],&surfS[3][i],&surfS[4][i]);
  }
  for(i=0;i<simpa.Nangle;i++){
    count += fscanf(fp,"%f %f %f %f\n",&surfP[0][i],&surfP[1][i],&surfP[2][i],&surfP[3][i]);
  }
  fclose(fp);
  if(count!=9*simpa.Nangle){
    printf("File %s: Number of entries did not match expected number.",simpa.SurfaceName);
    return -1;
  }
  return 0;
}

int surface(particle * pat_p, int Nangle, bodyparm body, float ** surfP_t, float ** surfS_t, double * path)
{
  int ind;              /* index in table */
  double Es, Ep;
  double frac;
  double rds, cds;    /* path  length of reflected and converted wave after leaving the surface */
  double x[2];
  ind = (int) (pat_p->dir[0]*Nangle/PI2);
  if(pat_p->mode){  /* incident S wave */
    Es = surfS_t[3][ind] * pat_p->stokes[2] + surfS_t[2][ind]*pat_p->stokes[1];
    Ep = surfS_t[4][ind] * pat_p->stokes[1];
    if(Ep/(Ep + Es)>(double)rand()/(RAND_MAX)){  /* Conversion into P-wave */
      frac = pat_p->pos[2]/pat_p->vec[2];
      /* new path length after leaving the surface */
      rds = frac * path[1];
      cds = frac * path[0];
      /* surface position */
      x[0] = pat_p->pos[0] - frac * pat_p->vec[0];
      x[1] = pat_p->pos[1] - frac * pat_p->vec[1];
      /* new direction and position */
      pat_p->dir[0] = surfS_t[1][ind];
      pat_p->vec[0] = path[0] * sin(pat_p->dir[0]) * cos(pat_p->dir[1]);
      pat_p->vec[1] = path[0] * sin(pat_p->dir[0]) * sin(pat_p->dir[1]);
      pat_p->vec[2] = - path[0] * cos(pat_p->dir[0]);
      pat_p->pos[0] = x[0] + cds * sin(pat_p->dir[0]) * cos(pat_p->dir[1]);
      pat_p->pos[1] = x[1] + cds * sin(pat_p->dir[0]) * sin(pat_p->dir[1]);
      pat_p->pos[2] = - cds * cos(pat_p->dir[0]);
      /* convert energy and correct for intrinsic atternuation */
      pat_p->mode = 0;
      pat_p->stokes[0] = (pat_p->stokes[1] + pat_p->stokes[2]);
      pat_p->stokes[1] = pat_p->stokes[2] = pat_p->stokes[3] = pat_p->stokes[4] = 0.0;
      //pat_p->stokes[0] *= exp(body.Qex[1]*rds - body.Qex[0]*cds);
      pat_p->stokes[0] *= pow(body.Qex[1],rds/body.path[1]-1.) * pow(body.Qex[0],cds/body.v[0]);

      /* correct travel time of particle */
      // pat_p->time += (cds/body.v[0] - rds/body.v[1]); this is bullshit too
    }
    else{  /* Reflection as S-wave */
      pat_p->pos[2] = -pat_p->pos[2];
      pat_p->dir[0] = PI - pat_p->dir[0];
      pat_p->vec[2] = -pat_p->vec[2];
    }
  }
  else{    /* incident P-wave */
    if(surfP_t[3][ind]/(surfP_t[3][ind] + surfP_t[2][ind])>(double)rand()/(RAND_MAX)){  /* Conversion into S-wave */
      frac = pat_p->pos[2]/pat_p->vec[2];
      /* new path length after leaving the surface */
      rds = frac * path[0];
      cds = frac * path[1];
      /* surface position */
      x[0] = pat_p->pos[0] - frac * pat_p->vec[0];
      x[1] = pat_p->pos[1] - frac * pat_p->vec[1];
      /* new direction and position */
      pat_p->dir[0] = surfP_t[1][ind];
      pat_p->vec[0] = path[1] * sin(pat_p->dir[0]) * cos(pat_p->dir[1]);
      pat_p->vec[1] = path[1] * sin(pat_p->dir[0]) * sin(pat_p->dir[1]);
      pat_p->vec[2] = - path[1] * cos(pat_p->dir[0]);
      pat_p->pos[0] = x[0] + cds * sin(pat_p->dir[0]) * cos(pat_p->dir[1]);
      pat_p->pos[1] = x[1] + cds * sin(pat_p->dir[0]) * sin(pat_p->dir[1]);
      pat_p->pos[2] = - cds * cos(pat_p->dir[0]);
      /* convert energy and correct for intrinsic atternuation */
      pat_p->mode = 1;
      pat_p->stokes[1] = pat_p->stokes[0];
      pat_p->stokes[0] = pat_p->stokes[2] = pat_p->stokes[3] = pat_p->stokes[4] = 0.0;
      //pat_p->stokes[1] *= exp(body.Qex[0]*rds - body.Qex[1]*cds);
      pat_p->stokes[1] *= pow(body.Qex[0],rds/body.path[0]-1.) * pow(body.Qex[1],cds/body.v[1]);
      /* correct travel time of particle */
      // pat_p->time += (cds/body.v[1] - rds/body.v[0]); this is bullshit
    }
    else{   /* reflection as P-wave */
      pat_p->pos[2] = -pat_p->pos[2];
      pat_p->dir[0] = PI - pat_p->dir[0];
      pat_p->vec[2] = -pat_p->vec[2];
    }
  }
  return 0;
}
