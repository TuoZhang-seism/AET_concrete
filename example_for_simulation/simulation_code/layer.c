#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "edefinition.h"
#include "scat_coeff_2D.h"
#include "layer.h"
#include "incidence.h"
#include "eutility.h"


int MkInterface(simulationparameters simpa, float ** ptS, float ** ptP)
{
  int mode = 1;
  int i;
  double al1,al2,be1,be2;
  double ESVR, ESVT, ESHR, ESHT, EPT, EPR;
  FILE * fp;
  fp = fopen(simpa.InterfaceName,"w+");
  fprintf(fp,"Reflection and conversion coefficients for incidence on the interface.\n");
  fprintf(fp,"Rows 6 to %d contain:\n",(int) 2*simpa.Nangle+6);
  fprintf(fp,"theta of reflected S, of transmitted S, of reflected P, of transmitted P, SVR, SVT, SHR, SHT, PR, PT for incident S-wave\n");
  fprintf(fp,"Rows %d to %d contain:\n",(int) 2.*simpa.Nangle+7,(int) 4*simpa.Nangle+7);
  fprintf(fp,"theta of reflected S, of transmitted S, of reflected P, of transmitted P, SVR, SVT, SHR, SHT, PR, PT for incident P-wave\n");
  for(i=0;i<simpa.Nangle;i++){    /* S-incidence from below */
    be1 = PI2/simpa.Nangle*(i+0.5);
    incidence(&mode, &ESVR, &ESVT, &ESHR, &ESHT, &EPT, &EPR, &simpa.body[1].v[0], &simpa.body[1].v[1], &simpa.body[0].v[0], &simpa.body[0].v[1], &simpa.body[1].rho, &simpa.body[0].rho, &be1, &be2, &al1, &al2);
    ptS[0][i] = PI - be1;         /* Theta angles of reflected S */
    ptS[1][i] = be2;              /* transmitted S */
    ptS[2][i] = PI - al1;         /* reflected P */
    ptS[3][i] = al2;              /* transmitted P */
    ptS[4][i] = (float)ESVR;
    ptS[5][i] = (float)ESVT;
    ptS[6][i] = (float)ESHR;
    ptS[7][i] = (float)ESHT;
    ptS[8][i] = (float)EPR;
    ptS[9][i] = (float)EPT;
  }
  for(i=0;i<simpa.Nangle;i++){    /* S-incidence from above */
    be1 = PI2 - PI2/simpa.Nangle*(i+0.5);
    incidence(&mode, &ESVR, &ESVT, &ESHR, &ESHT, &EPT, &EPR, &simpa.body[0].v[0], &simpa.body[0].v[1], &simpa.body[1].v[0], &simpa.body[1].v[1], &simpa.body[0].rho, &simpa.body[1].rho, &be1, &be2, &al1, &al2);
    ptS[0][i + simpa.Nangle] = be1;              /* Theta angles of reflected S */
    ptS[1][i + simpa.Nangle] = PI - be2;         /* transmitted S */
    ptS[2][i + simpa.Nangle] = al1;              /* reflected P */
    ptS[3][i + simpa.Nangle] = PI - al2;         /* transmitted P */
    ptS[4][i + simpa.Nangle] = (float)ESVR;
    ptS[5][i + simpa.Nangle] = (float)ESVT;
    ptS[6][i + simpa.Nangle] = (float)ESHR;
    ptS[7][i + simpa.Nangle] = (float)ESHT;
    ptS[8][i + simpa.Nangle] = (float)EPR;
    ptS[9][i + simpa.Nangle] = (float)EPT;
  }

  mode = 0;
  for(i=0;i<simpa.Nangle;i++){    /* P-incidence from below */
    al1 = PI2/simpa.Nangle*(i+0.5);
    incidence(&mode, &ESVR, &ESVT, &ESHR, &ESHT, &EPT, &EPR, &simpa.body[1].v[0], &simpa.body[1].v[1], &simpa.body[0].v[0], &simpa.body[0].v[1], &simpa.body[1].rho, &simpa.body[0].rho, &be1, &be2, &al1, &al2);
    ptP[0][i] = PI - be1;         /* Theta angles of reflected S */
    ptP[1][i] = be2;              /* transmitted S */
    ptP[2][i] = PI - al1;         /* reflected P */
    ptP[3][i] = al2;              /* transmitted P */
    ptP[4][i] = (float)ESVR;
    ptP[5][i] = (float)ESVT;
    ptP[6][i] = (float)ESHR;
    ptP[7][i] = (float)ESHT;
    ptP[8][i] = (float)EPR;
    ptP[9][i] = (float)EPT;
  }
  for(i=0;i<simpa.Nangle;i++){    /* P-incidence from above */
    al1 = PI2 - PI2/simpa.Nangle*(i+0.5);
    incidence(&mode, &ESVR, &ESVT, &ESHR, &ESHT, &EPT, &EPR, &simpa.body[0].v[0], &simpa.body[0].v[1], &simpa.body[1].v[0], &simpa.body[1].v[1], &simpa.body[0].rho, &simpa.body[1].rho, &be1, &be2, &al1, &al2);
    ptP[0][i + simpa.Nangle] = be1;              /* Theta angles of reflected S */
    ptP[1][i + simpa.Nangle] = PI - be2;         /* transmitted S */
    ptP[2][i + simpa.Nangle] = al1;              /* reflected P */
    ptP[3][i + simpa.Nangle] = PI - al2;         /* transmitted P */
    ptP[4][i + simpa.Nangle] = (float)ESVR;
    ptP[5][i + simpa.Nangle] = (float)ESVT;
    ptP[6][i + simpa.Nangle] = (float)ESHR;
    ptP[7][i + simpa.Nangle] = (float)ESHT;
    ptP[8][i + simpa.Nangle] = (float)EPR;
    ptP[9][i + simpa.Nangle] = (float)EPT;
  }
  for(i=0;i<simpa.Nangle*2;i++){
    fprintf(fp,"%e %e %e %e %e %e %e %e %e %e \n",ptS[0][i], ptS[1][i], ptS[2][i], ptS[3][i], ptS[4][i], ptS[5][i], ptS[6][i], ptS[7][i], ptS[8][i], ptS[9][i]);
  }
  for(i=0;i<simpa.Nangle*2;i++){
    fprintf(fp,"%e %e %e %e %e %e %e %e %e %e \n",ptP[0][i], ptP[1][i], ptP[2][i], ptP[3][i], ptP[4][i], ptP[5][i], ptP[6][i], ptP[7][i], ptP[8][i], ptP[9][i]);
  }
  fclose(fp);
  return 0;
}



int ReadInterface(simulationparameters simpa, float ** ptS, float ** ptP)
{
  FILE * fp;
  char dump[301];
  int i;
  int count =0;
  for(i=0;i<10;i++){
    ptS[i] = (float *) calloc((size_t) 2*simpa.Nangle, sizeof(float));
    ptP[i] = (float *) calloc((size_t) 2*simpa.Nangle, sizeof(float));
  }
  if((fp = fopen(simpa.InterfaceName,"r"))==NULL){
    printf("%s:%d Error opening file %s\n",__FILE__,__LINE__,simpa.InterfaceName);
    return -1;
  }
  fgets(dump,300,fp);fgets(dump,300,fp);fgets(dump,300,fp);fgets(dump,300,fp);fgets(dump,300,fp);
  for(i=0;i<2*simpa.Nangle;i++){
    count += fscanf(fp,"%f %f %f %f %f %f %f %f %f %f\n",&ptS[0][i],&ptS[1][i],&ptS[2][i],&ptS[3][i],&ptS[4][i],&ptS[5][i],&ptS[6][i],&ptS[7][i],&ptS[8][i],&ptS[9][i]);
  }
  for(i=0;i<2*simpa.Nangle;i++){
    count += fscanf(fp,"%f %f %f %f %f %f %f %f %f %f\n",&ptP[0][i],&ptP[1][i],&ptP[2][i],&ptP[3][i],&ptP[4][i],&ptP[5][i],&ptP[6][i],&ptP[7][i],&ptP[8][i],&ptP[9][i]);
  }
  fclose(fp);
  if(count!=10*4*simpa.Nangle){
    printf("File %s: Number of entries did not match expected number.",simpa.InterfaceName);
    return -1;
  }
  return 0;
}



int layer(particle * pat_p, simulationparameters simpa, float ** layP_t, float ** layS_t, int bfrom, int bto)
{
  /* here the argument b is the body in which the partical is expected to appear if it crosses the interface */
  int ind;
  double r;
  double ESVR,ESVT,ESHR,ESHT,EPT,EPR;   /* transmitted and reflected energies */
  double SE;
  double frac;
  double rds, cds, tds, tcds;
  double tmp;
  double x[2];
  r = (double)rand()/(RAND_MAX);
  ind = MIN((int) (pat_p->dir[0]*2.*simpa.Nangle/PI),2*simpa.Nangle-1);  /* in case pat_p->dir[0] is == PI */
  if(pat_p->mode){  /* incident S-wave */
    ESVR = layS_t[4][ind] * pat_p->stokes[1];
    ESVT = layS_t[5][ind] * pat_p->stokes[1];
    ESHR = layS_t[6][ind] * pat_p->stokes[2];
    ESHT = layS_t[7][ind] * pat_p->stokes[2];
    EPR = layS_t[8][ind] * pat_p->stokes[1];
    EPT = layS_t[9][ind] * pat_p->stokes[1];
    SE = ESVR+ESVT+ESHR+ESHT+EPR+EPT;
    if(r<(ESVR+ESHR)/SE){ /* Reflection as S-wave */
      pat_p->dir[0] = PI - pat_p->dir[0];
      pat_p->stokes[1] = (pat_p->stokes[1]+pat_p->stokes[2])*ESVR/(ESVR+ESHR);
      pat_p->stokes[2] = (pat_p->stokes[1]+pat_p->stokes[2])*ESHR/(ESVR+ESHR);
      pat_p->stokes[3] = 2.*sqrt(pat_p->stokes[1] * pat_p->stokes[2]);
      pat_p->stokes[4] = 0.;
      pat_p->pol = atan(sqrt(pat_p->stokes[2]/pat_p->stokes[1]));
      pat_p->pos[2] = 2.*simpa.body[0].z[1] - pat_p->pos[2];
      pat_p->vec[2] = - pat_p->vec[2];
      return bfrom;
    }

    else{ /* no simple reflection */
      frac = fabs((pat_p->pos[2]-simpa.body[0].z[1])/pat_p->vec[2]);    /* ratio of path that will be changed */
      /* new path length after leaving the surface */
      rds = frac * simpa.body[bfrom].path[1];
      /* interface position */
      x[0] = pat_p->pos[0] - frac * pat_p->vec[0];
      x[1] = pat_p->pos[1] - frac * pat_p->vec[1];
      if(r<(ESVR+ESHR+EPR)/SE){ /* Reflection as P-wave */
        pat_p->mode = 0;
        cds = frac * simpa.body[bfrom].path[0];
        /* new theta */
        pat_p->dir[0] = layS_t[2][ind];
        /* correct travel time of particle */
        //pat_p->time += (cds/simpa.body[bfrom].v[0] - rds/simpa.body[bfrom].v[1]);  No correction needed because TimeStepFrac = const
        /* convert energy and correct for intrinsic atternuation */
        pat_p->stokes[0] = (pat_p->stokes[1]+pat_p->stokes[2]);
        pat_p->stokes[1] = pat_p->stokes[2] = pat_p->stokes[3] = pat_p->stokes[4] = 0.0;
        //pat_p->stokes[0] *= exp(simpa.body[bfrom].Qex[0]*(-cds)+simpa.body[bfrom].Qex[1]*rds);
        pat_p->stokes[0] *= pow(simpa.body[bfrom].Qex[0]/simpa.body[bfrom].Qex[1],frac);
        /* new direction and position */
        pat_p->pos[0] = x[0] + cds * sin(pat_p->dir[0]) * cos(pat_p->dir[1]);
        pat_p->pos[1] = x[1] + cds * sin(pat_p->dir[0]) * sin(pat_p->dir[1]);
        pat_p->pos[2] = simpa.body[0].z[1] - cds * cos(pat_p->dir[0]);
        pat_p->vec[0] = simpa.body[bfrom].path[0] * sin(pat_p->dir[0]) * cos(pat_p->dir[1]);
        pat_p->vec[1] = simpa.body[bfrom].path[0] * sin(pat_p->dir[0]) * sin(pat_p->dir[1]);
        pat_p->vec[2] = - simpa.body[bfrom].path[0] * cos(pat_p->dir[0]);
        return bfrom;
      }
      else if(r<(ESVR+ESHR+EPR+ESVT+ESHT)/SE){ /* Transmission as S-wave */
        tds = frac * simpa.body[bto].path[1];
        /* new theta */
        pat_p->dir[0] = layS_t[1][ind];
        /* correct travel time of particle */
	pat_p->time += frac*simpa.TimeStep*(1./simpa.body[bto].TimeStepFrac - 1./simpa.body[bfrom].TimeStepFrac);
        /* convert energy and correct for intrinsic atternuation */
        tmp = pow(simpa.body[bto].Qex[1]/simpa.body[bfrom].Qex[1],frac);
        pat_p->stokes[1] *= tmp;
        pat_p->stokes[2] *= tmp;
        pat_p->stokes[3] *= tmp;
        /* new direction and position */
        pat_p->pos[0] = x[0] + tds * sin(pat_p->dir[0]) * cos(pat_p->dir[1]);
        pat_p->pos[1] = x[1] + tds * sin(pat_p->dir[0]) * sin(pat_p->dir[1]);
        pat_p->pos[2] = simpa.body[0].z[1] - tds * cos(pat_p->dir[0]);
        pat_p->vec[0] = simpa.body[bto].path[1] * sin(pat_p->dir[0]) * cos(pat_p->dir[1]);
        pat_p->vec[1] = simpa.body[bto].path[1] * sin(pat_p->dir[0]) * sin(pat_p->dir[1]);
        pat_p->vec[2] = - simpa.body[bto].path[1] * cos(pat_p->dir[0]);
        return bto;
      }
      else{ /* Transmission as P-wave */
        pat_p->mode = 0;
        tcds = frac * simpa.body[bto].path[0];
        /* new theta */
        pat_p->dir[0] = layS_t[3][ind];
        /* correct travel time of particle */
        //pat_p->time += (tcds/simpa.body[bto].v[0] - rds/simpa.body[bfrom].v[1]);
	pat_p->time += frac*simpa.TimeStep*(1./simpa.body[bto].TimeStepFrac - 1./simpa.body[bfrom].TimeStepFrac);
        /* convert energy and correct for intrinsic atternuation */
        pat_p->stokes[0] = (pat_p->stokes[1]+pat_p->stokes[2]);
        pat_p->stokes[1] = pat_p->stokes[2] = pat_p->stokes[3] = pat_p->stokes[4] = 0.0;
	pat_p->stokes[0] *= pow(simpa.body[bto].Qex[0]/simpa.body[bfrom].Qex[1],frac);
        /* new direction and position */
        pat_p->pos[0] = x[0] + tcds * sin(pat_p->dir[0]) * cos(pat_p->dir[1]);
        pat_p->pos[1] = x[1] + tcds * sin(pat_p->dir[0]) * sin(pat_p->dir[1]);
        pat_p->pos[2] = simpa.body[0].z[1] - tcds * cos(pat_p->dir[0]);
        pat_p->vec[0] = simpa.body[bto].path[0] * sin(pat_p->dir[0]) * cos(pat_p->dir[1]);
        pat_p->vec[1] = simpa.body[bto].path[0] * sin(pat_p->dir[0]) * sin(pat_p->dir[1]);
        pat_p->vec[2] = - simpa.body[bto].path[0] * cos(pat_p->dir[0]);
        return bto;
      }
    }
  }

  else{ /* incident P-wave */
/* 
    ESVR = layP_t[ind][4] * pat_p.stokes[0];
    ESVT = layP_t[ind][5] * pat_p.stokes[0];
    EPR = layP_t[ind][8] * pat_p.stokes[0];
    EPT = layP_t[ind][9] * pat_p.stokes[0];

    SE = ESVR+ESVT+EPR+EPT; assume this is one
*/

    if(r<layP_t[8][ind]){ /* Reflection as P-wave */
      pat_p->dir[0] = PI - pat_p->dir[0];
      pat_p->pos[2] = 2.*simpa.body[0].z[1] - pat_p->pos[2];
      pat_p->vec[2] = - pat_p->vec[2];
      return bfrom;
    }
    else{ /* no simple reflection */
      frac = fabs((pat_p->pos[2]-simpa.body[0].z[1])/pat_p->vec[2]);
      /* new path length after leaving the surface */
      rds = frac * simpa.body[bfrom].path[0];
      /* interface position */
      x[0] = pat_p->pos[0] - frac * pat_p->vec[0];
      x[1] = pat_p->pos[1] - frac * pat_p->vec[1];

      if(r<(layP_t[4][ind] + layP_t[8][ind])){ /* Reflection as S-wave */
        pat_p->mode = 1;
        cds = frac * simpa.body[bfrom].path[1];
        /* new theta */
        pat_p->dir[0] = layP_t[0][ind];
        /* correct travel time of particle */
        // pat_p->time += (cds/simpa.body[bfrom].v[1] - rds/simpa.body[bfrom].v[0]);  No correction needed
        /* convert energy and correct for intrinsic atternuation */
        pat_p->stokes[1] = pat_p->stokes[0]*pow(simpa.body[bfrom].Qex[1]/simpa.body[bfrom].Qex[0],frac);
        pat_p->stokes[0] = pat_p->stokes[2] = pat_p->stokes[3] = pat_p->stokes[4] = 0.0;
        /* new direction and position */
        pat_p->pos[0] = x[0] + cds * sin(pat_p->dir[0]) * cos(pat_p->dir[1]);
        pat_p->pos[1] = x[1] + cds * sin(pat_p->dir[0]) * sin(pat_p->dir[1]);
        pat_p->pos[2] = simpa.body[0].z[1] - cds * cos(pat_p->dir[0]);
        pat_p->vec[0] = simpa.body[bfrom].path[1] * sin(pat_p->dir[0]) * cos(pat_p->dir[1]);
        pat_p->vec[1] = simpa.body[bfrom].path[1] * sin(pat_p->dir[0]) * sin(pat_p->dir[1]);
        pat_p->vec[2] = - simpa.body[bfrom].path[1] * cos(pat_p->dir[0]);
        return bfrom;
      }
      else if(r<(layP_t[4][ind] + layP_t[8][ind] + layP_t[5][ind])){ /* Transmission as S-wave */
        pat_p->mode = 1;
        tcds = frac * simpa.body[bto].path[1];
        /* new theta */
        pat_p->dir[0] = layP_t[1][ind];
        /* correct travel time of particle */
	pat_p->time += frac*simpa.TimeStep*(1./simpa.body[bto].TimeStepFrac - 1./simpa.body[bfrom].TimeStepFrac);
        /* convert energy and correct for intrinsic atternuation */
        //pat_p->stokes[1] = pat_p->stokes[0]*exp(rds*simpa.body[bfrom].Qex[0]-tcds*simpa.body[bto].Qex[1]);
        pat_p->stokes[1] = pat_p->stokes[0] * pow(simpa.body[bto].Qex[1]/simpa.body[bfrom].Qex[0],frac);
        pat_p->stokes[0] = pat_p->stokes[2] = pat_p->stokes[3] = pat_p->stokes[4] = 0.0;
        /* new direction and position */
        pat_p->pos[0] = x[0] + tcds * sin(pat_p->dir[0]) * cos(pat_p->dir[1]);
        pat_p->pos[1] = x[1] + tcds * sin(pat_p->dir[0]) * sin(pat_p->dir[1]);
        pat_p->pos[2] = simpa.body[0].z[1] - tcds * cos(pat_p->dir[0]);
        pat_p->vec[0] = simpa.body[bto].path[1] * sin(pat_p->dir[0]) * cos(pat_p->dir[1]);
        pat_p->vec[1] = simpa.body[bto].path[1] * sin(pat_p->dir[0]) * sin(pat_p->dir[1]);
        pat_p->vec[2] = - simpa.body[bto].path[1] * cos(pat_p->dir[0]);
        return bto;
      }
      else{ /* Transmission as P-wave */
        tds = frac * simpa.body[bto].path[1];
        /* new theta */
        pat_p->dir[0] = layP_t[3][ind];
        /* correct travel time of particle */
	pat_p->time += frac*simpa.TimeStep*(1./simpa.body[bto].TimeStepFrac - 1./simpa.body[bfrom].TimeStepFrac);
        /* convert energy and correct for intrinsic atternuation */
        //pat_p->stokes[0] *= exp(rds*simpa.body[bfrom].Qex[0]-tds*simpa.body[bto].Qex[0]);
        pat_p->stokes[0] *= pow(simpa.body[bto].Qex[0]/simpa.body[bfrom].Qex[0],frac);
        /* new direction and position */
        pat_p->pos[0] = x[0] + tds * sin(pat_p->dir[0]) * cos(pat_p->dir[1]);
        pat_p->pos[1] = x[1] + tds * sin(pat_p->dir[0]) * sin(pat_p->dir[1]);
        pat_p->pos[2] = simpa.body[0].z[1] - tds * cos(pat_p->dir[0]);
        pat_p->vec[0] = simpa.body[bto].path[0] * sin(pat_p->dir[0]) * cos(pat_p->dir[1]);
        pat_p->vec[1] = simpa.body[bto].path[0] * sin(pat_p->dir[0]) * sin(pat_p->dir[1]);
        pat_p->vec[2] = - simpa.body[bto].path[0] * cos(pat_p->dir[0]);
        return bto;
      }
    }
  }
}
