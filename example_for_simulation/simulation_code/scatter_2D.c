#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "edefinition.h"
#include "eutility.h"
#include "scatter_2D.h"


int scatterP(bodyparm body, particle * pat_p, double gpp0, float * wahr_pp, float * wahr_ps)
{

  /* Declaration */
  double r;            /* random number */
  double the;          /* scattaring angle in local coordinate system */

  /* Code */
  r = (double)rand()/RAND_MAX;
  if(1){ /*acoustic: only scattering into P-wave */
    the = (double)wahr_pp[(int)ceil(body.q *(double)rand()/RAND_MAX)-1];   /* -0.5 because round() always round away from zero */
    // 2D propagation
    if(rand()%1000 >500) pat_p->dir[1] = pat_p->dir[1] + the;
    else pat_p->dir[1] = pat_p->dir[1] - the;
    pat_p->dir[0] = PI2;
    /* stokesvector remains the same */
  }
  else{ /*scattering into SV-wave in 2D */
    pat_p->mode = 1;
    the = (double)wahr_ps[(int)ceil(body.q * (double)rand()/RAND_MAX)-1];
    if(rand()%1000 >500) pat_p->dir[1] = pat_p->dir[1] + the;
    else pat_p->dir[1] = pat_p->dir[1] - the;
    pat_p->dir[0] = PI2;
    /* stokesvector only SV-wave */
    pat_p->stokes[2] = pat_p->stokes[0];
    pat_p->stokes[0] = pat_p->stokes[1] = pat_p->stokes[3] = pat_p->stokes[4] =0.;
    pat_p->pol = PI2;
  }
  return 0;
}




int scatterS(bodyparm body, particle * pat_p, double gsp0, double gssf0, float * wahr_sp, float * wahr_sp_phi, float * wahr_ssf, float * wahr_ssf_phi, float * wahr_ssp, float * wahr_ssp_phi)
{

  /* Declaration */
  double r;            /* random number */
  double the;

  /* Code */
  r = (double)rand()/RAND_MAX;
    if (pat_p->stokes[1]-0<1e-10){/*The particle is SV-wave*/
        if(r <= gsp0/body.pg[1]){ /* scattering into P-wave */
            pat_p->mode = 0;
            the = (double)wahr_sp[(int)ceil(body.q * (double)rand()/RAND_MAX)-1];   /* -0.5 because round() always rounds away from zero */
            // 2D propagation
            if(rand()%1000 >500) pat_p->dir[1] = pat_p->dir[1] + the;
            else pat_p->dir[1] = pat_p->dir[1] - the;
            pat_p->dir[0] = PI2;
            /* new stokes parameters */
            pat_p->stokes[0] = pat_p->stokes[2];
            pat_p->stokes[1] = pat_p->stokes[2] = pat_p->stokes[3] = pat_p->stokes[4] = 0.0;
        }
        else { /* scattering into SV-wave */
            the = (double)wahr_ssf[(int)ceil(body.q * (double)rand()/RAND_MAX)-1];   /* -0.5 because round() always rounds away from zero */
            // 2D propagation
            if(rand()%1000 >500) pat_p->dir[1] = pat_p->dir[1] + the;
            else pat_p->dir[1] = pat_p->dir[1] - the;
            pat_p->dir[0] = PI2;
            /* stokesvector remains the same */
            }
    }
    else{  /* The particle is SH-wave, for 2D only scattering into SH-wave */
        the = (double)wahr_ssp[(int)ceil(body.q * (double)rand()/RAND_MAX)-1];   /* -0.5 because round() always round away from zero */
        // 2D propagation
        if(rand()%1000 >500) pat_p->dir[1] = pat_p->dir[1] + the;
        else pat_p->dir[1] = pat_p->dir[1] - the;
        pat_p->dir[0] = PI2;
        /* stokesvector remains the same */
    }
  return 0;
}

