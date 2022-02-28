
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
#include "mk_scat_pattern.h"
#include "eutility.h"





int main(int ARGC, char *ARGV[])
{
  /* Declaration */
  bodyparm body;
  float * fltpt[6];
  double g0[10];
  int status;



  /* Code */
  if(ARGC != 9){
    printf("Please give 8 arguments:\n\tS-wave velocity\n\tRatio of P- to S-velocity\n\tDensity\n\tS-Wave-number\n\tCorrelation length\n\tFluctuation strength\n\tBirch ny\n\tNumber of lines in the probability table\n");
    return -1;
  }

  body.v[1] = atof(ARGV[1]);
  body.gam = atof(ARGV[2]);
  body.v[0] = body.gam * body.v[1];
  body.rho = atof(ARGV[3]);
  body.l = atof(ARGV[4]);
  body.a = atof(ARGV[5]);
  body.e = atof(ARGV[6]);
  body.ny = atof(ARGV[7]);
  body.q = atof(ARGV[8]);
  sprintf(body.type,"%s", "expo");
  sprintf(body.ProbTabFileName,"tables/%s_%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_%1.2e_scat",body.type,body.gam,body.l,body.a,body.e,body.ny, body.q);
  status = ReadScatPattern(body, fltpt, g0);
  if(status != 0){
    printf("Recalculating propability tables for scattering pattern.\n");
    status = MkScatPattern(body, fltpt, g0);
  }

  return 0;
}
