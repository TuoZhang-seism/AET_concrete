#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include "edefinition.h"


int incidence(int * mode, double* ESVR, double* ESVT, double* ESHR, double* ESHT, double* EPT, double* EPR, double* v1p, double* v1s,  double* v2p, double* v2s, double* rho1, double* rho2, double* be1, double* be2, double* al1, double* al2)
{

  complex double z, cbe1, cbe2, cal1, cal2;
  complex double a,b,c,d,E,F,G,H,K,Delta;
  complex double svsr,svst,shsr,shst,spr,spt,ppr,ppt,psr,pst;
  complex double cbe1v,cbe2v,cal1v,cal2v;
  complex double esvsr,esvst,eshsr,eshst,espr,espt,epsr,epst,eppt,eppr;
  double p;

  /* calculate ray parameter */
  if(*mode == 1){  // Incident S-wave
    p = sin(*be1)/ *v1s;
  }
  else{  // Incident P-wave
    p = sin(*al1)/ *v1p;
  }


  /* calculate global parameters  */
  /* cos of angle divided by velocity (e.g. cos(be1)/v1s) */
  cbe1v = csqrt(cpow(*v1s,-2)-cpow(p,2));
  cbe2v = csqrt(cpow(*v2s,-2)-cpow(p,2));
  cal1v = csqrt(cpow(*v1p,-2)-cpow(p,2));
  cal2v = csqrt(cpow(*v2p,-2)-cpow(p,2));

  /* Parameters of Aki 1980 p.149 eq.5.38 */
  a = *rho2 * (1.-2.*cpow(*v2s *p,2)) - *rho1 *(1.-2.*cpow(*v1s*p,2));
  b = *rho2 * (1.-2.*cpow(*v2s*p,2)) + 2 * *rho1 * cpow(*v1s*p,2);
  c = *rho1 * (1.-2.*cpow(*v1s*p,2)) + 2 * *rho2 * cpow(*v2s*p,2);
  d = 2*(*rho2 * *v2s * *v2s - *rho1 * *v1s * *v1s);
  E = b * cal1v + c * cal2v;
  F = b * cbe1v + c * cbe2v;
  G = a - d * cal1v * cbe2v;
  H = a - d * cal2v * cbe1v;
  K = E * F + G * H * p * p;
  Delta = *rho1 * cpow(*v1s,2) * cbe1v + *rho2 * cpow(*v2s,2) * cbe2v;

  /*###################################################
  ## Incident S-wave ##################################
  ###################################################*/
  if(*mode == 1){  // Incident S-wave
    /* P-wave transmission */
    if(cimag(cal2v) == 0){
      /* SV\P\ */
      spt =  (-2. * *rho1 * *v1s * cbe1v * G * p / (*v2p * K));
      espt = spt*csqrt((*rho2 * cpow(*v2p,2) * cal2v) / (*rho1 * cpow(*v1s,2) * cbe1v));
      *EPT = (double)(espt * conj(espt));
      *al2 = (double)casin(p* *v2p);
    }
    else{
      *EPT = 0.;
      *al2 = 111;
    }
    /* S-wave transmission */
    if (cimag(cbe2v) == 0){
      /* SV\SV\  v1[1] gekürzt */
      svst = (2 * *rho1 * *v1s * cbe1v * E / (*v2s * K));
      /* SH\SH\ */
      shst = (2. * *rho1 * pow(*v1s,2) * cbe1v / Delta);
      esvst = svst*csqrt((* rho2 * pow(*v2s,2) * cbe2v) / (*rho1 * pow(*v1s,2) * cbe1v));
      eshst = shst * csqrt((*rho2 * pow(*v2s,2) * cbe2v) / (*rho1 * pow(*v1s,2) * cbe1v));
      *ESVT = esvst * conj(esvst);
      *ESHT = eshst * conj(eshst);
      *be2 = casin(p * *v2s);
    }
    else{
      *ESVT = 0.;
      *ESHT = 0.;
      *be2 = 111;
    }
    /* P-wave reflection */
    if (cimag(cal1v) == 0){
      /* SV\P/ */
      spr = (-2. * cbe1v * (a * b + c * d * cal2v * cbe2v) * p * *v1s / (*v1p * K));
      espr = spr * csqrt((cpow(*v1p,2) * cal1v) / (pow(*v1s,2) * cbe1v));
      *EPR = espr*conj(espr);
      *al1 = casin(p * *v1p);
    }
    else{
      *EPR = 0.;
      *al1 = 111;
    }
    /* S-wave reflection */
    /* SH\SH/ */
    shsr = (*rho1 * cpow(*v1s,2) * cbe1v - *rho2 * cpow(*v2s,2) * cbe2v) / Delta;
    /* SV\SV/ */
    svsr = (-((b * cbe1v - c * cbe2v) * E - (a + d* cal2v * cbe1v) * G * p * p) / K);
    *ESVR = svsr * conj(svsr);
    *ESHR = shsr * conj(shsr);
  }

  /*####################################################
  ## Incident P-wave ###################################
  ####################################################*/
  else if(*mode == 0){
    /* P-wave transmission */
    if (cimag(cal2v) == 0){
      /* P\P\ */
      ppt = 2. * *rho1 * cal1v * F * *v1p / (*v2p * K);
      eppt = ppt * csqrt((*rho2 * cpow(*v2p,2) * cal2v) / (*rho1 * cpow(*v1p,2) * cal1v));
      *EPT = eppt*conj(eppt);
      *al2 = casin(p * *v2p);
    }
    else{
      *EPT = 0.;
      *al2 = 111;
    }
    /* S-wave transmission */
    if (cimag(cbe2v) == 0){
      /* P\S\ */
      pst = 2. * *rho1 * cal1v * H * p * *v1p / (*v2s * K);
      epst = pst * csqrt((*rho2 * cpow(*v2s,2) * cbe2v) / (*rho1 * cpow(*v1p,2) * cal1v));
      *ESVT = epst * conj(epst);
      *ESHT = 0.;
      *be2 = casin(p * *v2s);
    }
    else{
      *ESVT = 0.;
      *ESHT = 0.;
      *be2 = 111;
    }
    /* P-wave reflection */
    /* P\P/ */
    ppr = ((b * cal1v - c * cal2v) * F - (a + d * cal1v * cbe2v) * H * p * p) / K;
    *EPR = ppr * conj(ppr);
    /* S-wave reflection */
    /* P\S/ */
    psr = -2. * cal1v * (a * b + c * d * cal2v * cbe2v) * p * *v1p / (*v1s * K);
    epsr = psr * csqrt((cpow(*v1s,2) * cbe1v) / (cpow(*v1p,2) * cal1v));
    *ESVR = epsr * conj(epsr);
    *ESHR = 0.;
    *be1 = casin(p * *v1s);
  }
  
  /* Check the velocity of the transmitted S-wave. S-velocity below 0.2 is assumed to designate a liquid layer.
  Put ESHT and ESVT to zero in this case and add their value to EPT */
  if(*v2s<0.2){
    if(*mode){
      *ESVR += *ESVT;
      *ESVT = 0.;
      *ESHR += *ESHT;
      *ESHT = 0.;
    }
    else{
      *EPR += *ESVT;
      *ESVT = 0.;
    }
  }
  /* Check the velocity of the reflected S-wave. S-velocity 0.2 is assumed to designate a liquid layer.
  Put ESVR to zero and add its values to EPT */
  if(*v1s<0.2){
    *EPT += *ESVR;
    *ESVR = 0.;
  }
  
  return 0;
}
