#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "edefinition.h"
#include "scat_coeff_2D.h"
#include "eutility.h"


int MkScatPattern(bodyparm body,  float ** fltpt, double * g0)
{
  char tmpc[121];
  FILE  *tot_file, *pat_file;
  float *pp, *ps, *ssl, *ssr, *sins, *coss;
  double  *pp_b, *ps_b, *ssl_b, *ssr_b, *phi_b;
  double  *pp_b_save, *ps_b_save, *ssl_b_save, *ssr_b_save, *phi_b_save;
  double pp_max=0, ps_max=0, ssl_max=0, ssr_max=0, phi_max=0;
  double pp_max_save=0, ps_max_save=0, ssl_max_save=0, ssr_max_save=0, phi_max_save=0;
  double pp_sum=0, ps_sum=0, ssl_sum=0, ssr_sum=0, phi_sum=0;
  double pp_cos=0, ps_cos=0, ssl_cos=0, ssr_cos=0, * cosr;
  float *r, *r2;
  float *pp_r2, *ps_r2, *ssl_r2, *ssr_r2, *phi_r2;
  int pp_num=0,  ps_num=0, ssl_num=0, ssr_num=0, phi_num=0;
  double gpp_0, gps_0, gsp_0, gssl_0, gssr_0;
  double sinr, tmpa[3], m[3], tmp;
  int i, tmpi, count=0;

  pp = (float *) calloc((size_t) body.q, (size_t) sizeof(float));
  ps = (float *) calloc((size_t) body.q, (size_t) sizeof(float));
  ssl = (float *) calloc((size_t) body.q, (size_t) sizeof(float));
  ssr = (float *) calloc((size_t) body.q, (size_t) sizeof(float));
  coss = (float *) calloc((size_t) body.q, (size_t) sizeof(float));
  sins = (float *) calloc((size_t) body.q, (size_t) sizeof(float));

  pp_b = (double *) calloc((size_t) body.q, (size_t) sizeof(double));
  ps_b = (double *) calloc((size_t) body.q, (size_t) sizeof(double));
  ssl_b = (double *) calloc((size_t) body.q, (size_t) sizeof(double));
  ssr_b = (double *) calloc((size_t) body.q, (size_t) sizeof(double));
  phi_b = (double *) calloc((size_t) body.q, (size_t) sizeof(double));

  pp_b_save = (double *) calloc((size_t) body.q, (size_t) sizeof(double));
  ps_b_save = (double *) calloc((size_t) body.q, (size_t) sizeof(double));
  ssl_b_save = (double *) calloc((size_t) body.q, (size_t) sizeof(double));
  ssr_b_save = (double *) calloc((size_t) body.q, (size_t) sizeof(double));
  phi_b_save = (double *) calloc((size_t) body.q, (size_t) sizeof(double));

  r = (float *) calloc((size_t) body.q, (size_t) sizeof(float));
  r2 = (float *) calloc((size_t) body.q, (size_t) sizeof(float));
  cosr = (double *) calloc((size_t) body.q, (size_t) sizeof(double));

  pp_r2 = (float *) calloc((size_t) body.q, (size_t) sizeof(float));
  ps_r2 = (float *) calloc((size_t) body.q, (size_t) sizeof(float));
  ssl_r2 = (float *) calloc((size_t) body.q, (size_t) sizeof(float));
  ssr_r2 = (float *) calloc((size_t) body.q, (size_t) sizeof(float));
  phi_r2 = (float *) calloc((size_t) body.q, (size_t) sizeof(float));

  while((pp_num<body.q) | (ps_num<body.q) | (ssl_num<body.q) | (ssr_num<body.q) | (phi_num<body.q)){
    count ++;
    for(i=0;i<body.q;i++){           /* scan pattern */
      r[i] = PI*rand()/RAND_MAX;
      sinr = sin(r[i]);           /* differential of spherical coordinates */
      cosr[i] = cos(r[i]);        /* cos(theta) for calculation of transport mean free path */
      m[0] = 2.*body.l/body.gam*sin(r[i]/2.);
      m[1] = body.l/body.gam*sqrt(1+pow(body.gam,2)-2.*body.gam*cos(r[i]));
      m[2] = 2.*body.l*sin(r[i]/2.);
      if (!strcmp(body.type,"gauss")){
        tmpa[0] = PSDF_GAUSS_ANG(m[0],body.a);
        tmpa[1] = PSDF_GAUSS_ANG(m[1],body.a);
        tmpa[2] = PSDF_GAUSS_ANG(m[2],body.a);
      }
      else if (!strcmp(body.type,"expo")){
        tmpa[0] = PSDF_EXP_ANG(m[0],body.a);
      	tmpa[1] = PSDF_EXP_ANG(m[1],body.a);
      	tmpa[2] = PSDF_EXP_ANG(m[2],body.a);
      }
      else if (!strcmp(body.type,"karman")){
        tmpa[0] = PSDF_KARMAN_ANG(m[0],body.a,body.kap);
        tmpa[1] = PSDF_KARMAN_ANG(m[1],body.a,body.kap);
        tmpa[2] = PSDF_KARMAN_ANG(m[2],body.a,body.kap);
      }
      else if (!strcmp(body.type,"isotrop")){
        tmpa[0] = 1.;
        tmpa[1] = 1.;
        tmpa[2] = 1.;
      }
      else{
        printf("Wrong type of ACF given in parameter file.\n");
        exit(-1);
      }
      if (!strcmp(body.type,"isotrop")){          /* Isotropy has no basic scattering pattern */
        pp_b[i] = ps_b[i] = ssl_b[i] = ssr_b[i] = tmpa[0]*sinr;
        phi_b[i] = 1.;
      }
      else{
        //pp_b[i] = pow(XPPP_THETA(body.gam, body.ny, r[i]),2)*tmpa[0]*sinr;
        //ps_b[i] = pow(XPSL_THETA(body.gam, body.ny, r[i]),2)*tmpa[1]*sinr;
        //ssl_b[i] = pow(XSSL_THETA(body.gam, body.ny, r[i]),2)*tmpa[2]*sinr;
        //ssr_b[i] = pow(XSSR_THETA(body.gam, body.ny, r[i]),2)*tmpa[2]*sinr;
          
        /*for the 2D*/
          pp_b[i] = tmpa[0];
          ps_b[i] = pow(XPSL_THETA(body.gam, body.ny, r[i]),2)*tmpa[1];
          ssl_b[i] = pow(XSSL_THETA(body.gam, body.ny, r[i]),2)*tmpa[2];
          ssr_b[i] = pow(XSSR_THETA(body.gam, body.ny, r[i]),2)*tmpa[2];
        phi_b[i] = pow(XSSR_PHI(r[i]*2.),2);
      }
      r2[i] = ((float)rand())/RAND_MAX;
    }
    for(i=0;i<body.q;i++){          /* find maximum and sum */
      if(pp_b[i]>pp_max)pp_max = pp_b[i];
      if(ps_b[i]>ps_max)ps_max = ps_b[i];
      if(ssl_b[i]>ssl_max)ssl_max = ssl_b[i];
      if(ssr_b[i]>ssr_max)ssr_max = ssr_b[i];
      if(phi_b[i]>phi_max)phi_max = phi_b[i];
      pp_sum += pp_b[i];
      ps_sum += ps_b[i];
      ssl_sum += ssl_b[i];
      ssr_sum += ssr_b[i];
      phi_sum += phi_b[i];
      pp_cos += pp_b[i]*cosr[i];
      ps_cos += ps_b[i]*cosr[i];
      ssl_cos += ssl_b[i]*cosr[i];
      ssr_cos += ssr_b[i]*cosr[i];
    }
    if(pp_max != pp_max_save){    /* in case maximum has changed, i.e. probability distribution is better sampled */
      tmpi = pp_num;
      pp_num = 0;
      pp_max_save = pp_max;
      for(i=0;i<tmpi;i++){
        if(pp_b_save[i]/pp_max/pp_r2[i]>=1){
          pp[pp_num] = pp[i];
          pp_b_save[pp_num] = pp_b_save[i];
          pp_r2[pp_num] = pp_r2[i];
          pp_num ++;
        }
      }
    }
    if(ps_max != ps_max_save){    /* in case maximum has changed, i.e. probability distribution is better sampled */
      tmpi = ps_num;
      ps_num = 0;
      ps_max_save = ps_max;
      for(i=0;i<tmpi;i++){
        if(ps_b_save[i]/ps_max/ps_r2[i]>=1){
          ps[ps_num] = ps[i];
          ps_b_save[ps_num] = ps_b_save[i];
          ps_r2[ps_num] = ps_r2[i];
          ps_num ++;
        }
      }
    }
    if(ssl_max != ssl_max_save){    /* in case maximum has changed, i.e. probability distribution is better sampled */
      tmpi = ssl_num;
      ssl_num = 0;
      ssl_max_save = ssl_max;
      for(i=0;i<tmpi;i++){
        if(ssl_b_save[i]/ssl_max/ssl_r2[i]>=1){
          ssl[ssl_num] = ssl[i];
          ssl_b_save[ssl_num] = ssl_b_save[i];
          ssl_r2[ssl_num] = ssl_r2[i];
          ssl_num ++;
        }
      }
    }
    if(ssr_max != ssr_max_save){    /* in case maximum has changed, i.e. probability distribution is better sampled */
      tmpi = ssr_num;
      ssr_num = 0;
      ssr_max_save = ssr_max;
      for(i=0;i<tmpi;i++){
        if(ssr_b_save[i]/ssr_max/ssr_r2[i]>=1){
          ssr[ssr_num] = ssr[i];
          ssr_b_save[ssr_num] = ssr_b_save[i];
          ssr_r2[ssr_num] = ssr_r2[i];
          ssr_num ++;
        }
      }
    }
    if(phi_max != phi_max_save){    /* in case maximum has changed, i.e. probability distribution is better sampled */
      tmpi = phi_num;
      phi_num = 0;
      phi_max_save = phi_max;
      for(i=0;i<tmpi;i++){
        if(phi_b_save[i]/phi_max/phi_r2[i]>=1){
          sins[phi_num] = sins[i];
          phi_b_save[phi_num] = phi_b_save[i];
          phi_r2[phi_num] = phi_r2[i];
          phi_num ++;
        }
        else{
          tmpi=tmpi;
          printf("%d\n",phi_num);
        }
      }
    }
    for(i=0;i<body.q;i++){          /* find new directions */
      if((pp_b[i]/pp_max/r2[i]>=1) & (pp_num<body.q)){
        pp[pp_num] = r[i];
        pp_b_save[pp_num] = pp_b[i];
        pp_r2[pp_num] = r2[i];
        pp_num ++;
      }
      if((ps_b[i]/ps_max/r2[i]>=1) & (ps_num<body.q)){
        ps[ps_num] =  r[i];
        ps_b_save[ps_num] = ps_b[i];
        ps_r2[ps_num] = r2[i];
        ps_num ++;
      }
      if((ssl_b[i]/ssl_max/r2[i]>=1) & (ssl_num<body.q)){
        ssl[ssl_num] = r[i];
        ssl_b_save[ssl_num] = ssl_b[i];
        ssl_r2[ssl_num] = r2[i];
        ssl_num ++;
      }
      if((ssr_b[i]/ssr_max/r2[i]>=1) & (ssr_num<body.q)){
        ssr[ssr_num] = r[i];
        ssr_b_save[ssr_num] = ssr_b[i];
        ssr_r2[ssr_num] = r2[i];
        ssr_num ++;
      }
      if((phi_b[i]/phi_max/r2[i]>=1) & (phi_num<body.q)){
        sins[phi_num] = r[i] * 2.;
        phi_b_save[phi_num] = phi_b[i];
        phi_r2[phi_num] = r2[i];
        phi_num ++;
      }
    }
    for(i=0;i<body.q;i++){
      if(sins[i]>P3I2)coss[i] = sins[i]-P3I2;
      else coss[i] = sins[i] + PI2;
    }
  }

  /* Estimate total scattering coefficient */
  /* 1/(4pi)*int_(4pi)g {ED4PI}; {PI} because of numerical Theta integration (sum/N*interval);
     {count * body.q} number of THETA samples; {pow(body.l,4)*ED4PI} originates from g definition;
     {PSDF_GAUSS_FAC} is a constant factor of the PSDF which is not taken into account while calculating
     probability tables; the same applies to {X***_FAC} terms; finally the {PI} or {P2I} factors
     originate from integration over PHI */
  if (!strcmp(body.type,"gauss")){
    tmpa[0] = PSDF_GAUSS_FAC(1.,body.a);
  }
  else if (!strcmp(body.type,"expo")){
    tmpa[0] = PSDF_EXP_FAC(1.,body.a);
  }
  else if (!strcmp(body.type,"karman")){
    tmpa[0] = PSDF_KARMAN_FAC(1.,body.a,body.kap);
  }
//  tmp = ED4PI*PI/(count * body.q)*pow(body.l,4)*ED4PI*tmpa[0];
//  g0[0] = gpp_0 = pp_sum * tmp * pow(XPPP_FAC(body.gam),2) * P2I;
//  g0[1] = gps_0 = ps_sum * tmp / body.gam * P2I;
//  g0[2] = gsp_0 = ps_sum * tmp * pow(XSPP_FAC(body.gam),2) * body.gam * PI;
//  g0[3] = gssl_0 = ssl_sum * tmp * PI;
//  g0[4] = gssr_0 = ssr_sum * tmp * PI;
// 2D
  /* The following line should be replaced */
  //tmp = ED4PI/2*PI/(count * body.q)*pow(body.l,3)*2*ED4PI*tmpa[0];
  /* The following line takes into account the definition of g_xx in 2D with the factor "k**3/(8PI)" the number of samples in the sum and the PSDF_FAC.
     There is no further integration over phi that brings in a 2PI or Pi factor and according to the definition of g0_xx with the 1/(2PI) in front of the integral
     there is no further PI from the integral boarders. */
  tmp = 1.0/(count * body.q)*pow((body.l/body.gam),3)*tmpa[0];
  g0[0] = gpp_0 = pp_sum * tmp;
  g0[1] = gps_0 = ps_sum * tmp;
  g0[2] = gsp_0 = ps_sum * tmp * pow(XSPP_FAC(body.gam),2) * body.gam;
  g0[3] = gssl_0 = ssl_sum * tmp ;
  g0[4] = gssr_0 = ssr_sum * tmp ;

//  g0[5] = pp_cos * tmp * pow(XPPP_FAC(body.gam),2) * P2I;
//  g0[6] = ps_cos * tmp / body.gam * P2I;
//  g0[7] = ps_cos * tmp * pow(XSPP_FAC(body.gam),2) * body.gam * PI;
//  g0[8] = ssl_cos * tmp * PI;
//  g0[9] = ssr_cos * tmp * PI;
// 2D
  g0[5] = pp_cos * tmp * body.gam  * pow(XPPP_FAC(body.gam),2);
  g0[6] = ps_cos * tmp;
  g0[7] = ps_cos * tmp * pow(XSPP_FAC(body.gam),2) * body.gam;
  g0[8] = ssl_cos * tmp;
  g0[9] = ssr_cos * tmp;

//  g0[10] = ((gssl_0+gssr_0+gsp_0) - (g0[8]+g0[9]) + g0[6]) / ((gpp_0 + gps_0 - g0[5]) * ((gssl_0+gssr_0+gsp_0) - (g0[8]+g0[9])) - g0[6] * g0[7]);
//  g0[11] = ((gpp_0 + gps_0) - g0[5] + g0[7]) / ((gpp_0 + gps_0 - g0[5]) * ((gssl_0+gssr_0+gsp_0) - (g0[8]+g0[9])) - g0[6] * g0[7]);
// 2D
  g0[10] = ((gssl_0+gsp_0) - (g0[8]) + g0[6]) / ((gpp_0 + gps_0 - g0[5]) * ((gssl_0+gsp_0) - (g0[8])) - g0[6] * g0[7]);
  g0[11] = ((gpp_0 + gps_0) - g0[5] + g0[7]) / ((gpp_0 + gps_0 - g0[5]) * ((gssl_0+gsp_0) - (g0[8])) - g0[6] * g0[7]);

  if (!strcmp(body.type,"isotrop")){
    g0[0] = gpp_0 = 1./body.a;   /* The parameter body.a (correlation length) is interpreted as the scattering mean free path of the isotropic scattering. Only g0[0] is meaningful, all other parameters are nonsense in the case of isotropic acoustic scattering. */
    g0[1] = gps_0 = g0[2] = gsp_0 = g0[3] = gssl_0 = g0[4] = gssr_0 = g0[5] = g0[6] = g0[7] = g0[8] = g0[9] = g0[10] = g0[11] = 0.;
  }
  /* Write total scattering coefficients to file */
  sprintf(tmpc,"%s.tot",body.ProbTabFileName);
  if((tot_file = fopen(tmpc,"w"))==NULL){
    printf("%s:%d Error creating file %s\n",__FILE__,__LINE__,tmpc);
    return -1;
  }

  fprintf(tot_file,"Total scattering coefficients for %s gam=%1.4e l=%1.4e a=%1.4e e=%1.4e ny=%1.4e q=%1.2e.\n",body.type,body.gam,body.l, body.a, body.e, body.ny, body.q);
  fprintf(tot_file,"gpp_0 = %1.6e\n",gpp_0);
  fprintf(tot_file,"gps_0 = %1.6e\n",gps_0);
  fprintf(tot_file,"gsp_0 = %1.6e\n",gsp_0);
  fprintf(tot_file,"gssl_0 = %1.6e\n",gssl_0);
  fprintf(tot_file,"gssr_0 = %1.6e\n",gssr_0);
  fprintf(tot_file,"<cos(theta)>_pp = %1.6e\n",g0[5]);
  fprintf(tot_file,"<cos(theta)>_ps = %1.6e\n",g0[6]);
  fprintf(tot_file,"<cos(theta)>_sp = %1.6e\n",g0[7]);
  fprintf(tot_file,"<cos(theta)>_ssl = %1.6e\n",g0[8]);
  fprintf(tot_file,"<cos(theta)>_ssr = %1.6e\n",g0[9]);
  fprintf(tot_file,"l^*_p = %e\n",g0[10]);
  fprintf(tot_file,"l^*_s = %e\n",g0[11]);

  fclose(tot_file);

  
  /* Write probability tables to files */
  sprintf(tmpc,"%s.dat",body.ProbTabFileName);
  if((pat_file = fopen(tmpc,"w"))==NULL){
    printf("%s:%d Error creating file %s\n",__FILE__,__LINE__,tmpc);
    return -1;
  }

  fprintf(pat_file,"Probability table for THETA angle of P-P, P-S, S-Sl, S-Sr scattering and for (cos(phi))^2 and (sin(phi))^2 (THETA pattern for S-P is the same as for P-S\n");

  for(i=0;i<body.q;i++){
    if(sins[i]>P3I2)coss[i] = sins[i]-P3I2;
    else coss[i] = sins[i] + PI2;
    fprintf(pat_file,"%1.12f\t%1.12f\t%1.12f\t%1.12f\t%1.12f\t%1.12f\n",pp[i],ps[i],ssl[i],ssr[i],coss[i],sins[i]);
  }
  fclose(pat_file);
  *(fltpt+0) = pp;
  *(fltpt+1) = ps;
  *(fltpt+2) = ssl;
  *(fltpt+3) = ssr;
  *(fltpt+4) = coss;
  *(fltpt+5) = sins;

  return 0;
}







int ReadScatPattern(bodyparm body,  float ** fltpt, double * g0)
{
  FILE * pat_file, *tot_file;
  float *pp, *ps, *ssl, *ssr, *sins, *coss;
  char dump[501];
  char tmpc[121];
  int count = 0, i;

  /* Read prabability tables */
  sprintf(tmpc,"%s.dat",body.ProbTabFileName);
  if((pat_file = fopen(tmpc,"r"))==NULL){
    printf("%s:%d Error opening file %s\n",__FILE__,__LINE__,tmpc);
    return -1;
  }
  pp = (float *) calloc((size_t) body.q, (size_t) sizeof(float));
  ps = (float *) calloc((size_t) body.q, (size_t) sizeof(float));
  ssl = (float *) calloc((size_t) body.q, (size_t) sizeof(float));
  ssr = (float *) calloc((size_t) body.q, (size_t) sizeof(float));
  coss = (float *) calloc((size_t) body.q, (size_t) sizeof(float));
  sins = (float *) calloc((size_t) body.q, (size_t) sizeof(float));

  if((pp == NULL) | (ps == NULL) | (ssl == NULL) | (ssr == NULL) | (coss == NULL) | (sins == NULL)){
    printf("%s:%d Error allocating memory.\n",__FILE__,__LINE__);
    return -2;
    }
  fgets(dump,500,pat_file);
  for(i=0;i<body.q;i++){
    count += fscanf(pat_file,"%f\t%f\t%f\t%f\t%f\t%f\n",&pp[i],&ps[i],&ssl[i],&ssr[i],&coss[i],&sins[i]);
  }
  fclose(pat_file);
  if(count!=6*body.q){
    printf("File %s: Number of entries did not match expected number.",tmpc);
    return -1;
  }
  *(fltpt+0) = pp;
  *(fltpt+1) = ps;
  *(fltpt+2) = ssl;
  *(fltpt+3) = ssr;
  *(fltpt+4) = coss;
  *(fltpt+5) = sins;

  /* Read total scattering coefficients */
  sprintf(tmpc,"%s.tot",body.ProbTabFileName);
  if((tot_file = fopen(tmpc,"r"))==NULL){
    printf("%s:%d Error opening file %s\n",__FILE__,__LINE__,tmpc);
    return -1;
  }
  fgets(dump,500,tot_file);
  count = 0;
  count += fscanf(tot_file,"gpp_0 = %le\n",&g0[0]);
  count += fscanf(tot_file,"gps_0 = %le\n",&g0[1]);
  count += fscanf(tot_file,"gsp_0 = %le\n",&g0[2]);
  count += fscanf(tot_file,"gssl_0 = %le\n",&g0[3]);
  count += fscanf(tot_file,"gssr_0 = %le\n",&g0[4]);
  count += fscanf(tot_file,"<cos(theta)>_pp = %le\n",&g0[5]);
  count += fscanf(tot_file,"<cos(theta)>_ps = %le\n",&g0[6]);
  count += fscanf(tot_file,"<cos(theta)>_sp = %le\n",&g0[7]);
  count += fscanf(tot_file,"<cos(theta)>_ssl = %le\n",&g0[8]);
  count += fscanf(tot_file,"<cos(theta)>_ssr = %le\n",&g0[9]);
  count += fscanf(tot_file,"l^*_p = %le\n",&g0[10]);
  count += fscanf(tot_file,"l^*_s = %le\n",&g0[11]);
  fclose(tot_file);
  if(count!=12){
    printf("File %s: Number of entries did not match expected number.",tmpc);
    return -1;
  }
  return 0;
}
