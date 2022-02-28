
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
#if ISMPI
  #include "mpi.h"
#endif
#include "edefinition.h"
#include "eutility.h"
#include "mk_scat_pattern.h"
#include "mk_source_pattern.h"
#include "scatter_2D.h"
#include "surface.h"
#include "layer.h"






int main(int ARGC, char *ARGV[])
{
  /* DECLARATION */
  simulationparameters simpa;       /* stores all parameters that define the simulation */
  particle pat;                     /* stores all attributes of the particle that is being propagated */  
  
  /* DEFINE GRIDS FOR STORING COMPONENTS */
  gridtype ***** P;
  gridtype ***** S;

  double eP[3],eSV[3],eSH[3],eL[3],eQ[3],eT[3];                                                            
  double d1,d2,p1,p2,distance;
  double tpl,tpq,tpt,tsl_sv,tsq_sv,tst_sv,tsl_sh,tsq_sh,tst_sh;
  
  gridtype ***** tmpgrd;             /* Stores the energy grid during reduction of MPI results */
  int status;                       /* return value of most of the functions */
  float *fltpt[6];
  double g0[12];                    /* gpp_0, gps_0, gsp_0, gssl_0 gssr_0, <cos(theta)>_pp, <cos(theta)>_ps(sp), <cos(theta)>_ssl, <cos(theta)>_ssr, l*_p, l*_s */
  char call[201];                   /* general purpose string */

  /* parameters describing the scattering coefficients */
  bodytable * scattab;                 /* This structure contains the scattering probability tables of the different bodies */
                                       /*total scattering coefficients (1/km)*/

  /* parameters describing the source radiation */
  float *pr_rad_t, *st_rad_t, *sf_rad_t, *cos_phi_rad, *sin_phi_rad;  /* probability tables for theta of P, S_theta and S_phi radiation and probability tables for (cos(2*phi)) and (sin(2*phi))^2 */
  double rl[3], rm[3], rn[3];          /* rotation matrix that rotates a vector from the local coordinates
				       of the source with n=(1 0 0) and s=(0 1 0) into the global system
				       with source orientation given by simpa.SourceOri */

  /* parameters describing the interface/surface reflections/conversions */
  float * layerP_t[10];
  float * layerS_t[10];
  float * surfP_t[4];
  float * surfS_t[5];

  /* variables used in the particle propagaion */
  register long count;             /* particle counter */
  register int tcount;             /* timestep counter */
  register int tfraccount;         /* counter of timeStep fractions */
  register long lnpart;
  int tmpa;
  double epsilon_now;             /* the value of e in the current grid*/
  double Qexp_now[2];
  register int stepnum;            /* number of time steps to simulate */
  double tmpd;
  int b = 0,oldb;                  /* number of body in which the particel currently is */

  /* parametes describing the scattering event */
  double the;                      /* scattering angel */
  double phi;                      /* angle betwenn polarization and scattering plane */

  /* parameters for particle counting */
  int x,y,z,theta,t;
  int temp_t;
  float percentage;  
  /* parameters for MPI */
  int my_rank;
  int psize;
  char my_node_name[101];
  int node_name_length = 101;
  gridtype * timesP;
  gridtype * timesS;
  int i,j,k,l;
  char tmpc[101];
  time_t starttime;                /* beginning of execution */
  time_t loopstart;                /* beginning of particle loop */
  double TimeStepPerGrid;
  /* Initialize MPI */
#if ISMPI
  starttime = time(0);
  DBG(INITIALIZE MPI);
  MPI_Init(&ARGC, &ARGV); 
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &psize);
  MPI_Get_processor_name(my_node_name, &node_name_length);
  srand((unsigned int) (starttime+my_rank+2));
#else
  starttime = time(0);
  srand((unsigned int) starttime);
#endif

  /* Read parameter file */
  status = ReadParFile(ARGC, ARGV, &simpa);
  if(status==0) {   
                    printf("***************************\n");
                    printf("Name of parameter file: %s\n",simpa.ParFileName);
                    printf("***************************\n");
                    printf("\n"); 
                }
  else return -1;

    
  status = ReadModelFiles(&simpa);

  sprintf(tmpc,"sims/%s_%d.log",simpa.dstring,my_rank);
  simpa.local_out = fopen(tmpc,"w");
  fprintf(simpa.local_out,"Process %d on machine %s with first random number %d\n",my_rank,my_node_name,rand());
#if ISMPI
  lnpart = (long)((double)simpa.npart/psize);
  if (my_rank == 0)lnpart = (long)(simpa.npart - lnpart*(psize - 1));
#else
  lnpart = simpa.npart;
#endif
  /* SCATTERING PATTERN */
  scattab = (bodytable *) malloc((size_t) simpa.Nbody * sizeof(bodytable));
  for(i=0;i<simpa.Nbody;i++){
    status = ReadScatPattern(simpa.body[i], fltpt, g0);
    if(status != 0){
      printf("Recalculating propability tables for scattering pattern.\n");
      status = MkScatPattern(simpa.body[i], fltpt, g0);
    }
    scattab[i].pp=fltpt[0];scattab[i].ps = fltpt[1]; scattab[i].sp = fltpt[1]; scattab[i].ssl = fltpt[2]; scattab[i].ssr = fltpt[3]; scattab[i].sp_phi = fltpt[4]; scattab[i].ssl_phi = fltpt[4]; scattab[i].ssr_phi = fltpt[5];
    scattab[i].gpp0=g0[0]; scattab[i].gps0=g0[1]; scattab[i].gsp0=g0[2]; scattab[i].gssl0=g0[3]; scattab[i].gssr0=g0[4];
    scattab[i].cos_pp=g0[5]; scattab[i].cos_ps=g0[6]; scattab[i].cos_sp=g0[7]; scattab[i].cos_ssl=g0[8]; scattab[i].cos_ssr=g0[9];
    scattab[i].l_p=g0[10]; scattab[i].l_s=g0[11];
  }

  /* Calculate auxillary parameters */
  for(i=0;i<simpa.Nbody;i++){
    simpa.body[i].pg[0] = (scattab[i].gpp0);                             /*acoustic: only scattering coefficient from p to p */
    simpa.body[i].pg[1] = (scattab[i].gssl0+scattab[i].gsp0);           /* scattering coefficient from s */
    simpa.body[i].tau[0] = 1./(simpa.body[i].pg[0]*simpa.body[i].v[0]);                  /* mean free time of P-particle */
    simpa.body[i].tau[1] = 1./(simpa.body[i].pg[1]*simpa.body[i].v[1]);                  /* mean free time of S-particle */
  }

  /* check input parameters */
  status = checkSimulationParameters(&simpa);

  /* Calculate other auxillary parameters */
  for(i=0;i<simpa.Nbody;i++){
    simpa.body[i].Qex[0] = exp(-simpa.frc * (simpa.body[i].Qi[0]-simpa.Qref) * (simpa.TimeStep/simpa.body[b].TimeStepFrac));       /* 2*\pi*f*Q_p^{-1}/\alpha (exponent from 5.1, Sato p.110) (factor to multiplyed on the energy each time interval */
    simpa.body[i].Qex[1] = exp(-simpa.frc * (simpa.body[i].Qi[1]-simpa.Qref) * (simpa.TimeStep/simpa.body[b].TimeStepFrac));       /* 2*\pi*f*Q_s^{-1}/\alpha Sato p.110 */
    simpa.body[i].pgexp[0] = exp(-simpa.body[i].pg[0]*simpa.body[i].path[0]);            /* probability of P-wave to be scattered in one path */
    simpa.body[i].pgexp[1] = exp(-simpa.body[i].pg[1]*simpa.body[i].path[1]);            /* probability of S-wave to be scattered in one path */
  }
  simpa.spart = 1.5 * pow(simpa.body[simpa.sourcebody].gam,5)/(1.5 * pow(simpa.body[simpa.sourcebody].gam,5) + 1);     /* fraction energy radiated as S-wave assuming source is in body 1 */
 /* ALLOCATE GRIDS */
  P     = GridAlloc(simpa); 
  S     = GridAlloc(simpa); 
  DBG(After the GridAlloc); 
#if ISMPI
  if(psize>1){   /* Grid for reduction of results */
    tmpgrd = GridAlloc(simpa);
  }
#endif

  /* SOURCE RADIATION PATTERN */
  status = source_rot_mat(simpa,rl,rm,rn);
  status =  ReadRadPat(simpa, fltpt);
  if(status){
    printf("\n"); 
    printf("***************************\n");
    printf("Recalculating propability tables for source radiation pattern.\n");
    printf("***************************\n");
    printf("\n");  
    status = MkSourcePattern(simpa, fltpt);
  }
  pr_rad_t = fltpt[0];  st_rad_t = fltpt[1]; sf_rad_t = fltpt[2]; cos_phi_rad = fltpt[3]; sin_phi_rad = fltpt[4];

  /* SURFACE COEFFICIENTS */
  status = ReadSurface(simpa,surfS_t, surfP_t);
  if(status != 0){
    printf("Recalculating reflection/conversion/propability tables for the interface.\n");
    status = MkSurface(simpa, surfS_t, surfP_t);
  }

  /* INTERFACE COEFFICIENTS */
  status = ReadInterface(simpa,layerS_t, layerP_t);
  if(status != 0){
    printf("Recalculating reflection/conversion/propability tables for the interface.\n");
    status = MkInterface(simpa, layerS_t, layerP_t);
  }

  /* Calculate parameters of the modeling calculations */
  stepnum = (int)ceil((float)simpa.SimLen/simpa.TimeStep);
  TimeStepPerGrid = round(simpa.GridSpacing[4]/simpa.TimeStep);

  /* ## PARTICLE LOOP ################################################# */
  printf("***************************\n");
  DBG(Starting simulation);
  printf("***************************\n");
  
  printf("\n");
  printf("-----------------------------------\n");
  printf("Source Type: %d\n",simpa.SourceType);
  printf("Source position x: %f\n",simpa.SourcePos[0]);
  printf("Source position y: %f\n",simpa.SourcePos[1]);
  printf("Source position z: %f\n",simpa.SourcePos[2]);
  printf("-----------------------------------\n");
  printf("GridSize x: %d\n",simpa.GridSize[0]);
  printf("GridSize y: %d\n",simpa.GridSize[1]);
  printf("GridSize z: %d\n",simpa.GridSize[2]);
  printf("GridSize theta: %d\n",simpa.GridSize[3]);
  printf("GridSize t: %d\n",simpa.GridSize[4]);
  printf("-----------------------------------\n");
  printf("GridSpacing x: %f\n",simpa.GridSpacing[0]);
  printf("GridSpacing y: %f\n",simpa.GridSpacing[1]);
  printf("GridSpacing z: %f\n",simpa.GridSpacing[2]);
  printf("GridSpacing theta: %f\n",simpa.GridSpacing[3]);
  printf("GridSpacing t: %f\n",simpa.GridSpacing[4]);
  printf("-----------------------------------\n");
  printf("\n");

#if ISMPI
  if (my_rank == 0)loopstart = time(0);
#else
  loopstart = time(0);
#endif

  /* Initialize particle */
  for(count=0;count<lnpart;count++){
    percentage = count/simpa.npart;

   /* printf("\b");*/

   /* printf("\rParticle number %ld of %ld",count,simpa.npart);*/
    tcount = 0;
    pat.time=0.0;
    pat.pos[0] = simpa.SourcePos[0];
    pat.pos[1] = simpa.SourcePos[1];
    pat.pos[2] = simpa.SourcePos[2];
    tfraccount = simpa.body[simpa.sourcebody].TimeStepFrac;
    b = simpa.sourcebody;

    /*************************************************/
    /* Isotropic P-source */
    if(!simpa.SourceType){
        pat.mode = 0;
        pat.stokes[0] = 1.0;
        pat.stokes[1] = pat.stokes[2] = pat.stokes[3] = pat.stokes[4] = 0.0;
        pat.dir[0]=PI2;//acos(1.0-2.0*rand()/(RAND_MAX));
        pat.dir[1]=P2I*rand()/RAND_MAX;
    /*pat.dir[1]=0;
    pat.dir[0]=1.4;*/
        pat.vec[0] = simpa.body[b].path[0] * sin(pat.dir[0]) * cos(pat.dir[1]);
        pat.vec[1] = simpa.body[b].path[0] * sin(pat.dir[0]) * sin(pat.dir[1]);
        pat.vec[2] = 0;//- simpa.body[b].path[0] * cos(pat.dir[0]);
    }
    /*************************************************/

    /*************************************************/
    /* S-source */
if(simpa.SourceType){
    pat.mode = 1;
    pat.dir[0]=PI2;//acos(1.0-2.0*rand()/(RAND_MAX));
    pat.dir[1]=P2I*rand()/RAND_MAX; 
/* pat.dir[1]=0;
 pat.dir[0]=0.5;*/
    if((double) rand()/RAND_MAX < -0.1){      /* S theta particle */
      the = pat.dir[0]+PI2;
      phi = pat.dir[1];
      pat.pol = PI-asin(sin(the)*sin(pat.dir[1]-phi));
      status = stokesRot(pat.pol, 1, 0.0, 1.0, 0.0, 0.0, 0.0, &pat.stokes[0], &pat.stokes[1], &pat.stokes[2], &pat.stokes[3], &pat.stokes[4]);
    }
    else{                      /* S phi particle */
      the = PI2;
      phi = pat.dir[1]+PI2;
      pat.pol = PI-asin(sin(the)*sin(pat.dir[1]-phi));
      status = stokesRot(pat.pol, 1, 0.0, 1.0, 0.0, 0.0, 0.0, &pat.stokes[0], &pat.stokes[1], &pat.stokes[2], &pat.stokes[3], &pat.stokes[4]);
    }
    pat.vec[0] = simpa.body[b].path[1] * sin(pat.dir[0]) * cos(pat.dir[1]);
    pat.vec[1] = simpa.body[b].path[1] * sin(pat.dir[0]) * sin(pat.dir[1]);
    pat.vec[2] = 0;//- simpa.body[b].path[0] * cos(pat.dir[0]);
}
    /*************************************************/
 fprintf(simpa.local_out,"%1.4f\t%1.4f\n",pat.dir[0],pat.dir[1]);

    /*************************************************/
    /* Nonisotropic source */
if(0)
{
tmpd = (double) rand()/RAND_MAX;

    if (0){//(((double) rand()/RAND_MAX) > simpa.spart){    /* P particle */
    //if (tmpd > simpa.spart){    /* P particle */
//printf("H\n",tmpd);
      pat.mode = 0;
      tmpa = (int)ceil((double)simpa.q * rand()/RAND_MAX)-1;
      pat.dir[0] = pr_rad_t[tmpa];
      tmpa = (int)ceil((double)simpa.q * rand()/RAND_MAX)-1;
      pat.dir[1] = sin_phi_rad[tmpa];
//printf("P-part the:%e phi:%e -->",pat.dir[0],pat.dir[1]);
      status = SourceRot(&pat.dir[0],&pat.dir[1],rl,rm,rn);
//printf("%e %e\n",pat.dir[0],pat.dir[1]);
      pat.stokes[0] = 1.0;
      pat.stokes[1] = pat.stokes[2] = pat.stokes[3] = pat.stokes[4] = 0.0;
      pat.vec[0] = simpa.body[b].path[0] * sin(pat.dir[0]) * cos(pat.dir[1]);
      pat.vec[1] = simpa.body[b].path[0] * sin(pat.dir[0]) * sin(pat.dir[1]);
      pat.vec[2] = - simpa.body[b].path[0] * cos(pat.dir[0]);
//printf("%e %e %e\n",pat.vec[0],pat.vec[1],pat.vec[2]);
    }
    else{                /* S particle */
      pat.mode = 1;
      if((double) rand()/RAND_MAX < 0.1666666666666667){      /* S theta particle */
        tmpa = (int)ceil((double)simpa.q * rand()/RAND_MAX)-1;
        pat.dir[0] = st_rad_t[tmpa];
        tmpa = (int)ceil((double)simpa.q * rand()/RAND_MAX)-1;
        pat.dir[1] = sin_phi_rad[tmpa];
        the = pat.dir[0]+PI2;
        phi = pat.dir[1];
        status = SourceRot(&pat.dir[0],&pat.dir[1],rl,rm,rn);
        status = SourceRot(&the,&phi,rl,rm,rn);
        pat.pol = PI-asin(sin(the)*sin(pat.dir[1]-phi));
        status = stokesRot(pat.pol, 1, 0.0, 1.0, 0.0, 0.0, 0.0, &pat.stokes[0], &pat.stokes[1], &pat.stokes[2], &pat.stokes[3], &pat.stokes[4]);
      }
      else{                      /* S phi particle */
        tmpa = (int)ceil((double)simpa.q * rand()/RAND_MAX)-1;
        pat.dir[0] = sf_rad_t[tmpa];
        tmpa = (int)ceil((double)simpa.q * rand()/RAND_MAX)-1;
        pat.dir[1] = cos_phi_rad[tmpa];
        the = PI2;
        phi = pat.dir[1]+PI2;
        status = SourceRot(&pat.dir[0],&pat.dir[1],rl,rm,rn);
        status = SourceRot(&the,&phi,rl,rm,rn);
        pat.pol = PI-asin(sin(the)*sin(pat.dir[1]-phi));
        status = stokesRot(pat.pol, 1, 0.0, 1.0, 0.0, 0.0, 0.0, &pat.stokes[0], &pat.stokes[1], &pat.stokes[2], &pat.stokes[3], &pat.stokes[4]);
      }
      pat.vec[0] = simpa.body[b].path[1] * sin(pat.dir[0]) * cos(pat.dir[1]);
      pat.vec[1] = simpa.body[b].path[1] * sin(pat.dir[0]) * sin(pat.dir[1]);
      pat.vec[2] = - simpa.body[b].path[1] * cos(pat.dir[0]);
    }
}
    /*************************************************/
    //temp_t = -1;
    while(pat.time<simpa.SimLen){
      /************************************************/
      /* count */
      if(tfraccount == simpa.body[b].TimeStepFrac){
        if((0<=pat.pos[0]) & (0<=pat.pos[1]) & (0<=pat.pos[2])){
          x = (int)(pat.pos[0]/simpa.GridSpacing[0]);
          y = (int)(pat.pos[1]/simpa.GridSpacing[1]);
          z = (int)(pat.pos[2]/simpa.GridSpacing[2]);
          while (pat.dir[1]< 0) { pat.dir[1] = pat.dir[1] + P2I;}
          theta = (int)((pat.dir[1]/PI*180.00)/simpa.GridSpacing[3]);
          while (theta >= simpa.GridSize[3]) {theta = theta - simpa.GridSize[3];}
          t = (int)(tcount/TimeStepPerGrid);
          /*t = (int)(pat.time/simpa.GridSpacing[3]);*/
/* 
          distance = sqrt( pow(simpa.SourcePos[0]-pat.pos[0],2) + pow(simpa.SourcePos[1]-pat.pos[1],2) + pow(simpa.SourcePos[2]-pat.pos[2],2) );  
          p1 = PI - acos((simpa.SourcePos[2]-pat.pos[2])/distance);
          p2 = atan2(simpa.SourcePos[1]-pat.pos[1],simpa.SourcePos[0]-pat.pos[0]);  
            if (p2<0){
            p2 = p2 + 2*PI;
            }
          d1 = pat.dir[0];
          d2 = pat.dir[1];  

          eP[0]  = sin(d1)*cos(d2);
          eP[1]  = sin(d1)*sin(d2);
          eP[2]  = cos(d1);
          eSV[0] = cos(d1)*cos(d2);
          eSV[1] = cos(d1)*sin(d2);
          eSV[2] = -sin(d1);
          eSH[0] = -sin(d2);
          eSH[1] = cos(d2);
          eSH[2] = 0;
          eL[0]  = sin(p1)*cos(p2);
          eL[1]  = sin(p1)*sin(p2);
          eL[2]  = cos(p1);
     	  eQ[0]  = cos(p1)*cos(p2);
          eQ[1]  = cos(p1)*sin(p2);
          eQ[2]  = -sin(p1);
          eT[0]  = -sin(p2);
          eT[1]  = cos(p2);
          eT[2]  = 0;
	
          tpl    = pow(((eP[0]  * eL[0]) + (eP[1]  * eL[1]) + (eP[2]  * eL[2])),2);
          tpq    = pow(((eP[0]  * eQ[0]) + (eP[1]  * eQ[1]) + (eP[2]  * eQ[2])),2);
          tpt    = pow(((eP[0]  * eT[0]) + (eP[1]  * eT[1]) + (eP[2]  * eT[2])),2);
          tsl_sv = pow(((eSV[0] * eL[0]) + (eSV[1] * eL[1]) + (eSV[2] * eL[2])),2);
          tsq_sv = pow(((eSV[0] * eQ[0]) + (eSV[1] * eQ[1]) + (eSV[2] * eQ[2])),2);
          tst_sv = pow(((eSV[0] * eT[0]) + (eSV[1] * eT[1]) + (eSV[2] * eT[2])),2);
          tsl_sh = pow(((eSH[0] * eL[0]) + (eSH[1] * eL[1]) + (eSH[2] * eL[2])),2);
          tsq_sh = pow(((eSH[0] * eQ[0]) + (eSH[1] * eQ[1]) + (eSH[2] * eQ[2])),2);
          tst_sh = pow(((eSH[0] * eT[0]) + (eSH[1] * eT[1]) + (eSH[2] * eT[2])),2);
  */     
       // printf("x=%d ,y=%d, z=%d",x,y,z); 
	 if((x<simpa.GridSize[0]) & (y<simpa.GridSize[1]) & (z<simpa.GridSize[2]-REC_FLUX) & (t<simpa.GridSize[4])){  /* -REC_FLUX in z condition is for flux recording through the layer */
          //if((int)(t - temp_t!=0))
        //{
            //temp_t = t ;
          if(pat.mode)
		{
		/* S TOTAL*/
		S[x][y][z][theta][t]  += (gridtype)  ((pat.stokes[1]+pat.stokes[2])*simpa.emodulo);
        }
      else{
		/* P TOTAL */
		P[x][y][z][theta][t]  += (gridtype)  (pat.stokes[0]*simpa.emodulo);
         /* if(P[x][y][z][t]!=0) printf("P is %llu and  %llu,x= %d,y= %d, Z=%d t=%d pat.time = %.15f; tcount = %d, %.15f \n",P[x][y][z][t],simpa.emodulo,x,y,z,t,pat.time,tcount,(tcount/TimeStepPerGrid));*/
	//	printf("pat.stokes[0] = %lf, simpa.emodulo= %llu and P[x][y][z][theta][t]=%llu \n",pat.stokes[0], simpa.emodulo, P[x][y][z][theta][t]);
      }
		//}
          }
          	}
        
        tfraccount = 0;
        tcount ++;
      }
      /************************************************/

      /* propagation */
      if(b==1){  /* bend the ray if in the mantle */
        tmpd = 1./300*(pat.pos[2]-simpa.body[0].z[1]);  /* velocity increase compared to Moho */
        pat.dir[0] -= sin(PI-pat.dir[0])/(simpa.body[1].v[pat.mode]+tmpd) / 300 * simpa.body[1].path[pat.mode];
        pat.vec[0] = simpa.body[b].path[pat.mode] * sin(pat.dir[0]) * cos(pat.dir[1]);
        pat.vec[1] = simpa.body[b].path[pat.mode] * sin(pat.dir[0]) * sin(pat.dir[1]);
        pat.vec[2] = - simpa.body[b].path[pat.mode] * cos(pat.dir[0]);
        pat.time -= (tmpd*simpa.body[b].path[pat.mode])/(pow(simpa.body[b].v[pat.mode],2) + simpa.body[b].v[pat.mode]*tmpd);
      }
      pat.pos[0] += pat.vec[0];
      pat.pos[1] += pat.vec[1];
      pat.pos[2] += pat.vec[2];
      if(0 > pat.pos[0]){
            pat.pos[0] = - pat.pos[0];
            if(pat.dir[0] >=  PI){
                pat.dir[1] = 3*PI - pat.dir[1];
            }
            if(pat.dir[0] <  PI){
                pat.dir[1] = PI - pat.dir[1];
            }
        }
      if(4 < pat.pos[0]){
            pat.pos[0] = 8 - pat.pos[0];
            if(pat.dir[1] >=  PI){
                pat.dir[1] = 3*PI - pat.dir[1];
            }
            if(pat.dir[1] <  PI){
                pat.dir[1] = PI - pat.dir[1];
            }
        }
      if(0 > pat.pos[1]){
            pat.pos[1] = - pat.pos[1];
            pat.dir[1] = 2*PI - pat.dir[1];
        }
      if(5 < pat.pos[1]){
            pat.pos[1] = 10 - pat.pos[1];
            pat.dir[1] = 2*PI - pat.dir[1];
        }
      pat.vec[0] = simpa.body[b].path[pat.mode] * sin(pat.dir[0]) * cos(pat.dir[1]);
      pat.vec[1] = simpa.body[b].path[pat.mode] * sin(pat.dir[0]) * sin(pat.dir[1]);
  
      pat.time += simpa.TimeStep/simpa.body[b].TimeStepFrac;
      tfraccount ++;
      /* Intrinsic attenuation */
    /* get the exact value of Qex in current grid */
      if(((x<simpa.GridSize[0]) & (x > -1)) & ((y<simpa.GridSize[1]) & (y > -1)) & ((z<simpa.GridSize[2])  & (z > -1)))
        {
            Qexp_now[0]=exp(-simpa.frc * (simpa.body[b].anQp[x][y][z]-simpa.Qref) * (simpa.TimeStep/simpa.body[b].TimeStepFrac));
            Qexp_now[1]=exp(-simpa.frc * (simpa.body[b].anQs[x][y][z]-simpa.Qref) * (simpa.TimeStep/simpa.body[b].TimeStepFrac));
        }
        else
        {
            Qexp_now[0]=simpa.body[b].Qex[0];
            Qexp_now[1]=simpa.body[b].Qex[1];
        }
        
        
      if (pat.mode){
        pat.stokes[1] *= Qexp_now[pat.mode];
        pat.stokes[2] *= Qexp_now[pat.mode];
      }
      else pat.stokes[0] *= Qexp_now[pat.mode];
      /************************************************/
      /* check geometry */
#if SURFACE
      if (pat.pos[2]<0){
        status =surface(&pat,simpa.Nangle,simpa.body[0],surfP_t,surfS_t,simpa.body[b].path);
      }
#endif
#if INTERFACE
      status = 0;
      if (pat.pos[2] >= simpa.body[0].z[1]){     /* particle is below the interface*/
        if (pat.pos[2]-pat.vec[2] < simpa.body[0].z[1]){    /* and was above before */
	  oldb = b;
          b = layer(&pat,simpa,layerP_t,layerS_t,b,1);
          status = 1;   /* mark whether particle encountered the interface */
#if REC_FLUX
	  if(oldb!=b){  /* transmission occured */
	    if((0<=pat.pos[0]) & (0<=pat.pos[1]) & (0<=pat.pos[2])){
	      x = (int)(pat.pos[0]/simpa.GridSpacing[0]);y = (int)(pat.pos[1]/simpa.GridSpacing[1]);theta = (int)((pat.dir[1]/PI*180.00)/simpa.GridSpacing[3]);t = (int)(pat.time/simpa.GridSpacing[4]);
	      if((x<simpa.GridSize[0]) & (y<simpa.GridSize[1]) & (t<simpa.GridSize[4])){
		if(pat.mode) S[x][y][simpa.GridSize[2]-2][theta][t] += (gridtype)((pat.stokes[1]+pat.stokes[2])*simpa.emodulo);
		else P[x][y][simpa.GridSize[2]-2][theta][t] += (gridtype)(pat.stokes[0]*simpa.emodulo);
	      }
	    }
	  }
#endif
        }
      }

      if ((pat.pos[2] <= simpa.body[0].z[1]) & (status == 0)){    /* particle is above the interface and did not pass through it before */
        if (pat.pos[2]-pat.vec[2] > simpa.body[0].z[1]){    /* and was below before */ 
	  oldb = b;
          if(whichbody(&pat, &simpa) & 0) b = 2; else b = 0;  /* set the body of expected appeareance */  /*warum steht da eine 0 in der Bedingung 15.6.07*/
          b = layer(&pat,simpa,layerP_t,layerS_t,1,b);
#if REC_FLUX
	  if(oldb!=b){  /* transmission occured */
	    if((0<=pat.pos[0]) & (0<=pat.pos[1]) & (0<=pat.pos[2])){
	      x = (int)(pat.pos[0]/simpa.GridSpacing[0]);y = (int)(pat.pos[1]/simpa.GridSpacing[1]);theta = (int)((pat.dir[1]/PI*180.00)/simpa.GridSpacing[3]);t = (int)(pat.time/simpa.GridSpacing[4]);
	      if((x<simpa.GridSize[0]) & (y<simpa.GridSize[1]) & (t<simpa.GridSize[4])){
		if(pat.mode) S[x][y][simpa.GridSize[2]-1][theta][t] += (gridtype)((pat.stokes[1]+pat.stokes[2])*simpa.emodulo);
		else P[x][y][simpa.GridSize[2]-1][theta][t] += (gridtype)(pat.stokes[0]*simpa.emodulo);
	      }
	    }
	  }
#endif
        }
      }
#endif

      /* geometry */
#if BODY>0
      if(whichbody(&pat, &simpa) & pat.pos[2]<simpa.body[0].z[1]){
	if (b!=2){
	  pat.vec[0] *= simpa.body[2].path[pat.mode]/simpa.body[b].path[pat.mode];
	  pat.vec[1] *= simpa.body[2].path[pat.mode]/simpa.body[b].path[pat.mode];
	  pat.vec[2] *= simpa.body[2].path[pat.mode]/simpa.body[b].path[pat.mode];
	  tfraccount = (int)floor(tfraccount * simpa.body[2].TimeStepFrac / simpa.body[b].TimeStepFrac);
	  b = 2;
	}
      }
      else{
	if (b==2){
	  pat.vec[0] *= simpa.body[0].path[pat.mode]/simpa.body[b].path[pat.mode];
	  pat.vec[1] *= simpa.body[0].path[pat.mode]/simpa.body[b].path[pat.mode];
	  pat.vec[2] *= simpa.body[0].path[pat.mode]/simpa.body[b].path[pat.mode];
	  tfraccount = (int)floor(tfraccount * simpa.body[0].TimeStepFrac / simpa.body[b].TimeStepFrac);
	  b = 0;
        }
      }
#endif
        
        /************************************************/
        /* get the exact value of e in current grid */
        if(((x<simpa.GridSize[0]) & (x > -1)) & ((y<simpa.GridSize[1]) & (y > -1)) & ((z<simpa.GridSize[2])  & (z > -1)))
        {
            epsilon_now=simpa.body[b].ane[x][y][z];
        }
        else
        {
            epsilon_now=simpa.body[b].e;
        }
        /************************************************/
        /* check for scattering */
        if (epsilon_now>0){
        tmpd = (double)rand()/RAND_MAX;
        if(!pat.mode){
            if(pow(simpa.body[b].pgexp[0],pow(epsilon_now,2))<tmpd){
                status = scatterP(simpa.body[b], &pat, scattab[b].gpp0, scattab[b].pp, scattab[b].ps);
                pat.vec[0] = simpa.body[b].path[pat.mode] * sin(pat.dir[0]) * cos(pat.dir[1]);
                pat.vec[1] = simpa.body[b].path[pat.mode] * sin(pat.dir[0]) * sin(pat.dir[1]);
                pat.vec[2] = 0;//simpa.body[b].path[pat.mode] * cos(pat.dir[0]); 2D
                }
            }
            else{
                if(pow(simpa.body[b].pgexp[1],pow(epsilon_now,2))<tmpd){
                    status = scatterS(simpa.body[b], &pat, scattab[b].gsp0, scattab[b].gssl0, scattab[b].sp, scattab[b].sp_phi, scattab[b].ssl, scattab[b].ssl_phi, scattab[b].ssr, scattab[b].ssr_phi);
                    pat.vec[0] = simpa.body[b].path[pat.mode] * sin(pat.dir[0]) * cos(pat.dir[1]);
                    pat.vec[1] = simpa.body[b].path[pat.mode] * sin(pat.dir[0]) * sin(pat.dir[1]);
                    pat.vec[2] = 0;//- simpa.body[b].path[pat.mode] * cos(pat.dir[0]); 2D
                }
            }
	}
    }  /* end of propagation loop */
  }    /* end of particle loop */
printf("\n");
printf("********************************************\n");
printf("END OF PARTICLE LOOP\n");
printf("********************************************\n");
printf("\n");
#if ISMPI
  if (my_rank == 0)simpa.looptime = time(0) - loopstart;
#else
  simpa.looptime = time(0) - loopstart;
#endif

  /* reduce the MPI results */
#if ISMPI
  fprintf(simpa.local_out,"PROCESS %d beginning with grid reduction\n",my_rank);
  if(psize>1){
    DBG(Reduce arrays);
    for(i=0;i<simpa.GridSize[0];i++){
      for(j=0;j<simpa.GridSize[1];j++){
        for(k=0;k<simpa.GridSize[2];k++){
            for(l=0;l<simpa.GridSize[3];l++){
          timesP = (gridtype *) calloc((size_t) simpa.GridSize[4], sizeof(gridtype));
          timesS = (gridtype *) calloc((size_t) simpa.GridSize[4], sizeof(gridtype));
     
          MPI_Reduce((void *)P[i][j][k][l],(void *)timesP ,simpa.GridSize[4],MPI_GRIDTYPE,MPI_SUM,0,MPI_COMM_WORLD);
	  MPI_Reduce((void *)S[i][j][k][l],(void *)timesS ,simpa.GridSize[4],MPI_GRIDTYPE,MPI_SUM,0,MPI_COMM_WORLD);
	 
          if(my_rank==0){
            free(P[i][j][k][l]);
	    free(S[i][j][k][l]);
          
            P[i][j][k][l] = timesP;
	    S[i][j][k][l] = timesS;
            }
          }
        }
      }
    }
  } 
  fprintf(simpa.local_out,"PROCESS %d finished grid reduction\n",my_rank);
#endif
  fclose(simpa.local_out);

  /* WRITE OUTPUT */
#if ISMPI
  if(my_rank==0) status = WriteOutput(simpa,P,S,scattab);
#else
  status = WriteOutput(simpa,P,S,scattab);
#endif

  /* FREE GRIDS */
  status = GridFree(simpa,P);
  status = GridFree(simpa,S);

  
    printf("\n");
    printf("***********************************\n");
    printf("**SIMULATION FINISHED SUCCESFULLY**\n");
    printf("***********************************\n");
    printf("\n");
#if ISMPI
  MPI_Finalize();
#endif
  sprintf(call, "cp %s sims",simpa.ParFileName);
  system(call);
  return 0;
}
