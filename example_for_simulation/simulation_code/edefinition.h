#if defined(DEFINITIONS_H)
#else
#define DEFINITIONS_H

/*  coordinate system ***********************
 *                                          *
 *             X,North                      *
 *             /                            *
 *            /___ Y,East                   *
 *            |                             *
 *            |                             *
 *           Z,Depth                        *
 *                                          *
 * phi: is measured from north towards east *
 * the: is measured from up towards down    *   
 ********************************************/


#define PI 3.14159265358979
#define PI2 1.57079632679490   /* pi/2 */
#define PI4 0.78539816339745   /* pi/4 */
#define P2I 6.28318530717959   /* 2*pi */
#define P4I 12.56637061435917  /* 4*pi */
#define P3I2 4.71238898038469  /* 3/2*pi */
#define P7I4 5.49778714378214  /* 7/4*pi */
#define ED4PI 0.07957747154595 /* 1/(4*pi) */

#define DB 1
#define DBG(x) if(DB==1) fprintf(stderr,"%s:%d No error before:" #x "\n",__FILE__,__LINE__)
#define MAX(A,B)((A)>(B)?(A):(B))
#define MIN(A,B)((A)<(B)?(A):(B))

#define ISMPI 0       /* switch to toggle on and off the mpi related code. 0 is off */
#define BODY 0        /* switch to toggle shape of inclusion: 0 no inclusion
                         1 cylinderring around source with distance body[2].x[0] < r^2 < body[2].x[1]
                         2 box */
#define SURFACE 1     /* switch to ignore the surface conversion and reflection
                         0 ignore the surface
                         1 take the surface into account */
#define INTERFACE 1   /* switch to ignore the interface to halfspace
                         0 ignore the interface
                         1 take it into account */
#define REC_FLUX 0    /* switch to switch on and off recording of flux through interface in the lowest two layers
                         2 to use the two lowest layers to record the flux through the interface N-2 is downgoing N-1 is upgoing flux
                         0 to switch recording off*/

#include <time.h>


/* Redefinition of the function 'acos' to account for numerical uncertainties that can cause the argument to be out off [-1,1] */
#define ACOS2(I)((I)<=-1.?(3.14159265358979):((I)>=1.?(0.):(acos(I))))


typedef unsigned long gridtype;
#define MPI_GRIDTYPE MPI_UNSIGNED_LONG


typedef struct
{
  int tm_sec;
  int tm_min;
  int tm_hour;
  int tm_mday;
  int tm_mon;
  int tm_year;
  int tm_wday;
  int tm_yday;
  int tm_isdst;
}timestruct;

typedef struct
{
  double x[2];                   /* minimum, maximum x-coordinate */
  double y[2];                   /* minimum, maximum x-coordinate */
  double z[2];                   /* minimum, maximum x-coordinate */
  double v[2];                   /* P-wave velocity, S-wave velocity */
  double gam;
  char type[20];                 /* Type of ACF (gauss, exponential,...) */
  double rho;                    /* Density */
  double l;                      /* S-wave wavenumber */
  double a;                      /* correlation distance */
  double e;                      /* Fractional fluctuation */
  double ***ane;                 /* Array for anomaly epsilon in recording region */
  double ny;                     /* Birch ny */
  double kap;                    /* von Karman kappa*/
  double Qi[2];                  /* Intrinsic attenuation 1/Q_p and 1/Q_s */
  double *** anQp;                /* Array for anomaly 1/Q_p in recording region */
  double *** anQs;                /* Array for anomaly 1/Q_s in recording region */
  double Qex[2];                 /* axillary parameter describing P and S absobtion */
  double q;                      /* Number of intervals in probability table */
  double pg[2];                  /* Inverse of free path of P- and S-waves */
  double pgexp[2];               /* probabilities of P- and S-waves to be scattered on one path */
  double tau[2];                 /* mean free times  of P- and S-waves */
  double path[2];                /* length of pathes in one time step for P- and S-energy */
  int TimeStepFrac;              /* Fraction of simulationparameters.TimeStep in current body */
  char ProbTabFileName[100];     /* Name of file contianing Probability tables for scattering coefficients */
}bodyparm;


typedef struct
{
  char tag[300];                  /* identification tag perceeding the output files */
  int outputtype;                /* type of output file (0 for 4-dimensional data and 1 for seismic section (time-distance) of surface layer */
  long npart;                    /* number of particles to shoot */
  char ParFileName[300];          /* name of parameter file */
  time_t SimulationTime;         /* random number identifying the simulation */
  time_t looptime;               /* seconds spent in the particle loop */
  double SourcePos[3];           /* x, y, z coordinates of source position */
  int SourceType;                /* Type of source, 0 for P and 1 for S-Source */
  double SourceOri[3];           /* Source orientation, (strike, dip and rake) of one nodal plane */
  int sourcebody;                /* Number of body that contains the source */
  double frc;                    /* angular frequency */
  double spart;                  /* number of S particels for number of S plus number of P particle for double couple radiation */
  gridtype emodulo;                   /* = 1/EnergyResolution (used to model damping or other things that cause the weight of a particle to be non integer) */
  float Qref;                    /* reference value of (intrinsic) attenuation 1/Q (a body with this Q appears to have no attenuation in the simulation) */
  int GridSize[5];               /* Size of the grid (number of grid points in x, y, z direction ,theta and time) */
  double GridSpacing[5];          /* Spacing of grid points in km (x, y, z) and theta (degree) s for time */
  double TimeStep;                /* Time step size (interval) of simulation */
  float SimLen;                  /* Length of simulation in s */
  double q;                      /* Number of intervals in probability table */
  int Nangle;                    /* Number of partitions of the [0 PI/2] at which the reflection/transmission coefficients are evaluated */
  char dstring[100];             /* Tag_data_time string of the output files */
  FILE * local_out;              /* File pointer to individual output files of each process */
  int Nbody;                     /* Number of bodies with different scattering properties */
  bodyparm * body;
  char Espsilon_file[300];       /* name of file sedcribing the spatial distribution od epsilon */
  char Qp_file[300];             /* name of file sedcribing the spatial distribution od 1/Q_p */
  char Qs_file[300];             /* name of file sedcribing the spatial distribution od 1/Q_s */
  char InterfaceName[100];       /* Name of file contianing reflection and transmission coefficients of the interface */
  char SurfaceName[100];         /* Name of file contianing reflection and conversion coefficients of the surface */
}simulationparameters;

typedef struct
{
  float *pp;
  float *ps;
  float *sp;
  float *ssl;
  float *ssr;
  float *sp_phi;
  float *ssl_phi;
  float *ssr_phi;
  double gpp0;
  double gps0;
  double gsp0;
  double gssl0;
  double gssr0;
  double cos_pp;
  double cos_ps;
  double cos_sp;
  double cos_ssl;
  double cos_ssr;
  double l_p;
  double l_s;
}bodytable;



typedef struct
{
  double pos[3];       /* cartesian x, y, z coordinates of previous and current location: x points south, y to east and z upwards*/
  double vec[3];       /* cartesian x, y, z components of the propagation vector of previous and current steps */
  /*double ds[2];*/       /* previous and current free path lengths */
  /*int cur;*/             /* either 0 or 1 to flip indices of previous and current location */
  double dir[2];       /* previous and current  directions. dir[0]: (theta) angle between ez and direction vector; dir[1]: (phi) angle between ex and projection of dircection vector onto x-y plane */
  double stokes[5];    /* Ip, Itheta, Iphi, U, V */
  double time;         /* lapse time */
  int mode;            /* 0 for P-wave, 1 for S-wave */
  double pol;          /* polarization angle: looking against the ray, the clockwise angle between vertical plane containing the ray and the plane containing the ray and the polarization direction (same as 2pi-phi from Ishimaru[1987] p.35, rotation of stokes parameters) */
}particle;


#endif
