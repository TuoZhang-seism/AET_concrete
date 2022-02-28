#if defined(SCAT_COEFF_H)
#else
#define SCAT_COEFF_H



/* Scattering Pattern */
#define XPPP(GAM, NY, THE)(pow(GAM,-2)*(NY*(cos(THE)+2/pow(GAM,2)*pow(sin(THE),2)-1)-2+4/pow(GAM,2)*pow(sin(THE),2)))
#define XPSL(GAM, NY, THE)(-sin(THE)*(NY*(1-2/GAM*cos(THE))-4/GAM*cos(THE)))
#define XSPP(GAM, NY, THE, PHI)(pow(GAM,-2)*sin(THE)*cos(PHI)*(NY*(1-2/GAM*cos(THE))-4/GAM*cos(THE)))
#define XSSL(GAM, NY, THE, PHI)(cos(PHI)*(NY*(cos(THE)-cos(2*THE))-2*cos(2*THE)))
#define XSSR(GAM, NY, THE, PHI)(sin(PHI)*(NY*(cos(THE)-1)+2*cos(THE)))

/* Proportional to THETA part of scattering pattern */
#define XPPP_THETA(GAM, NY, THE)(NY*(cos(THE)+2/pow(GAM,2)*pow(sin(THE),2)-1)-2+4/pow(GAM,2)*pow(sin(THE),2))
#define XPSL_THETA(GAM, NY, THE)(sin(THE)*(NY*(1-2/GAM*cos(THE))-4/GAM*cos(THE)))
#define XSPP_THETA(GAM, NY, THE)(sin(THE)*(NY*(1-2/GAM*cos(THE))-4/GAM*cos(THE)))
#define XSSL_THETA(GAM, NY, THE)(NY*(cos(THE)-cos(2*THE))-2*cos(2*THE))
#define XSSR_THETA(GAM, NY, THE)(NY*(cos(THE)-1)+2*cos(THE))

/* Proportional to PHI part of scattering pattern */
#define XPPP_PHI 1
#define XPSL_PHI 1
#define XSPP_PHI(PHI)(cos(PHI))
#define XSSL_PHI(PHI)(cos(PHI))
#define XSSR_PHI(PHI)(sin(PHI))

/* Proprotionality factor (XSSL(GAM, NY, THE, PHI) = XSSL_THETA(GAM, NY, THE) * XSSL_PHI(PHI) * XSSL_FAC) */
#define XPPP_FAC(GAM)(pow(GAM,-2))
#define XPSL_FAC -1.
#define XSPP_FAC(GAM)(pow(GAM,-2))
#define XSSL_FAC 1.
#define XSSR_FAC 1.


/* Power spectral density functions */
/* 3D
#define PSDF_GAUSS(M,E,A)(pow(E,2)*sqrt(pow(PI,3))*pow(A,3)*exp(-pow(M,2)*pow(A,2)/4))
#define PSDF_EXP(M,E,A)(8*PI*pow(E,2)*pow(A,3)/pow(1+pow(A,2)*pow(M,2),2))
#define PSDF_KARMAN(M,E,A,KAP)(8*pow(PI,1.5)*E*E*A*A*A/pow(1+A*A*M*M,KAP+1.5) * exp(gammln(KAP+1.5) - gammln(KAP)))*/
/* 2D  */
#define PSDF_GAUSS(M,E,A)(pow(E,2)*PI*pow(A,2)*exp(-pow(M,2)*pow(A,2)/4))
#define PSDF_EXP(M,E,A)(2*PI*pow(E,2)*pow(A,2)/pow(1+pow(A,2)*pow(M,2),1.5))
#define PSDF_KARMAN(M,E,A,KAP)(4*PI*E*E*A*A/pow(1+A*A*M*M,KAP+1))

/* Proportional to angle dependent part of power spectral density function */
/* 3D
#define PSDF_GAUSS_ANG(M,A)(exp(-pow(M,2)*pow(A,2)/4))
#define PSDF_EXP_ANG(M,A)(1./pow(1+pow(A,2)*pow(M,2),2))
#define PSDF_KARMAN_ANG(M,A,KAP)(1/pow(1+A*A*M*M,KAP+1.5))*/
/* 2D */
#define PSDF_GAUSS_ANG(M,A)(exp(-pow(M,2)*pow(A,2)/4))
#define PSDF_EXP_ANG(M,A)(1./pow(1+pow(A,2)*pow(M,2),1.5))
#define PSDF_KARMAN_ANG(M,A,KAP)(1/pow(1+A*A*M*M,KAP+1))

/* Proprotionality factor (PSDF_GAUSS = PSDF_GAUSS_ANG * PSDF_GAUSS_FAC) */
/*3D
#define PSDF_GAUSS_FAC(E,A)(pow(E,2)*sqrt(pow(PI,3))*pow(A,3))
#define PSDF_EXP_FAC(E,A)(8*PI*pow(E,2)*pow(A,3))
#define PSDF_KARMAN_FAC(E,A,KAP)(8*pow(PI,1.5)*E*E*A*A*A * exp(gammln(KAP+1.5) - gammln(KAP)))*/
#define PSDF_GAUSS_FAC(E,A)(pow(E,2)*sqrt(pow(PI,3))*pow(A,3))
#define PSDF_EXP_FAC(E,A)(2*PI*pow(E,2)*pow(A,2))
#define PSDF_KARMAN_FAC(E,A,KAP)(4*PI*E*E*A*A)

/* Source radiation pattern (point shear dislocation) */
#define BPR_THETA(THE)(pow(sin(THE),4))                           // THETA part of P-wave radiation pattern
#define BPR_PHI(PHI)(pow(sin(2.*PHI),2))                          // PHI part of P-wave radiation pattern
#define BPR(THE, PHI)(3.75*pow(sin(THE),4)*pow(sin(2.*PHI),2))    // Complete P-wave radiation pattern normalized to 4PI

#define BST_THETA(THE)(pow(sin(2.*THE),2))                        // THETA part of S_theta-wave radiation pattern
#define BST_PHI(PHI)(pow(sin(2.*PHI),2))                          // PHI part of S_theta-wave radiation pattern
#define BST(THE, PHI)(0.625*pow(sin(2.*THE)*sin(2.*PHI),2))       // Complete S_theta-wave radiation pattern normalized to 2/3PI
  
#define BSF_THETA(THE)(pow(sin(THE),2))                           // THETA part of S_phi-wave radiation pattern
#define BSF_PHI(PHI)(pow(cos(2.*PHI),2))                          // PHI part of S_phi-wave radiation pattern
#define BSF(THE,PHI)(2.5*pow(sin(THE)*cos(2.*PHI),2))             // Complete S_theta-wave radiation pattern normalized to 10/3PI


/* Source radiation pattern (isotropic) */
// #define BPR_THETA(THE)(1)   /* THETA part of P-wave radiation pattern */
// #define BPR_PHI(PHI)(1)     /* PHI part of P-wave radiation pattern */
// #define BPR(THE, PHI)(1)    /* Complete P-wave radiation pattern normalized to 4PI */

// #define BST_THETA(THE)(0.5)   /* THETA part of S_theta-wave radiation pattern */
// #define BST_PHI(PHI)(1)     /* PHI part of S_theta-wave radiation pattern */
// #define BST(THE, PHI)(1)    /*Complete S_theta-wave radiation pattern normalized to 2/3PI */

// #define BSF_THETA(THE)(0.5)   /* THETA part of S_phi-wave radiation pattern */
// #define BSF_PHI(PHI)(1)     /* PHI part of S_phi-wave radiation pattern */
// #define BSF(THE,PHI)(1)     /* Complete S_theta-wave radiation pattern normalized to 10/3PI */


#endif
