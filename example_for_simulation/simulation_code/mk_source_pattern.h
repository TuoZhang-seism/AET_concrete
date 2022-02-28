#if defined(MK_SOURCE_PATTERN_H)
#else
#define MK_SOURCE_PATTERN_H
int MkSourcePattern(simulationparameters, float **);
int ReadRadPat(simulationparameters, float **);
int source_rot_mat(simulationparameters, double *, double *, double *);
int SourceRot(double *, double *, double *, double *, double *);


#endif
