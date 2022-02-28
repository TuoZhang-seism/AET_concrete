#if defined(SURFACE_H)
#else
#define SURFACE_H
int MkSurface(simulationparameters, float **, float **);
int ReadSurface(simulationparameters, float **, float **);
int surface(particle *, int, bodyparm, float **, float **, double *);
#endif
