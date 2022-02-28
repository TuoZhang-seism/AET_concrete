#if defined(LAYER_H)
#else
#define LAYER_H
int MkInterface(simulationparameters, float **, float **);
int ReadInterface(simulationparameters, float **, float **);
int layer_down(particle *, simulationparameters, float **, float **, int);
int layer_up(particle *, simulationparameters, float **, float **, int);
int layer(particle *, simulationparameters, float **, float **, int, int);
void incidence_(int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
#endif
