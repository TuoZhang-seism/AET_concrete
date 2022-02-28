#if defined(EUTILITY_H)
#else
#define EUTILITY_H

int ReadParFile(int, char **, simulationparameters*);
int ReadModelFiles(simulationparameters*);
gridtype ***** GridAlloc(simulationparameters);
int GridFree(simulationparameters, gridtype *****);
int checkSimulationParameters(simulationparameters *);
int WriteOutput(simulationparameters, gridtype *****,gridtype *****,bodytable *);
double gammln(double);
int stokesRot(double, int, double, double, double, double, double, double *, double *,double *, double *, double *);
int whichbody(particle *, simulationparameters *);
#endif
