#include <stdio.h>

#ifndef SAVE_OUTPUTS_H
#define SAVE_OUTPUTS_H

char *outfilename(int los, float Mag, double t_on, double fduty, double z, double tq);
void save_outputs(FILE *fptr, int ncells, double *nh, double **nh2, double **T, double **nhe1, double **nhe2,
double **nhe3, double *Tau, int los, float Mag, double zqso, int storej);


#endif
