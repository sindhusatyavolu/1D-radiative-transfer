#ifndef UTILS_H 
#define UTILS_H

double *convolve(double h[], double x[], int lenH, int lenX, int *lenY);

double cumsum(double *nH, double **nH2, int i, int j);

double hecumsum(double **nHe, int i, int j);

double mag_to_ndot(double mag);

double Hubble(double t);

double ascale(double t);

double getnh(double redshift);

#endif