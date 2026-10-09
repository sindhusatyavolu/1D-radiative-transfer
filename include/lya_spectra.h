#ifndef LYA_SPECTRA_H
#define LYA_SPECTRA_H

double absorb_along_lines_of_sight(double *Temp, double *NHI, double *Vpec,
                                   double *Tau, int correct_spectral_units);
double dmod(double x, double y);

#endif
