#ifndef THERMOCHEMISTRY_H
#define THERMOCHEMISTRY_H

double xh2solve(double Tu, double xh2, double xhe1, double xhe2, double xhe3, double dti, 
double NH1, double Nhe1, double Nhe2, double R, double nhi, double *xrate, double t, double gbkg);

double xhe1solve(double Tu, double xh2i, double xhe1, double xhe2, double xhe3, double dti,
double NH1, double Nhe1, double Nhe2, double R, double nhi, double t, double gbkg_he1);

double xhe2solve(double Tu, double xh2i, double xhe1, double xhe2, double xhe3, double dti,
double NH1, double Nhe1, double Nhe2, double R, double nhi, double *hrate, double t, double gbkg_he1, double gbkg_he2);

double xhe3solve(double Tu, double xh2i, double xhe1, double xhe2, double xhe3, double dti,
double NH1, double Nhe1, double Nhe2, double R, double nhi, double *hrate, double t, double gbkg_he2);

double Tsolve(double Tu, double NH1, double Nhe1, double Nhe2,
double xh2i, double xhe1, double xhe2, double xhe3,
double nhi, double dti, double R, double tsim,
double *trate, double *xrate, double t, double gbkg, double gbkg_he1, double gbkg_he2, double heatbg);

long int xTiter(double dti, double NH1, double Nhe1, double Nhe2, double R, double nhi,
double tsim, double ti, double *Tt, double *xh2t, double *xhe1t, double *xhe2t,
double *xhe3t, double *dtrec, double gbkg, double gbkg_he1, double gbkg_he2, double heatbg);

#endif
