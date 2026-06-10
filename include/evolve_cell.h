#ifndef EVOLVE_CELL_H
#define EVOLVE_CELL_H


double evolve_eqns(double *nh, double **nh2, double **NH, double **nhe1, double **nhe2,
double **nhe3, double **NHe1, double **NHe2, double **T,
int ncells, double ti, double tsim, double **dt, int j, double Dt);

#endif