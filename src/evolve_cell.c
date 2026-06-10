#include "../include/evolve_cell.h"
#include "../include/utils.h"
#include "../include/rates.h"
#include "../include/photoion.h"
#include "../include/thermochemistry.h"
#include "../include/params.h"
#include "../include/constants.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>




double evolve_eqns(double *nh, double **nh2, double **NH, double **nhe1, double **nhe2,
				   double **nhe3, double **NHe1, double **NHe2, double **T,
				   int ncells, double ti, double tsim, double **dt, int j, double Dt)
{
	//      double dt = Dt;
	//	double X =1 ;
	//	double Y =0;
	//	double xhe2 =0;
	//	double xhe3 =0;
	long int out;
	double *Tt = (double *)malloc(max_subcycling_steps * sizeof(double));
	double *xh2t = (double *)malloc(max_subcycling_steps * sizeof(double));
	double *xhe1t = (double *)malloc(max_subcycling_steps * sizeof(double));
	double *xhe2t = (double *)malloc(max_subcycling_steps * sizeof(double));
	double *xhe3t = (double *)malloc(max_subcycling_steps * sizeof(double));
	double dt_rec;
	double *Rcell = (double *)malloc(ncells * sizeof(double));
	for (int k = 0; k < ncells; k++)
	{
		//		dt[k][j+1] = dt[k][j];
		//		printf("cell %d\n",k);
		double nu;
		xh2t[0] = nh2[k][j] / nh[k];
		xhe2t[0] = nhe2[k][j] / ((Y / (4.0 * ((1.0 - Y)))) * nh[k]);
		xhe3t[0] = nhe3[k][j] / ((Y / (4.0 * ((1.0 - Y)))) * nh[k]);
		xhe1t[0] = 1.0 - xhe2t[0] - xhe3t[0];
		nu = (1.0 + 4.0 * Y) / (1.0 + Y + xh2t[0] + xhe2t[0] + 2 * xhe3t[0]);
		Tt[0] = (3.0 / 2.0) * T[k][j] * kb / (nu * PROTONMASS);
		Rcell[k] = (Roffset + k * cellsize) * hubble_inv * ascale(tsim);
		out = xTiter(dt[k][j], NH[k][j], NHe1[k][j], NHe2[k][j], Rcell[k], nh[k], tsim, ti, Tt, xh2t, xhe1t, xhe2t, xhe3t, &dt_rec, Gammabkg[k], Gammabkg_he1[k], Gammabkg_he2[k], Heat_bkg[k]);
		// equivalent of cool_step routine of Ramses RT (Rosdahl et al 2013)
		while ((out == 100000000000 || out == 100000000001 || out == 100000000002 || out == 100000000003))
		{
			dt[k][j] = (double)max(dt_min, dt[k][j] / 2);
			//      printf("%ld %e\n",out-100000000000, dt[k][j]);
			out = xTiter(dt[k][j], NH[k][j], NHe1[k][j], NHe2[k][j], Rcell[k], nh[k], tsim, ti, Tt, xh2t, xhe1t, xhe2t, xhe3t, &dt_rec, Gammabkg[k], Gammabkg_he1[k], Gammabkg_he2[k], Heat_bkg[k]);
		}

		nh2[k][j + 1] = pow(ascale(ti), 3) * nh[k] * xh2t[out] / pow(ascale(tsim), 3);
		nhe2[k][j + 1] = (Y / (4.0 * ((1.0 - Y)))) * pow(ascale(ti), 3) * nh[k] * xhe2t[out] / pow(ascale(tsim), 3);
		nhe3[k][j + 1] = (Y / (4.0 * ((1.0 - Y)))) * pow(ascale(ti), 3) * nh[k] * xhe3t[out] / pow(ascale(tsim), 3);
		nhe1[k][j + 1] = (Y / (4.0 * ((1.0 - Y)))) * pow(ascale(ti), 3) * nh[k] * xhe1t[out] / pow(ascale(tsim), 3); //(Y/(4*((1.0-Y))))*nh[k] - nhe2[k][j+1] - nhe3[k][j+1];
		nu = (1.0 + 4.0 * Y) / (1.0 + Y + xh2t[out] + xhe2t[out] + 2 * xhe3t[out]);
		T[k][j + 1] = (2.0 / 3.0) * Tt[out] * (nu * PROTONMASS) / (kb);
		//		dt[k][j+1] = Dt;//(double) min((dt_rec),tsim-ti); // check this criterion
	}
	for (int k = 0; k < ncells; k++)
	{
		NH[k][j + 1] = cumsum(nh, nh2, k, j + 1) * cellsize * hubble_inv * kpc * nm2 * ascale(tsim); // in cm^-2
		NHe1[k][j + 1] = hecumsum(nhe1, k, j + 1) * cellsize * hubble_inv * kpc * nm2 * ascale(tsim);
		NHe2[k][j + 1] = hecumsum(nhe2, k, j + 1) * cellsize * hubble_inv * kpc * nm2 * ascale(tsim);
		nh[k] = pow(ascale(ti), 3) * nh[k] / pow(ascale(tsim), 3);
	}
	free(Tt);
	free(xh2t);
	free(xhe2t);
	free(xhe3t);
	free(xhe1t);
	return dt_rec;
}
