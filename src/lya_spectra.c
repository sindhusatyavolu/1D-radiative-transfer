#include "../include/lya_spectra.h"
#include "../include/utils.h"
#include "../include/constants.h"
#include "../include/params.h"
#include "../include/proximity_zone.h"


#include <stdio.h>
#include <stdlib.h>
#include <math.h>



/* Routine for computing the Lyman alpha optical depth and spectrum */
double absorb_along_lines_of_sight(double *Temp, double *NHI, double *Vpec, double *Tau)
{
	int ncells = boxsize * hubble_inv / (cellsize * hubble_inv);
	int extra = 0; // number of cells to be included redward of the Lyman alpha line. Set to zero if you want the first cell of the sightline to be the location of the Lyman alpha line
	int encells = ncells + extra;

	double dz, dv, b, fac, aa;
	long long bin, k;
	double profile_H1, T0, T1, T2;
	double Hubble_h = 0.678;
	double Omega = 0.308;
	double Lambda = 0.692;
	double Time = 1.0 / (1.0 + z);
	double UnitLength_in_cm = 3.085678e21;
	double UnitMass_in_g = 1.989e43;
	double UnitVelocity_in_cm_per_s = 1e5;
	double UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
	double Hubble = HUBBLE * UnitTime_in_s;

	double H_a = Hubble * sqrt(Omega / (Time * Time * Time) +
							   Lambda + (1 - Omega - Lambda) / (Time * Time));

	double Wmax = Time * H_a * boxsize / Hubble_h;

	//	double *Tau = (double*)malloc(encells*sizeof(double));
	/*
		printf("Omega=%e\n",Omega);
		printf("Lambda=%e\n",Lambda);
		printf("Time=%e\n",Time);
		printf("z=%e\n",1.0/Time-1.0);
		printf("H_a=%e\n",H_a);
		printf("boxsize=%e\n",boxsize);

		printf("T :%0.12e\n",*Temp );
		printf("V : %0.12e\n",*Vpec );
		printf ("Wmax= %g  Time=%g\n", Wmax, Time);
	*/
	dz = boxsize / ncells;
	//	printf("cellsize=%e\n",dz);

	for (bin = 0; bin < encells; bin++)
	{
		Tau[bin] = 0;

		for (k = extra; k < encells; k++)
		{
			dv = (k - bin);

			/*            while (dv < -ncells / 2)
					  dv += ncells;
					  while (dv > ncells / 2)
					  dv -= ncells;
					  */
			dv = (dv * Wmax / ncells + Vpec[k - extra]) *
				 UnitVelocity_in_cm_per_s;

			b = sqrt(2 * BOLTZMANN * Temp[k - extra] / PROTONMASS);

			aa = GAMMA_LYA_H1 * LYMAN_ALPHA / (4.0 * PI) / b;

			T0 = dv * dv / b / b;
			T1 = exp(-T0);
			T2 = 1.5 / T0;

			profile_H1 = (T0 < 1.0e-6)
							 ? T1
							 : T1 - aa / sqrt(PI) / T0 * (T1 * T1 * (4.0 * T0 * T0 + 7.0 * T0 + 4.0 + T2) - T2 - 1.0);
			//   profile_H1 = exp (-dv * dv / (b * b)) ;
			//   profile_H1 = 1.0;
			Tau[bin] += NHI[k - extra] * profile_H1 / b * dz;
		}
	}

	/* multiply with correct prefactors */

	/*  to get things into cgs units */
	fac = 1.0 / pow(UnitLength_in_cm, 2);
	fac *= Hubble_h * Hubble_h;

	/* Ly-alpha cross section */
	fac *= OSCILLATOR_STRENGTH * PI * LYMAN_ALPHA *
		   sqrt(3 * THOMPSON / (8 * PI));

	fac *= Cv / (Time * Time) / sqrt(PI);

	/* Note: For HeII, the oscillator strength is equal to that of HI,
	   and the Lyman-alpha wavelength is 4 times shorter */

	for (bin = 0; bin < encells; bin++)
	{
		Tau[bin] *= fac;
		// printf("%le\n",Tau[bin]);
	}
	/* compute proximity zone size */
	int lenY;
	double h[2048 + extra];
	for (int k = 0; k < 2048 + extra; k++)
	{
		h[k] = exp(-Tau[k]);
	}
	double xs[60];
	for (int i = 0; i < 60; ++i)
	{
		xs[i] = 1.0;
	}
	double *yt = convolve(h, xs, 2048 + extra, 60, &lenY);

	//	free(yt);
	//	printf("rp %0.12e\n",Rp(yt,2048));

	//  if(Ng ==0)
	//   return 0.0;
	//  else
	return Rp(yt, 2048 + extra);
	//	return 1;
}

double dmod(double x, double y)
{
	return x - (int)(x / y) * y;
}
