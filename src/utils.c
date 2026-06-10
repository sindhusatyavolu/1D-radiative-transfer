#include "../include/utils.h"
#include "../include/constants.h"
#include "../include/params.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/* Convolution to return smoothed spectra */
double *convolve(double h[], double x[], int lenH, int lenX, int *lenY)
{
	int nconv = lenH + lenX - 1; // MAX(lenH,lenX);//
	int samecon = max(lenH, lenX);
	(*lenY) = samecon; // nconv;
	int i, j, h_start, x_start, x_end;

	double *y = (double *)calloc(nconv, sizeof(double));
	double *yt = (double *)calloc(samecon, sizeof(double));

	for (i = 0; i < nconv; i++)
	{
		x_start = max(0, i - lenH + 1);
		x_end = min(i + 1, lenX);
		h_start = min(i, lenH - 1);

		for (j = x_start; j < x_end; j++)
		{
			y[i] += h[h_start--] * x[j];
		}
	}
	for (int i = 0; i < samecon; ++i)
	{
		yt[i] = y[29 + i] / 60.0;
	}

	return yt;
}

/* Cumilative sum for H column density */
double cumsum(double *nH, double **nH2, int i, int j)
{
	double sum = 0;
	for (int l = 0; l < i; l++)
	{
		sum = sum + nH[l] - nH2[l][j];
	}
	return sum;
}

/* cumulative sum for He column density */
double hecumsum(double **nHe, int i, int j)
{
	double sum = 0;
	for (int l = 0; l < i; l++)
	{
		sum = sum + nHe[l][j];
	}
	return sum;
}

/* Convert magnitude to ionizing photons per sec */
double mag_to_ndot(double mag)
{
	double L1450 = pow(10, (51.6 - mag) / 2.5);
	double LHI = L1450 / pow(912.0 / 1450.0, alpha2);
	return -LHI / (6.62e-27 * alpha);
}


/* Hubble parameter at cosmic time t */
double Hubble(double t)
{
	double omega_m = 0.308;
	double omega_l = 1 - omega_m;
	double H0 = 67.8 * 1.05e-6;
	double a = pow(sqrt(omega_m / omega_l) * sinh(1.5 * H0 * sqrt(omega_l) * t), 2.0 / 3.0); // pow(1+z,-3);
	double zz = 1.0 / a - 1.0;
	//      return H0*sqrt(omega_m/pow(a,3) +1-omega_m);
	return sqrt(omega_m) * H0 * cosh(1.5 * H0 * sqrt(1 - omega_m) * t) * pow(1 + zz, 1.5); /// sinh(1.5*H0*sqrt(1-omega_m)*t); //sqrt(omega_m/a3 +1-omega_m);
}

/* scale factor */
double ascale(double t)
{
	double omega_m = 0.308;
	double omega_l = 1 - omega_m;
	double H0 = 67.8 * 1.05e-6;
	return 1.0 / (1.0 + z); // when input units are in comoving else change a = 1.0 , a = 1/1+z   //pow(sqrt(omega_m/omega_l)*sinh(1.5*H0*sqrt(omega_l)*t), 2.0/3.0);
}

/* Mean cosmological hydrogen number density at given redshift */
double getnh(double redshift)
{
	double omegam = 0.3;
	double omegab = 0.048;
	double omegal = 0.7;
	double h = 0.68;
	double H0 = 100 * h; // km s^-1 Mpc^-1
	double G = 6.67e-11; // Nm^2/kg^2 == m^3s^-2kg^-1?
	double kmtoMPc = 1.0 / 3.08e19;
	double rhocr = pow((kmtoMPc), 2) * 3.0 * pow(H0, 2) / (8.0 * PI * G); // km^2  mpc^-2/m^3kg^-1 = kg/m^3 *(km/Mpc)^2
	double Xfac = 1.0 / X;												  // #(1.0+4.0*(1.0/X -1.0));
	//	omegah^2 = rho,0/rhocr
	//	rho,0 = rho,z*(1+z)**3
	double mH = 1.67e-27; // kg
	double rhoz = rhocr * omegab * pow((1.0 + z), 3);
	return 1e-6 * rhoz / (mH * Xfac); // in comoving cm^-3
}
