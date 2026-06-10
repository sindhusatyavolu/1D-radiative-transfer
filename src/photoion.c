#include "../include/photoion.h"
#include "../include/constants.h"
#include "../include/params.h"
#include "../include/rates.h"
#include <stdio.h>
#include <math.h>


/* Photoionisation rates in energy units eV **/

double Gamma(double NH, double Nhe1, double Nhe2, double nh1, double nhe1l, double nhe2l, double t)
{
	double a = 13.6;
	double b = 13.6 * fend;
	double u[NFREQ];

	for (int i = 0; i < NFREQ; i++)
	{
		u[i] = a * pow(b / a, i / (NFREQ - 1.0));
		// printf("%0.12e\n",u[i]);
	}

	double fsum = 0;
	double x0 = u[0];
	int len = sizeof(u) / sizeof(u[0]);

	/*Integral for calculating photoionisation rate*/
	for (int i = 1; i < len; i++)
	{

		double df = (u[i] - u[i - 1]) * 0.5;
		fsum = fsum + df * (pow(u[i] / x0, alpha - 1) * exp(-sigma(u[i]) * NH - sigmahe1(u[i]) * Nhe1 - sigmahe2(u[i]) * Nhe2) * PH1(u[i], nh1, nhe1l, nhe2l, t) + pow(u[i - 1] / x0, alpha - 1) * exp(-sigma(u[i - 1]) * NH - sigmahe1(u[i - 1]) * Nhe1 - sigmahe2(u[i - 1]) * Nhe2) * PH1(u[i - 1], nh1, nhe1l, nhe2l, t));
	}

	return fsum;
}

double Gammahe1(double NH, double Nhe1, double Nhe2, double nh1, double nhe1l, double nhe2l, double t)
{
	double a = 13.6;
	double b = 13.6 * fend;
	double u[NFREQ];

	for (int i = 0; i < NFREQ; i++)
	{
		u[i] = a * pow(b / a, i / (NFREQ - 1.0));
	}

	double fsum = 0;
	double x0 = 13.6; // u[0];
	int len = sizeof(u) / sizeof(u[0]);

	for (int i = 1; i < len; i++)
	{
		double df = (u[i] - u[i - 1]) * 0.5;
		fsum = fsum + df * (pow(u[i] / x0, alpha - 1) * exp(-sigma(u[i]) * NH - sigmahe1(u[i]) * Nhe1 - sigmahe2(u[i]) * Nhe2) * PHe1(u[i], nh1, nhe1l, nhe2l, t) + pow(u[i - 1] / x0, alpha - 1) * exp(-sigma(u[i - 1]) * NH - sigmahe1(u[i - 1]) * Nhe1 - sigmahe2(u[i - 1]) * Nhe2) * PHe1(u[i - 1], nh1, nhe1l, nhe2l, t));
	}

	return fsum;
}

double Gammahe2(double NH, double Nhe1, double Nhe2, double nh1, double nhe1l, double nhe2l, double t)
{
	double a = 13.6;
	double b = 13.6 * fend;
	double u[NFREQ];

	for (int i = 0; i < NFREQ; i++)
	{
		u[i] = a * pow(b / a, i / (NFREQ - 1.0));
	}

	double fsum = 0;
	double x0 = 13.6; // u[0];
	int len = sizeof(u) / sizeof(u[0]);

	for (int i = 1; i < len; i++)
	{
		double df = (u[i] - u[i - 1]) * 0.5;
		fsum = fsum + df * (pow(u[i] / x0, alpha - 1) * exp(-sigma(u[i]) * NH - sigmahe1(u[i]) * Nhe1 - sigmahe2(u[i]) * Nhe2) * PHe2(u[i], nh1, nhe1l, nhe2l, t) + pow(u[i - 1] / x0, alpha - 1) * exp(-sigma(u[i - 1]) * NH - sigmahe1(u[i - 1]) * Nhe1 - sigmahe2(u[i - 1]) * Nhe2) * PHe2(u[i - 1], nh1, nhe1l, nhe2l, t));
	}

	return fsum;
}

/* Photoheating rates - ev^2 */

double Heat(double NH, double Nhe1, double Nhe2, double nh1, double nhe1l, double nhe2l, double t)
{
	double a = 13.6;
	double b = 13.6 * fend;
	double u[NFREQ];

	for (int i = 0; i < NFREQ; i++)
	{
		u[i] = a * pow(b / a, i / (NFREQ - 1.0));
	}

	double fsum = 0;
	double x0 = u[0];
	int len = sizeof(u) / sizeof(u[0]);

	for (int i = 1; i < len; i++)
	{
		double df = (u[i] - u[i - 1]) * 0.5;
		fsum = fsum + df * ((u[i] - x0) * pow(u[i] / x0, alpha - 1) * exp(-sigma(u[i]) * NH - sigmahe1(u[i]) * Nhe1 - sigmahe2(u[i]) * Nhe2) * PH1(u[i], nh1, nhe1l, nhe2l, t) + (u[i - 1] - x0) * pow(u[i - 1] / x0, alpha - 1) * exp(-sigma(u[i - 1]) * NH - sigmahe1(u[i - 1]) * Nhe1 - sigmahe2(u[i - 1]) * Nhe2) * PH1(u[i - 1], nh1, nhe1l, nhe2l, t));
		// printf("%0.12e\n",PH1(u[i],nh1,nhe1l,nhe2l));
	}

	return fsum;
}

double Heat2(double NH, double Nhe1, double Nhe2, double nh1, double nhe1l, double nhe2l, double t)
{
	double a = 13.6;
	double b = 13.6 * fend;
	double u[NFREQ];

	for (int i = 0; i < NFREQ; i++)
	{
		u[i] = a * pow(b / a, i / (NFREQ - 1.0));
	}

	double fsum = 0;
	double x0 = u[0];
	int len = sizeof(u) / sizeof(u[0]);

	for (int i = 1; i < len; i++)
	{
		double df = (u[i] - u[i - 1]) * 0.5;
		fsum = fsum + df * ((u[i] - x0) * pow(u[i] / 13.6, alpha - 1) * exp(-sigma(u[i]) * NH - sigmahe1(u[i]) * Nhe1 - sigmahe2(u[i]) * Nhe2) * PHe1(u[i], nh1, nhe1l, nhe2l, t) + (u[i - 1] - x0) * pow(u[i - 1] / 13.6, alpha - 1) * exp(-sigma(u[i - 1]) * NH - sigmahe1(u[i - 1]) * Nhe1 - sigmahe2(u[i - 1]) * Nhe2) * PHe1(u[i - 1], nh1, nhe1l, nhe2l, t));
	}

	return fsum;
}

double Heat3(double NH, double Nhe1, double Nhe2, double nh1, double nhe1l, double nhe2l, double t)
{
	double a = 13.6;
	double b = 13.6 * fend;
	double u[NFREQ];

	for (int i = 0; i < NFREQ; i++)
	{
		u[i] = a * pow(b / a, i / (NFREQ - 1.0));
	}

	double fsum = 0;
	double x0 = u[0];
	int len = sizeof(u) / sizeof(u[0]);

	for (int i = 1; i < len; i++)
	{
		double df = (u[i] - u[i - 1]) * 0.5;
		fsum = fsum + df * ((u[i] - x0) * pow(u[i] / 13.6, alpha - 1) * exp(-sigma(u[i]) * NH - sigmahe1(u[i]) * Nhe1 - sigmahe2(u[i]) * Nhe2) * PHe2(u[i], nh1, nhe1l, nhe2l, t) + (u[i - 1] - x0) * pow(u[i - 1] / 13.6, alpha - 1) * exp(-sigma(u[i - 1]) * NH - sigmahe1(u[i - 1]) * Nhe1 - sigmahe2(u[i - 1]) * Nhe2) * PHe2(u[i - 1], nh1, nhe1l, nhe2l, t));
	}

	return fsum;
}
