#include "../include/photoion.h"
#include "../include/constants.h"
#include "../include/params.h"
#include "../include/rates.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef RTE_CACHE_FREQUENCY_GRID
static _Thread_local double *frequency_grid_cache;
static _Thread_local double *frequency_spectral_cache;
static _Thread_local int frequency_grid_nfreq;
static _Thread_local double frequency_grid_fend;
static _Thread_local double frequency_grid_alpha;

static const double *frequency_grid(void)
{
	if (frequency_grid_cache == NULL || frequency_spectral_cache == NULL ||
		frequency_grid_nfreq != NFREQ ||
		memcmp(&frequency_grid_fend, &fend, sizeof(fend)) != 0 ||
		memcmp(&frequency_grid_alpha, &alpha, sizeof(alpha)) != 0)
	{
		double *new_grid = realloc(frequency_grid_cache, NFREQ * sizeof(*new_grid));
		double *new_spectral = realloc(frequency_spectral_cache,
									   NFREQ * sizeof(*new_spectral));
		if (new_grid == NULL || new_spectral == NULL)
		{
			fprintf(stderr, "Failed to allocate photoionization frequency grid\n");
			exit(1);
		}
		frequency_grid_cache = new_grid;
		frequency_spectral_cache = new_spectral;
		frequency_grid_nfreq = NFREQ;
		frequency_grid_fend = fend;
		frequency_grid_alpha = alpha;
		double a = 13.6;
		double b = 13.6 * fend;
		for (int i = 0; i < NFREQ; i++)
		{
			frequency_grid_cache[i] = a * pow(b / a, i / (NFREQ - 1.0));
			frequency_spectral_cache[i] = pow(frequency_grid_cache[i] / a, alpha - 1);
		}
	}
	return frequency_grid_cache;
}
#endif

#ifdef RTE_CACHE_FREQUENCY_GRID
#define DECLARE_FREQUENCY_GRID()                  \
	const double *u = frequency_grid();            \
	const double *spectral_weight = frequency_spectral_cache
#define SPECTRAL_WEIGHT(index, reference) spectral_weight[(index)]
#else
#define DECLARE_FREQUENCY_GRID()                     \
	double a = 13.6;                                   \
	double b = 13.6 * fend;                            \
	double u[NFREQ];                                   \
	for (int i = 0; i < NFREQ; i++)                   \
		u[i] = a * pow(b / a, i / (NFREQ - 1.0))
#define SPECTRAL_WEIGHT(index, reference) pow(u[(index)] / (reference), alpha - 1)
#endif

/* Photoionisation rates in energy units eV **/

double Gamma(double NH, double Nhe1, double Nhe2, double nh1, double nhe1l, double nhe2l, double t)
{
	DECLARE_FREQUENCY_GRID();

	double fsum = 0;
	double x0 = u[0];
	int len = NFREQ;

	/*Integral for calculating photoionisation rate*/
	for (int i = 1; i < len; i++)
	{

		double df = (u[i] - u[i - 1]) * 0.5;
		fsum = fsum + df * (SPECTRAL_WEIGHT(i, x0) * exp(-sigma(u[i]) * NH - sigmahe1(u[i]) * Nhe1 - sigmahe2(u[i]) * Nhe2) * PH1(u[i], nh1, nhe1l, nhe2l, t) + SPECTRAL_WEIGHT(i - 1, x0) * exp(-sigma(u[i - 1]) * NH - sigmahe1(u[i - 1]) * Nhe1 - sigmahe2(u[i - 1]) * Nhe2) * PH1(u[i - 1], nh1, nhe1l, nhe2l, t));
	}

	return fsum;
}

double Gammahe1(double NH, double Nhe1, double Nhe2, double nh1, double nhe1l, double nhe2l, double t)
{
	DECLARE_FREQUENCY_GRID();

	double fsum = 0;
	double x0 = 13.6; // u[0];
	int len = NFREQ;

	for (int i = 1; i < len; i++)
	{
		double df = (u[i] - u[i - 1]) * 0.5;
		fsum = fsum + df * (SPECTRAL_WEIGHT(i, x0) * exp(-sigma(u[i]) * NH - sigmahe1(u[i]) * Nhe1 - sigmahe2(u[i]) * Nhe2) * PHe1(u[i], nh1, nhe1l, nhe2l, t) + SPECTRAL_WEIGHT(i - 1, x0) * exp(-sigma(u[i - 1]) * NH - sigmahe1(u[i - 1]) * Nhe1 - sigmahe2(u[i - 1]) * Nhe2) * PHe1(u[i - 1], nh1, nhe1l, nhe2l, t));
	}

	return fsum;
}

double Gammahe2(double NH, double Nhe1, double Nhe2, double nh1, double nhe1l, double nhe2l, double t)
{
	DECLARE_FREQUENCY_GRID();

	double fsum = 0;
	double x0 = 13.6; // u[0];
	int len = NFREQ;

	for (int i = 1; i < len; i++)
	{
		double df = (u[i] - u[i - 1]) * 0.5;
		fsum = fsum + df * (SPECTRAL_WEIGHT(i, x0) * exp(-sigma(u[i]) * NH - sigmahe1(u[i]) * Nhe1 - sigmahe2(u[i]) * Nhe2) * PHe2(u[i], nh1, nhe1l, nhe2l, t) + SPECTRAL_WEIGHT(i - 1, x0) * exp(-sigma(u[i - 1]) * NH - sigmahe1(u[i - 1]) * Nhe1 - sigmahe2(u[i - 1]) * Nhe2) * PHe2(u[i - 1], nh1, nhe1l, nhe2l, t));
	}

	return fsum;
}

/* Photoheating rates - ev^2 */

double Heat(double NH, double Nhe1, double Nhe2, double nh1, double nhe1l, double nhe2l, double t)
{
	DECLARE_FREQUENCY_GRID();

	double fsum = 0;
	double x0 = u[0];
	int len = NFREQ;

	for (int i = 1; i < len; i++)
	{
		double df = (u[i] - u[i - 1]) * 0.5;
		fsum = fsum + df * ((u[i] - x0) * SPECTRAL_WEIGHT(i, x0) * exp(-sigma(u[i]) * NH - sigmahe1(u[i]) * Nhe1 - sigmahe2(u[i]) * Nhe2) * PH1(u[i], nh1, nhe1l, nhe2l, t) + (u[i - 1] - x0) * SPECTRAL_WEIGHT(i - 1, x0) * exp(-sigma(u[i - 1]) * NH - sigmahe1(u[i - 1]) * Nhe1 - sigmahe2(u[i - 1]) * Nhe2) * PH1(u[i - 1], nh1, nhe1l, nhe2l, t));
		// printf("%0.12e\n",PH1(u[i],nh1,nhe1l,nhe2l));
	}

	return fsum;
}

double Heat2(double NH, double Nhe1, double Nhe2, double nh1, double nhe1l, double nhe2l, double t)
{
	DECLARE_FREQUENCY_GRID();

	double fsum = 0;
	double x0 = u[0];
	int len = NFREQ;

	for (int i = 1; i < len; i++)
	{
		double df = (u[i] - u[i - 1]) * 0.5;
		fsum = fsum + df * ((u[i] - x0) * SPECTRAL_WEIGHT(i, 13.6) * exp(-sigma(u[i]) * NH - sigmahe1(u[i]) * Nhe1 - sigmahe2(u[i]) * Nhe2) * PHe1(u[i], nh1, nhe1l, nhe2l, t) + (u[i - 1] - x0) * SPECTRAL_WEIGHT(i - 1, 13.6) * exp(-sigma(u[i - 1]) * NH - sigmahe1(u[i - 1]) * Nhe1 - sigmahe2(u[i - 1]) * Nhe2) * PHe1(u[i - 1], nh1, nhe1l, nhe2l, t));
	}

	return fsum;
}

double Heat3(double NH, double Nhe1, double Nhe2, double nh1, double nhe1l, double nhe2l, double t)
{
	DECLARE_FREQUENCY_GRID();

	double fsum = 0;
	double x0 = u[0];
	int len = NFREQ;

	for (int i = 1; i < len; i++)
	{
		double df = (u[i] - u[i - 1]) * 0.5;
		fsum = fsum + df * ((u[i] - x0) * SPECTRAL_WEIGHT(i, 13.6) * exp(-sigma(u[i]) * NH - sigmahe1(u[i]) * Nhe1 - sigmahe2(u[i]) * Nhe2) * PHe2(u[i], nh1, nhe1l, nhe2l, t) + (u[i - 1] - x0) * SPECTRAL_WEIGHT(i - 1, 13.6) * exp(-sigma(u[i - 1]) * NH - sigmahe1(u[i - 1]) * Nhe1 - sigmahe2(u[i - 1]) * Nhe2) * PHe2(u[i - 1], nh1, nhe1l, nhe2l, t));
	}

	return fsum;
}
