
#include "../include/rates.h"
#include "../include/utils.h"
#include "../include/params.h"
#include "../include/constants.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* electron collisonal ionisation rates in Myr^-1 cm^3 */
double GE(double Ti)
{
	double Th1 = 157807;
	double lam = 2 * Th1 / Ti;
	double T5 = Ti / 100000;
	//	return timefac*5.85e-11*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-157809.1/Ti);
	//	return timefac*1.17e-10*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-157809.1/Ti);//bolton
	if (Ti >= 1e4)
	{
		return timefac * 21.11 * pow(Ti, -1.5) * exp(-lam / 2) * pow(lam, -1.089) / pow(1.0 + pow((lam / 0.354), 0.874), 1.101); // HnG 97
	}
	else
		return 0.0;
}


double betahe1(double Ti)
{
	double lam = 2 * 285335 / Ti;
	double T5 = Ti / 100000;
	//	return timefac*2.38e-11*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-285335.4/Ti); //ramses
	//	return timefac*4.76e-11*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-285335.4/Ti);//bolton
	if (Ti >= 1e4)
	{
		return timefac * 32.38 * pow(Ti, -1.5) * exp(-lam / 2) * pow(lam, -1.146) / pow(1.0 + pow((lam / 0.416), 0.987), 1.056); // from HnG
	}
	else
		return 0.0;
}

double betahe2(double Ti)
{
	double lam = 2 * 631515 / Ti;
	//	double T5 =Ti/100000 ;
	//	return timefac*5.68e-12*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-631515.0/Ti); //ramses
	//	return timefac*1.14e-11*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-631515.0/Ti); //bolton
	if (Ti >= 1e4)
	{
		return timefac * 19.95 * pow(Ti, -1.5) * exp(-lam / 2) * pow(lam, -1.089) / pow(1.0 + pow((lam / 0.553), 0.735), 1.275); // HnG 97
	}
	else
		return 0.0;
}

/* Case A Recombination rates in cm^3 Myr^-1 */

double ralpha(double Ti)
{
	double Th1 = 157807;
	double lam = 2 * Th1 / Ti;
	//	double Tl = kb*Ti/1.6e-12;
	double ar = 7.982e-11 * timefac;
	double b = 0.7480;
	double T0 = 3.148;
	double T1 = 7.036e5;
	//	return ar*pow((sqrt(Ti/T0)*(pow(1+sqrt(Ti/T0),1.0-b))*(pow(1+sqrt(Ti/T1),1.0+b))),-1);
	//	return timefac*(exp(-28.6130338-0.72411256*log(Tl)-2.02604473e-2*log(pow(Tl,2))-2.38086188e-3*log(pow(Tl,3))-3.21260521e-4*log(pow(Tl,4))-1.42150291e-5*log(pow(Tl,5))+4.98910892e-6*log(pow(Tl,6))+5.75561414e-7*log(pow(Tl,7))-1.85676704e-8*log(pow(Tl,8))-3.07113524e-9*log(pow(Tl,9))));
	return timefac * 1.269e-13 * pow(lam, 1.503) / pow((1.0 + pow(lam / 0.522, 0.47)), 1.923); //  HnG 97
}





double alphahe2(double Ti)
{
	double lam = 570670 / Ti;
	double Tl = kb * Ti / 1.6e-12;
	double ar = 3.294e-11 * timefac;
	double b = 0.6910;
	double T0 = 15.54;
	double T1 = 3.676e7;
	if (5e3 < Ti && Ti < 5e5)
	{
		return timefac * 3e-14 * pow(lam, 0.654); // HnG 97
	}
	else
		return 0.0;
	//	return timefac*(3.925e-13*pow(Tl,-0.6353)+1.544e-9*pow(Tl,-1.5)*exp(-48.596/Tl)*(0.3+exp(8.1/Tl)));
	//	return ar*pow((sqrt(Ti/T0)*(pow(1+sqrt(Ti/T0),1.0-b))*(pow(1+sqrt(Ti/T1),1.0+b))),-1);
}

double alphahe3(double Ti)
{
	double lam = 1263030 / Ti;
	double ar = 1.891e-10 * timefac;
	double b = 0.7524;
	double T0 = 9.37;
	double T1 = 2.774e6;
	return timefac * 2.538e-13 * pow(lam, 1.503) / pow((1.0 + pow(lam / 0.522, 0.47)), 1.923); // HnG 97
	//	return 2*ralpha(Ti/4);//bolton
	//	return ar*pow((sqrt(Ti/T0)*(pow(1+sqrt(Ti/T0),1.0-b))*(pow(1+sqrt(Ti/T1),1.0+b))),-1);
}


/* Photoionisation cross sections in nm^2 */

double sigma(double f)
{
	double sigma_0 = 5.475; // in nm^2
	double nu_0 = 0.4298;
	double y0 = 0;
	double x = f / nu_0 - y0;
	double yw = 0;
	double P = 2.963;
	double ya = 32.88;
	double y1 = 0;
	double y = sqrt(x * x + y1 * y1);
	if (f < 13.6)
		return 0.0;
	else
		return sigma_0 * (pow((x - 1), 2) + pow(yw, 2)) * pow(y, 0.5 * P - 5.5) / pow(1 + sqrt(y / ya), P);
	// return 6.3e-18*1e14*(1.34*pow(f/13.6,-2.99)-0.34*pow(f/13.6,-3.99)); //bolton cross sections
}

double sigmahe1(double f)
{
	double sigma_0 = 9.492e-2; // in nm^2
	double nu_0 = 13.61;
	double y0 = 0.4434;
	double x = f / nu_0 - y0;
	double yw = 2.039;
	double P = 3.188;
	double ya = 1.469;
	double y1 = 2.136;
	double y = sqrt(x * x + y1 * y1);
	if (f < 24.6)
		return 0.0;
	else
		return sigma_0 * (pow((x - 1), 2) + pow(yw, 2)) * pow(y, 0.5 * P - 5.5) / pow(1 + sqrt(y / ya), P);
	// return 7.03e-18*1e14*(1.66*pow(f/24.6,-2.05)-0.66*pow(f/24.6,-3.05)); //bolton cross sections
}

double sigmahe2(double f)
{
	double sigma_0 = 1.369; // in nm^2
	double nu_0 = 1.720;
	double y0 = 0;
	double x = f / nu_0 - y0;
	double yw = 0;
	double P = 2.963;
	double ya = 32.88;
	double y1 = 0;
	double y = sqrt(x * x + y1 * y1);
	if (f < 54.5)
		return 0.0;
	else
		return sigma_0 * (pow((x - 1), 2) + pow(yw, 2)) * pow(y, 0.5 * P - 5.5) / pow(1 + sqrt(y / ya), P);
	//  return 1.5e-18*1e14*(1.34*pow(f/54.5,-2.99)-0.34*pow(f/54.5,-3.99)); //bolton
}

/* Probabilities */

double PH1(double x, double nh1, double nhe1l, double nhe2l, double t)
{
	double qH = exp(-sigma(x) * nh1 * cellsize * hubble_inv * kpc * nm2 * ascale(t));
	double qHe1 = exp(-sigmahe1(x) * nhe1l * cellsize * hubble_inv * kpc * nm2 * ascale(t));
	double qHe2 = exp(-sigmahe2(x) * nhe2l * cellsize * hubble_inv * kpc * nm2 * ascale(t));
	double tau_tot = sigma(x) * nh1 * cellsize * hubble_inv * kpc * ascale(t) * nm2 + sigmahe1(x) * nhe1l * cellsize * hubble_inv * kpc * nm2 * ascale(t) + sigmahe2(x) * nhe2l * cellsize * hubble_inv * kpc * nm2 * ascale(t);
	double D = max(x_min_fix, (1 - qH) * qHe2 * qHe1 + (1 - qHe1) * qH * qHe2 + (1 - qHe2) * qH * qHe1);
	// printf("%0.12e\n",D);
	if (D == 0 || D != D)
	{
		printf("D %0.12e\n", D);
		printf("%0.12e %0.12e %0.12e\n", nh1, nhe1l, nhe2l);
		exit(0);
	}
	return (1 - qH) * qHe2 * qHe1 * (1 - exp(-tau_tot)) / D;
}

double PHe1(double x, double nh1, double nhe1l, double nhe2l, double t)
{
	double qH = exp(-sigma(x) * nh1 * cellsize * hubble_inv * kpc * nm2 * ascale(t));
	double qHe1 = exp(-sigmahe1(x) * nhe1l * cellsize * hubble_inv * kpc * nm2 * ascale(t));
	double qHe2 = exp(-sigmahe2(x) * nhe2l * cellsize * hubble_inv * kpc * nm2 * ascale(t));
	double tau_tot = sigma(x) * nh1 * cellsize * hubble_inv * kpc * nm2 * ascale(t) + sigmahe1(x) * nhe1l * cellsize * hubble_inv * kpc * nm2 * ascale(t) + sigmahe2(x) * nhe2l * cellsize * hubble_inv * kpc * nm2 * ascale(t);
	double D = max(x_min_fix, (1 - qH) * qHe2 * qHe1 + (1 - qHe1) * qH * qHe2 + (1 - qHe2) * qH * qHe1);
	return (qH)*qHe2 * (1 - qHe1) * (1 - exp(-tau_tot)) / D;
}

double PHe2(double x, double nh1, double nhe1l, double nhe2l, double t)
{
	double qH = exp(-sigma(x) * nh1 * cellsize * hubble_inv * kpc * nm2 * ascale(t));
	double qHe1 = exp(-sigmahe1(x) * nhe1l * cellsize * hubble_inv * kpc * nm2 * ascale(t));
	double qHe2 = exp(-sigmahe2(x) * nhe2l * cellsize * hubble_inv * kpc * nm2 * ascale(t));
	double tau_tot = sigma(x) * nh1 * cellsize * hubble_inv * kpc * nm2 * ascale(t) + sigmahe1(x) * nhe1l * cellsize * hubble_inv * kpc * nm2 * ascale(t) + sigmahe2(x) * nhe2l * cellsize * hubble_inv * kpc * nm2 * ascale(t);
	double D = max(x_min_fix, (1 - qH) * qHe2 * qHe1 + (1 - qHe1) * qH * qHe2 + (1 - qHe2) * qH * qHe1);
	return (qH) * (1 - qHe2) * qHe1 * (1 - exp(-tau_tot)) / D;
}
