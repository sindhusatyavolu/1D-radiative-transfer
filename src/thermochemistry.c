#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../include/thermochemistry.h"
#include "../include/utils.h"
#include "../include/rates.h"
#include "../include/photoion.h"
#include "../include/constants.h"
#include "../include/params.h"


/* Routines to solve chemical rate equation during global time step dx/dt = C*x + D */

double xh2solve(double Tu, double xh2, double xhe1, double xhe2, double xhe3, double dti,
				double NH1, double Nhe1, double Nhe2, double R, double nhi, double *xrate, double t, double gbkg)
{
	double nhe1l = max((Y / (4.0 * ((1.0 - Y)))) * nhi * (x_min_fix), (Y / (4.0 * ((1.0 - Y)))) * nhi * (xhe1));
	double nhe2l = max((Y / (4.0 * ((1.0 - Y)))) * nhi * xhe2, (Y / (4.0 * ((1.0 - Y)))) * nhi * x_min_fix);
	double nhe3l = max((Y / (4.0 * ((1.0 - Y)))) * nhi * xhe3, (Y / (4.0 * ((1.0 - Y)))) * nhi * x_min_fix);
	double ne = xh2 * nhi + (Y / (4.0 * ((1.0 - Y)))) * nhi * (xhe2 + 2 * xhe3);
	double xe = xh2 + (xhe2 + 2 * xhe3);
	double nu = (1.0 + 4.0 * Y) / (1.0 + Y + xe); // bolton
	double Ti = (2.0 / 3.0) * Tu * (nu * PROTONMASS) / (kb);
	double fac = -alpha * Ng / (pow(3e21, 3) * 13.6); // Myr^-1 eV^-1

	/*isothermal case*/
	// double fac = Ng/pow(3e21,3);

	double vol = 1.0 / (4 * PI * R * R * cellsize * hubble_inv * ascale(t)); // kpc^-3 converted to cm^-3 through fac
	double nh1 = max(nhi * x_min_fix, nhi * (1 - xh2));
	double C, D, dbdt, dalphadt, Tnu, dCdx, dDdx, J;
	double Gammab;

	// if (background==1)
	//{
	Gammab = gbkg * timefac; // from haardt madau 2012
							 //}
							 // else
	//	Gammab = 0;//0.23e-12*timefac; // from haardt madau 2012
	//		Gammab = timefac*1.27e-11*Jb*pow(-alpha+3,-1);

	C = GE(Ti) * ne + fac * vol * Gamma(NH1, Nhe1, Nhe2, nh1, nhe1l, nhe2l, t) / (nh1) + Gammab;
	D = ralpha(Ti) * ne;

	if (C != C)
	{
		printf("H1 gone wrong\n");
		printf("gamma %0.12e %0.12e %0.12e\n", fac * vol, Gamma(NH1, Nhe1, Nhe2, nh1, nhe1l, nhe2l, t), 1 / (nh1));
		printf("elec coll %0.12e %0.12e\n", GE(Ti), ne);
		exit(0);
	}

	dbdt = -2.91364 * exp(-157809.1 / Ti) / pow((1 + 0.00316228 * pow(Ti, 0.5)), 2) + 2.90803 * exp(-157809.1 / Ti) / ((1 + 0.00316228 * pow(Ti, 0.5)) * pow(Ti, 1.5)) + 921.375 * exp(-157809.1 / Ti) / ((1 + 0.00316228 * pow(Ti, 0.5)) * pow(Ti, 0.5));

	dalphadt = -(1.10726 * 1e9 * pow((1 / Ti), 2.503)) / (1 + 521.553 * pow((1 / Ti), 0.47)) * 1.923 + (3.4727 * 1e11 * pow((1 / Ti), 2.973)) / pow((1 + 521.553 * pow((1 / Ti), 0.47)), 2.923);

	nu = (1.0 + 4.0 * Y) / (1.0 + Y + ne); // bolton
	//	nu = pow((X*(1+xh2)+Y*(1+xhe2+2*xhe3)/4),(-1));
	Tnu = Ti / nu;
	dCdx = nhi * GE(Ti) + ne * Tnu * pow(nu, 2) * X * dbdt;
	dDdx = nhi * ralpha(Ti) + ne * Tnu * pow(nu, 2) * X * dalphadt;
	J = dCdx - (C + D) - xh2 * (dCdx + dDdx);
	*(xrate) = C - xh2 * (C + D);														 //(double) max(0.0,(C-xh2*(C+D)));
																						 /*
																							 if ((1-J*dti)==0 || *(xrate)<0)
																							 {
																								 printf("x denom is zero\n");
																								 printf("%0.12e %0.12e %0.12e %0.12e %0.12e negative dervative\n",(C-xh2*(C+D)),xh2,C,D,xh2*(C+D) );
																								 exit(0);
																							 }
																						 */
	return min(1.0 - x_min_fix, max(x_min_fix, xh2 + dti * (*(xrate)) / (1 - J * dti))); // RAMSES implicit
																						 // 	return  (xh2 + C*dti)/(1+D*dti);  // Euler implicit
}

double xhe1solve(double Tu, double xh2i, double xhe1, double xhe2, double xhe3, double dti,
				 double NH1, double Nhe1, double Nhe2, double R, double nhi, double t, double gbkg_he1)
{
	double nhe1l = max((Y / (4.0 * ((1.0 - Y)))) * nhi * (x_min_fix), (Y / (4.0 * ((1.0 - Y)))) * nhi * (xhe1));
	double nhe2l = max((Y / (4.0 * ((1.0 - Y)))) * nhi * xhe2, (Y / (4.0 * ((1.0 - Y)))) * nhi * x_min_fix); //(Y/(4*((1.0-Y))))*nhi*xhe2;
	double nhe3l = max((Y / (4.0 * ((1.0 - Y)))) * nhi * xhe3, (Y / (4.0 * ((1.0 - Y)))) * nhi * x_min_fix); //(Y/(4*((1.0-Y))))*nhi*xhe3;
	//	double xhe1 = 1- xhe2 -xhe3;
	double ne = xh2i * nhi + (Y / (4.0 * ((1.0 - Y)))) * nhi * (xhe2 + 2 * xhe3);
	double xe = xh2i + (xhe2 + 2 * xhe3);
	double nu = (1.0 + 4.0 * Y) / (1.0 + Y + xe); // bolton
	double Ti = (2.0 / 3.0) * Tu * (nu * PROTONMASS) / (kb);
	double fac = -alpha * Ng / (pow(3e21, 3) * 13.6);

	/*isothermal case */
	// double fac = Ng/pow(3e21,3);

	double vol = 1.0 / (4 * PI * R * R * cellsize * hubble_inv * ascale(t));
	double nh1 = max(nhi * (1 - xh2i), nhi * x_min_fix); // nhi*(1-xh2i);
	double C, D, Tnu;
	//	nu = pow((X*(1+xh2i)+Y*(1+xhe2+2*xhe3)/4),(-1));
	nu = (1.0 + 4.0 * Y) / (1.0 + Y + ne); // bolton
	//	Tnu = Ti/nu;
	double Gammab;
	// if (background==1)
	//{
	Gammab = gbkg_he1 * timefac; // 0.154e-12*timefac; // from haardt madau 2012
	//}
	// else
	//	Gammab = 0;//0.154e-12*timefac; // from haardt madau 2012

	C = xhe2 * alphahe2(Ti) * ne;
	//	D = betahe2(Ti)*ne + alphahe2(Ti)*ne + fac*vol*Gammahe2(NH1,Nhe1,Nhe2,nh1,nhe1,nhe2)/(nhe2) ;
	D = betahe1(Ti) * ne + fac * vol * Gammahe1(NH1, Nhe1, Nhe2, nh1, nhe1l, nhe2l, t) / nhe1l + Gammab;

	if (D != D)
	{
		printf("He1 gone wrong\n");
		printf("gamma %0.12e %0.12e %0.12e\n", fac * vol, Gammahe1(NH1, Nhe1, Nhe2, nh1, nhe1l, nhe2l, t), 1 / (nhe1l));
		printf("elec coll %0.12e,%0.12e, %0.12e\n", xhe2, betahe1(Ti), ne);
		exit(0);
	}

	// printf("C %e D %e\n",C,D );

	if ((nhe1l + D * dti) == 0)
	{
		printf("xhe1 denom is zero\n");
		exit(0);
	}

	return min(1.0 - x_min_fix, max(x_min_fix, (xhe1 + C * dti) / (1 + D * dti)));
	//      return xhe1 + (C-D*xhe1)*dti;
}

double xhe2solve(double Tu, double xh2i, double xhe1, double xhe2, double xhe3, double dti,
				 double NH1, double Nhe1, double Nhe2, double R, double nhi, double *hrate, double t, double gbkg_he1, double gbkg_he2)
{
	//	double nhe1l = (Y/(4*((1.0-Y))))*nhi*(xhe1);
	double nhe1l = max((Y / (4.0 * ((1.0 - Y)))) * nhi * (x_min_fix), (Y / (4.0 * ((1.0 - Y)))) * nhi * (xhe1));
	double nhe2l = max((Y / (4.0 * ((1.0 - Y)))) * nhi * xhe2, (Y / (4.0 * ((1.0 - Y)))) * nhi * x_min_fix);
	double nhe3l = max((Y / (4.0 * ((1.0 - Y)))) * nhi * xhe3, (Y / (4.0 * ((1.0 - Y)))) * nhi * x_min_fix); //(Y/(4*((1.0-Y))))*nhi*xhe3;
	//	double xhe1 = 1- xhe2 -xhe3;
	double ne = xh2i * nhi + (Y / (4.0 * ((1.0 - Y)))) * nhi * (xhe2 + 2 * xhe3);
	double xe = xh2i + (xhe2 + 2 * xhe3);
	double nu = (1.0 + 4.0 * Y) / (1.0 + Y + xe); // bolton
	double Ti = (2.0 / 3.0) * Tu * (nu * PROTONMASS) / (kb);
	double fac = -alpha * Ng / (pow(3e21, 3) * 13.6);

	/*isothermal case */
	// double fac = Ng/pow(3e21,3);

	double vol = 1.0 / (4 * PI * R * R * cellsize * hubble_inv * ascale(t));
	double nh1 = max(nhi * (1 - xh2i), nhi * x_min_fix);
	double C, D, Tnu;
	//	nu = pow((X*(1+xh2i)+Y*(1+xhe2+2*xhe3)/4),(-1));
	nu = (1.0 + 4.0 * Y) / (1.0 + Y + ne); // bolton
	//	Tnu = Ti/nu;
	double Gammab1, Gammab2;

	Gammab1 = gbkg_he1 * timefac; // 0.154e-12*timefac;
	Gammab2 = gbkg_he2 * timefac; // 0.442e-18*timefac; // from haardt madau 2012

	C = xhe1 * ne * betahe1(Ti) + fac * vol * Gammahe1(NH1, Nhe1, Nhe2, nh1, nhe1l, nhe2l, t) / ((Y / (4.0 * ((1.0 - Y)))) * nhi) + xhe3 * alphahe3(Ti) * ne + xhe1 * Gammab1;
	D = betahe2(Ti) * ne + alphahe2(Ti) * ne + fac * vol * Gammahe2(NH1, Nhe1, Nhe2, nh1, nhe1l, nhe2l, t) / (nhe2l) + Gammab2;
	//	D = alphahe2(Ti)*ne;

	if (C != C)
	{
		printf("He2 gone wrong\n");
		printf("gamma %0.12e %0.12e %0.12e\n", fac * vol, Gammahe1(NH1, Nhe1, Nhe2, nh1, nhe1l, nhe2l, t), 1 / (nhe2l));
		printf("elec coll %0.12e,%0.12e, %0.12e\n", xhe1, betahe1(Ti), ne);
		exit(0);
	}

	*(hrate) = (C - D * xhe2);

	if ((1 + D * dti) == 0)
	{
		printf("xhe2 denom is zero\n");
		exit(0);
	}
	return min(1.0 - x_min_fix, max(x_min_fix, (xhe2 + C * dti) / (1 + D * dti)));
	//	return xhe2 + (C-D*xhe2)*dti;
}

double xhe3solve(double Tu, double xh2i, double xhe1, double xhe2, double xhe3, double dti,
				 double NH1, double Nhe1, double Nhe2, double R, double nhi, double *hrate, double t, double gbkg_he2)
{
	double nhe1l = max((Y / (4.0 * ((1.0 - Y)))) * nhi * (x_min_fix), (Y / (4.0 * ((1.0 - Y)))) * nhi * (xhe1));
	double nhe2l = max((Y / (4.0 * ((1.0 - Y)))) * nhi * xhe2, (Y / (4.0 * ((1.0 - Y)))) * nhi * x_min_fix);
	double nhe3l = max((Y / (4.0 * ((1.0 - Y)))) * nhi * xhe3, (Y / (4.0 * ((1.0 - Y)))) * nhi * x_min_fix);
	double ne = xh2i * nhi + (Y / (4.0 * ((1.0 - Y)))) * nhi * (xhe2 + 2 * xhe3);
	double xe = xh2i + (xhe2 + 2 * xhe3);
	double nu = (1.0 + 4.0 * Y) / (1.0 + Y + xe); // bolton
	double Ti = (2.0 / 3.0) * Tu * (nu * PROTONMASS) / (kb);
	double fac = -alpha * Ng / (pow(3e21, 3) * 13.6);

	/*isothermal case */
	// double fac = Ng/pow(3e21,3);

	double vol = 1.0 / (4 * PI * R * R * cellsize * hubble_inv * ascale(t));
	//	double nh1 = nhi*(1-xh2i);
	double nh1 = max(nhi * (1 - xh2i), nhi * x_min_fix);
	double C, D, Tnu;
	//	nu = pow((X*(1+xh2i)+Y*(1+xhe2+2*xhe3)/4.0),(-1));
	nu = (1.0 + 4.0 * Y) / (1.0 + Y + ne); // bolton
	//	Tnu = Ti/nu;
	double Gammab;

	Gammab = gbkg_he2 * timefac; // 0.442e-18*timefac; // from haardt madau 2012

	C = xhe2 * ne * betahe2(Ti) + fac * vol * Gammahe2(NH1, Nhe1, Nhe2, nh1, nhe1l, nhe2l, t) / ((Y / (4.0 * ((1.0 - Y)))) * nhi) + xhe2 * Gammab;
	D = alphahe3(Ti) * ne;

	*(hrate) = (C - D * xhe3);

	if (C != C)
	{
		printf("He3 gone wrong\n");
		printf("gamma %0.12e %0.12e %0.12e \n", fac * vol, Gammahe2(NH1, Nhe1, Nhe2, nh1, nhe1l, nhe2l, t), xhe2 / (nhe2l));
		printf("elec coll %0.12e %0.12e, %0.12e\n", xhe2, ne, betahe2(Ti));
		exit(0);
	}

	if ((1 + D * dti) == 0)
	{
		printf("xhe3 denom is zero\n");
		exit(0);
	}

	return min(1.0 - x_min_fix, max(x_min_fix, (xhe3 + C * dti) / (1 + D * dti)));
	//	return xhe3 + (C-D*xhe3)*dti;
}

/* Routine for solving dT_nu/dt using an explicit scheme */

double Tsolve(double Tu, double NH1, double Nhe1, double Nhe2,
			  double xh2i, double xhe1, double xhe2, double xhe3,
			  double nhi, double dti, double R, double tsim,
			  double *trate, double *xrate, double t, double gbkg, double gbkg_he1, double gbkg_he2, double heatbg)
{
	//	 printf("NG %0.12e\n",Ng);
	double fac = -alpha * Ng / (pow(3e21, 3) * 13.6);
	/* isothermal case */
	//	double fac = Ng/(pow(3e21,3));
	double fh = fac * 1.6e-12; // 1.6e-12 is factor ev to erg
	double vol = 1.0 / (4 * PI * R * R * cellsize * hubble_inv * ascale(t));
	//	double K = (2.0/3)*(1.0/(nhi*kb)); // only H # cm^3erg^-1 K
	double nh1 = max(nhi * (1 - xh2i), nhi * x_min_fix);
	double nhe1l = max((Y / (4.0 * ((1.0 - Y)))) * nhi * (x_min_fix), (Y / (4.0 * ((1.0 - Y)))) * nhi * (xhe1));
	double nhe2l = max((Y / (4.0 * ((1.0 - Y)))) * nhi * xhe2, (Y / (4.0 * ((1.0 - Y)))) * nhi * x_min_fix);
	double nhe3l = max((Y / (4.0 * ((1.0 - Y)))) * nhi * xhe3, (Y / (4.0 * ((1.0 - Y)))) * nhi * x_min_fix);
	double ne = xh2i * nhi + (Y / (4.0 * ((1.0 - Y)))) * nhi * (xhe2 + 2 * xhe3);
	double xe = xh2i + (xhe2 + 2 * xhe3);
	double n = nhi + (Y / (4.0 * ((1.0 - Y)))) * nhi + ne; // total species number density
	//	double nu = pow((X*(1+xh2i)+Y*(1+xhe2+2*xhe3)/4),(-1));
	double nu = (1.0 + 4.0 * Y) / (1.0 + Y + xe); // bolton
	double Ti = (2.0 / 3.0) * Tu * (nu * PROTONMASS) / (kb);
	double K = (2.0 / 3.0) * (1.0 / (n * kb)); // From Chen in units cm^3 erg^-1 K
	double T5 = Ti / 100000.0;
	double f1 = 1.778e-29 * timefac; // #*6.24e11 # in erg cm^3 Myr^-1 K^-1 #4.05   #
	double f2 = 7.5e-19 * timefac;	 // #*6.24e11 # in erg cm^3 Myr^-1 #1.711e+11
	double f3 = 1.42e-27 * timefac;	 // #*6.24e11# in erg cm^3 K^-1 #489.62
	double gff = 1.5;				 // # from cen #1.1+ 0.34*(np.exp(-5.5+np.log(T)))**2
	double f4 = 5.65e-36 * timefac;	 // #*6.24e11 #in erg cm^3 Myr^-1 K^-1 #1.289e-6
	double f5 = 1.27e-21 * timefac;
	double f11 = 1.0;
	double f12 = 8 * 1.778e-29 * timefac;
	double f21 = timefac * 5.54e-17;
	double f51 = timefac * 9.38e-22;
	double f52 = timefac * 4.95e-22;
	/*
	   f11 = 1
	   f12 =8*1.778*1e-29*timefac/(1.38e-16)
	   f21 =5.54e-17*timefac/(1.38e-16)
	   f22 =9.1e-27*timefac/(1.38e-16)
	   */
	double Th1 = 157807;
	//	double The = 631515;
	double lamb = 2 * Th1 / Ti;
	double lamb2 = 1263030.0 / Ti;
	double omega_m = 0.308;
	double H0 = 67.8 * 1.05e-6;
	double rho = 1.0;
	//	fh = 1.04e+8
	double Hf = 1.05e-6;
	double omega_b = 0.0482;
	double rhocb = 3.0 * omega_b * pow(H0 / timefac, 2.0) / (8.0 * PI * GRAVITY);
	double atime = ascale(t);
	double DELTA = 0.0;
	double RHOb = rhocb / pow(atime, 3.0);
	double Xe = ne * PROTONMASS / (RHOb * (1.0 + DELTA));
	/* To turn off helium */
	//      double X = 1;
	//	double Y = 0;
	//	double xhe2 = 0;
	//	double xhe3 = 0;

	if (K * nu == 0 || nu != nu)
	{
		printf("Knu is zero\n");
		printf("nu %0.12e, %0.12e,%0.12e,%0.12e\n", nu, xh2i, xhe2, xhe3);
		exit(0);
	}

	double zz = 1.0 / ascale(t) - 1.0;
	double c1 = (f1)*Ti * pow(lamb, 1.965) / pow((1.0 + pow((lamb / 0.541), 0.502)), 2.697); // recombination cooling coeffecient erg cm^3 Myr^-1 from HnG
	double c2 = (f2)*pow((1.0 + pow(T5, 1.0 / 2.0)), -1) * exp(-118355.0 / Ti);				 // collisional excitation cooling coeffecient
	double c3 = gff * (f3)*pow(Ti, 0.5);													 // Brehmstrallung/ free free emission
	//	double c4 = (2.726*(1.0+zz)-Ti)*6.35e-41*timefac*omega_b*0.678*0.678*Xe*pow(1.0+zz,7);
	double c4 = (f4) * (Ti - 2.73 * (1.0 + z)) * pow(1.0 + z, 4.0); // inverse compton
	double c5 = kb * Th1 * GE(Ti);									// from HnG in erg cm^3 Myr^-1
	//      f5*pow(Ti,0.5)*pow((1.0+pow(T5,0.5)),-1)*exp(-157809.1/Ti);// collisional Ionisation cooling
	double c11;
	if (5e3 < Ti && Ti < 5e5)
	{
		c11 = f11 * kb * Ti * alphahe2(Ti); // recombination cooling He II erg cm^3 Myr^-1
	}
	else
		c11 = 0.0;

	double c12 = f12 * Ti * pow(lamb2, 1.965) / pow((1.0 + pow((lamb2 / 0.541), 0.502)), 2.697);	// recombination cooling He III
	double c21 = f21 * pow(Ti, -0.397) * pow((1.0 + pow(T5, 1.0 / 2.0)), -1) * exp(-473638.0 / Ti); // collisional excitation HeII erg cm^3 Myr^-1
	double c51 = kb * 285335.0 * betahe1(Ti);														// f51*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-285335.4/Ti);
	double c52 = kb * 631515.0 * betahe2(Ti);														// f52*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-631515.0/Ti);
	double C;
	/* Total cooling term Lambda */
	C = c1 * ne * xh2i * nhi + c11 * ne * nhe2l + c12 * ne * nhe3l + c2 * ne * nh1 + c21 * ne * nhe2l + c3 * (nhi * xh2i + nhe2l + 4 * nhe3l) * ne + c4 * ne + c5 * ne * nh1 + c51 * ne * nhe1l + c52 * ne * nhe2l + 2 * Hubble(tsim) * Ti / (K) + *(xrate) * (Ti / (n)) / (K);

	if (C != C)
	{
		printf("Cooling is gone bonkers\n");
		printf("c1 %0.12e\n", c1);
		printf("c11 %0.12e\n", c11);
		printf("c12 %0.12e\n", c12);
		printf("c2 %0.12e\n", c2);
		printf("c21 %0.12e\n", c21);
		printf("c3 %0.12e\n", c3);
		printf("c4 %0.12e\n", c4);
		printf("c5 %0.12e\n", c5);
		printf("c51 %0.12e\n", c51);
		printf("c52 %0.12e\n", c52);
		printf("H %0.12e, %0.12e, %0.12e\n", 2 * Hubble(tsim), Ti, (K * nu));
		exit(0);
	}
	double bheat1, bheat2, bheat3;
	if (background == 1)
	{
		bheat1 = 0; // 13.6*timefac*gbkg*1.6e-12;//0.989e-12*timefac*1.6e-12; // background photoheating from haardt madau 2012 erg Myr^-1
		bheat2 = 0; // 0.112e-11*timefac*1.6e-12;
		bheat3 = 0; // 0.269e-16*timefac*1.6e-12;
	}
	else
	{
		bheat1 = 0;
		bheat2 = 0;
		bheat3 = 0;
	}
	// double L = (fh*vol*Heat(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l,t)+ nh1*bheat1
	//		+fh*vol*Heat2(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l,t)+nhe1l*bheat2
	//		+fh*vol*Heat3(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l,t)+nhe3l*bheat3
	//		-C);
	double L = (fh * vol * Heat(NH1, Nhe1, Nhe2, nh1, nhe1l, nhe2l, t) + fh * vol * Heat2(NH1, Nhe1, Nhe2, nh1, nhe1l, nhe2l, t) + fh * vol * Heat3(NH1, Nhe1, Nhe2, nh1, nhe1l, nhe2l, t) + heatbg - C);

	// derivatives of cooling processes
	double A1 = -(6.1712 * 1e10 * f1 * pow((1 / Ti), 1.965)) / pow((1 + 784.353 * pow((1 / Ti), 0.502)), 2.697) + (6.79107 * 1e13 * f1 * pow((1 / Ti), 2.467)) / pow((1 + 784.353 * pow((1 / Ti), 0.502)), 3.697);
	double A2 = (118355 * exp(-118355 / Ti) * f2) / ((1 + 0.00316228 * pow(Ti, 0.5)) * pow(Ti, 2)) - (0.00158114 * exp(-118355 / Ti) * f2) / (pow((1 + 0.00316228 * pow(Ti, 0.5)), 2) * pow(Ti, 0.5));
	double A3 = (0.5 * f3 * gff) / pow(Ti, 0.5);
	double A4 = f4 * pow((1 + z), 4);
	double A11 = -4.96095e-13 * pow(1 / Ti, 1.654);
	double A12 = -(9.41467e11 * f12 * pow(1 / Ti, 1.965) / (pow((1 + 1573.42 * pow((1 / Ti), 0.502)), 2.697))) + 2.07829e15 * f12 * pow(1 / Ti, 2.467) / pow((1 + 1573.42 * pow((1 / Ti), 0.502)), 3.697);
	double A21 = 473638.0 * exp(-473638 / Ti) * f21 / ((1 + pow(T5, 1.0 / 2)) * pow(Ti, 2.397)) - 0.397 * exp(-473638 / Ti) * f21 / ((1 + pow(T5, 1.0 / 2)) * pow(Ti, 1.397)) - exp(-473638 / Ti) * f21 / (200 * pow(10, 0.5) * pow((1 + pow(T5, 1.0 / 2)), 2) * pow(Ti, 0.897));
	/*Add coll ion derivatives*/
	double dLdt = -(A1 + A2 + A3 + A4 + A11 + A12 + A21); // +2*Hubble(tsim)/(K*nu));//+(1/(nhi*xh2i))*(*xrate));
	//	*(trate) = L*K;
	*(trate) = L / (nhi * PROTONMASS);
	//      return Ti+ dti*L*K;
	//	return Tu + dti*(L/(nhi*PROTONMASS))/(1-(dLdt/(nhi*PROTONMASS))*dti) ;
	//	return (L*K*nu*dti + Ti)/(1+2*Hubble(tsim)*dti) ; // anninos implicit
	return Tu + dti * L * X / (nhi * PROTONMASS);
}

/* Routine for solving the rate equations for the chemistry and temperature over global timestep */
long int xTiter(double dti, double NH1, double Nhe1, double Nhe2, double R, double nhi,
				double tsim, double ti, double *Tt, double *xh2t, double *xhe1t, double *xhe2t, double *xhe3t, double *dtrec, double gbkg, double gbkg_he1, double gbkg_he2, double heatbg)
{
	int storei = 0;
	double frac_change = 0.0;
	// double ne_init;
	double dUU;
	double Trate = 0;
	double Xrate = 0;
	double he2rate = 0;
	double he3rate = 0;
	double nerate = 0;
	double ne_init, ne_up;
	double mu;
	double tcal;
	tcal = ti;
	for (int i = 0; i < max_subcycling_steps; i++)
	{
		Tt[i + 1] = (double)Tsolve(Tt[i], NH1, Nhe1, Nhe2, xh2t[i], xhe1t[i], xhe2t[i], xhe3t[i], nhi, dti, R, tsim, &Trate, &nerate, tcal, gbkg, gbkg_he1, gbkg_he2, heatbg);
		dUU = (double)fabs(max(T2_MIN_FIX, Tt[i] + Trate * dti) - Tt[i]);

		/* Sub cycling condition */
		if ((double)max(dUU, fabs(Tt[i + 1] - Tt[i])) / (Tt[i] + T_MIN) > lim && dti > dt_min)
		{
			return 100000000000;
		}
		else
		{
			frac_change = (double)max(frac_change, dUU);
		}

		xh2t[i + 1] = (double)max(x_min_fix, xh2solve(Tt[i + 1], xh2t[i], xhe1t[i], xhe2t[i], xhe3t[i], dti, NH1, Nhe1, Nhe2, R, nhi, &Xrate, tcal, gbkg));

		/*complete this exactly?*/
		if ((double)max(fabs(xh2t[i + 1] - xh2t[i]), Xrate * dti) / (xh2t[i] + X_FM) > lim && dti > dt_min)
		{
			return 100000000001;
		}
		else
		{
			frac_change = (double)max(frac_change, fabs(xh2t[i + 1] - xh2t[i]) / (xh2t[i])); // check
		}

		//		xhe1t[i+1] = (double) max(x_min_fix,1-xhe2t[i+1]-xhe3t[i+1]);
		xhe1t[i + 1] = (double)max(x_min_fix, xhe1solve(Tt[i + 1], xh2t[i + 1], xhe1t[i], xhe2t[i], xhe3t[i], dti, NH1, Nhe1, Nhe2, R, nhi, tcal, gbkg_he1));

		xhe2t[i + 1] = (double)max(x_min_fix, xhe2solve(Tt[i + 1], xh2t[i + 1], xhe1t[i], xhe2t[i], xhe3t[i], dti, NH1, Nhe1, Nhe2, R, nhi, &he2rate, tcal, gbkg_he1, gbkg_he2));

		xhe3t[i + 1] = (double)max(x_min_fix, xhe3solve(Tt[i + 1], xh2t[i + 1], xhe1t[i], xhe2t[i + 1], xhe3t[i], dti, NH1, Nhe1, Nhe2, R, nhi, &he3rate, tcal, gbkg_he2));

		if (xhe1t[i + 1] + xhe2t[i + 1] + xhe3t[i + 1] != 1e0)
		{
			if (max(max(xhe2t[i + 1], xhe3t[i + 1]), xhe1t[i + 1]) == xhe2t[i + 1])
			{
				xhe2t[i + 1] = 1 - xhe1t[i + 1] - xhe3t[i + 1];
			}
			else if (max(max(xhe2t[i + 1], xhe3t[i + 1]), xhe1t[i + 1]) == xhe3t[i + 1])
			{
				xhe3t[i + 1] = 1 - xhe1t[i + 1] - xhe2t[i + 1];
			}
		}

		ne_init = xh2t[i] * nhi + (Y / (4.0 * ((1.0 - Y)))) * nhi * (xhe2t[i] + 2 * xhe3t[i]);
		ne_up = xh2t[i + 1] * nhi + (Y / (4.0 * ((1.0 - Y)))) * nhi * (xhe2t[i + 1] + 2 * xhe3t[i + 1]);
		nerate = Xrate * nhi + (Y / (4.0 * ((1.0 - Y)))) * nhi * (he2rate + 2 * he3rate);

		if (fabs(ne_up - ne_init) / (ne_init + X_FM) > lim && dti > dt_min)
		{
			// printf("change %0.12e %0.12e\n", fabs(ne_up - ne_init)/(ne_init+X_FM), lim);
			return 100000000002;
		}
		else
		{
			frac_change = (double)max(frac_change, fabs(ne_up - ne_init) / (ne_init));
		}

		tcal = ti + dti;
		ti = tcal;

		/* Sub cycling time step for next loop */
		if (frac_change < 0.5 && dti > dt_min)
		{ // has value 0.05 in paper but 0.5 here
			*(dtrec) = 2 * dti;
		}
		else
		{
			*(dtrec) = dti;
		}

		if (tcal >= tsim)
		{
			storei = i + 1;
			break;
		}
		dti = (double)min(*(dtrec), tsim - tcal);
	}
	double nu = (1.0 + 4.0 * Y) / (1.0 + Y + xh2t[storei] + xhe2t[storei] + 2 * xhe3t[storei]);
	//      printf("T %lf %d\n",(2.0/3.0)*Tt[storei]*(nu*PROTONMASS)/(kb),storei);
	return storei;
}
