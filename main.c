#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#include "include/params.h"
#include "include/constants.h"
#include "include/proximity_zone.h"
#include "include/utils.h"
#include "include/rates.h"
#include "include/photoion.h"
#include "include/thermochemistry.h"
#include "include/evolve_cell.h"
#include "include/lya_spectra.h"
#include "include/save_outputs.h"

static FILE *open_required_file(const char *path)
{
	FILE *fp = fopen(path, "r");
	if (fp == NULL)
	{
		fprintf(stderr, "Error opening input file '%s'\n", path);
		exit(1);
	}
	return fp;
}

/*****************************************************************************/

/*
	Author: Sindhu Satyavolu
	Institution: TIFR, Mumbai
	Email: sindhusatyavolu@gmail.com
	Date: 2024

	1D Radiative Transfer Code for Proximity Zones around Quasars
	Takes as input the quasar luminosity, redshift, and the density field, temperature field and ionisation field of the IGM (includes both H and He). (Go to line 100 for input format)
	Optional inputs include the background ionising background.

	Solves the chemical rate equations for H and He, and the temperature equation. (See Section 3 of Satyavolu et al. 2023 for details)
	Computes the Lyman alpha forest spectra and the proximity zone size. (See Appendix A of Satyavolu et al. 2023 for details)

	Outputs the Lyman alpha forest spectra, ionisation and temperature fields of the IGM in the presence of quasar radiation. (Go to line 100 for output format)
	Optional outputs include the proximity zone size and the quasar lifetime.

*/

/*****************************************************************************/



/* Main routine */

int main(int argc, char **argv)
{
   
    if (argc < 5) {
        fprintf(stderr, "Usage: %s tq_Myr t_on_Myr fduty LOS_index\n", argv[0]);
        return 1;
    }

    clock_t t;
	t = clock();

	tq = atof(argv[1]);
	t_on = atof(argv[2]);
	fduty = atof(argv[3]);
	int los = atoi(argv[4]); // number of sightline specified in the command line
	if (tq <= 0.0)
	{
		fprintf(stderr, "tq_Myr must be positive, got %g\n", tq);
		return 1;
	}
	if (t_on <= 0.0)
	{
		fprintf(stderr, "t_on_Myr must be positive, got %g\n", t_on);
		return 1;
	}
	if (fduty <= 0.0 || fduty > 1.0)
	{
		fprintf(stderr, "fduty must be in (0, 1], got %g\n", fduty);
		return 1;
	}
	if (los < 0 || los >= NLOS)
	{
		fprintf(stderr, "LOS_index must be in [0, %d), got %d\n", NLOS, los);
		return 1;
	}

    printf("Quasar properties using built-in defaults\n");

	Ng = 3.15e13 * mag_to_ndot(M1450); // number of ionizing photons in Myr^-1, for e.g. 1.8444014679444228e+57 
	printf("Mag %lf Ng %e\n", M1450, Ng); 


    // age of the universe at redshift z in Myr
	double tz = (1.0 / H0) * (2.0 / 3.0) * pow(1.0 - omega_m, -0.5) * asinh(pow((1.0 - omega_m) / omega_m, 0.5) * pow(ai(z), 1.5)); 


    /* Set up the spatial grid */
	int ncells = boxsize * hubble_inv / (cellsize * hubble_inv);
	printf("ncells %d\n", ncells);


	double t_off = t_on * (1.0 / fduty - 1.0); // off time in Myr
	double *tarr = (double *)malloc(max_global_steps * sizeof(double));
	double ti = tz - tq;	   // start of simulation in Myr
	double tf = tz;			   // end of simulation in Myr
	srand((unsigned int)time(NULL));
	double tcycle = t_on + t_off;
	double t_qso = ti; //- toffsets[los];//(double)rand()/(double)(RAND_MAX) * tcycle;


	/* RT timestep using Courant condition */
	double Dt = t_on / 1e2;  // initial timestep in Myr // min(3261.6*(epsilon/0.1)*(cellsize*hubble_inv/10.0)*1e-6,t_on/1e1) ;
	

	int nt = (tf - ti) / Dt + 1; // number of timesteps

	/* Set up the temporal grid */
	double truns[nt];
	for (int i = 0; i < nt; i++)
	{
		double t0 = 0.0; // 3e-3;
		
		truns[i] = t0 + i * Dt;

	}

	/* initialize nh,nh2,NH and T grid */
	double *nh = (double *)malloc(ncells * sizeof(double));
	int r = ncells;
	int c = max_global_steps;
	double **nh2 = (double **)malloc(r * sizeof(double *));
	for (int i = 0; i < r; i++)
		nh2[i] = (double *)malloc(c * sizeof(double));
	double **nhe2 = (double **)malloc(r * sizeof(double *));
	for (int i = 0; i < r; i++)
		nhe2[i] = (double *)malloc(c * sizeof(double));
	double **nhe3 = (double **)malloc(r * sizeof(double *));
	for (int i = 0; i < r; i++)
		nhe3[i] = (double *)malloc(c * sizeof(double));
	double **NH = (double **)malloc(r * sizeof(double *));
	for (int i = 0; i < r; i++)
		NH[i] = (double *)malloc(c * sizeof(double));
	double **NHe1 = (double **)malloc(r * sizeof(double *));
	for (int i = 0; i < r; i++)
		NHe1[i] = (double *)malloc(c * sizeof(double));
	double **NHe2 = (double **)malloc(r * sizeof(double *));
	for (int i = 0; i < r; i++)
		NHe2[i] = (double *)malloc(c * sizeof(double));
	double **T = (double **)malloc(r * sizeof(double *));
	for (int i = 0; i < r; i++)
		T[i] = (double *)malloc(c * sizeof(double));
	double **dt = (double **)malloc(r * sizeof(double *));
	for (int i = 0; i < r; i++)
		dt[i] = (double *)malloc(c * sizeof(double));
	// intermediate
	double **nhe1 = (double **)malloc(r * sizeof(double *));
	for (int i = 0; i < r; i++)
		nhe1[i] = (double *)malloc(c * sizeof(double));

	double *Vpec = (double *)malloc(ncells * sizeof(double));

	/* Loading files containing densities, ionization, temperatures and background photoionization rates */
	FILE *myFile, *myFile2, *myFile3, *myFile4, *myFile5, *myFile6, *myFile7, *myFile8, *myFile9;
	char root[] = "/Users/ssatyavolu/Downloads/example_for_jahan/f1kz";
	int length, string_len1, string_len2;
	char *Bet;
	char *File1, *File2, *File3, *File4, *File5, *File6, *File7, *File8, *File9;
	length = snprintf(NULL, 0, "%0.2lf", z) + 1;
	Bet = malloc(length);
	sprintf(Bet, "%0.2lf", z);
	string_len1 = strlen(root) + strlen("M3LOS_nh.txt") + strlen(Bet) + 3;
	File1 = malloc(string_len1);
	sprintf(File1, "%s%sM3LOS_nh.txt", root, Bet);
	myFile = open_required_file(File1);
	string_len2 = strlen(root) + strlen("M3LOS_nh1.txt") + strlen(Bet) + 3;
	File2 = malloc(string_len2);
	sprintf(File2, "%s%sM3LOS_nh1.txt", root, Bet);
	myFile2 = open_required_file(File2);

	string_len2 = strlen(root) + strlen("M3LOS_T.txt") + strlen(Bet) + 3;
	File3 = malloc(string_len2);
	sprintf(File3, "%s%sM3LOS_T.txt", root, Bet);
	myFile3 = open_required_file(File3);

	string_len2 = strlen(root) + strlen("M3LOS_vpec.txt") + strlen(Bet) + 3;
	File4 = malloc(string_len2);
	sprintf(File4, "%s%sM3LOS_vpec.txt", root, Bet);
	myFile4 = open_required_file(File4);
	string_len2 = strlen(root) + strlen("M3LOS_rho.txt") + strlen(Bet) + 3;
	File5 = malloc(string_len2);
	sprintf(File5, "%s%sM3LOS_rho.txt", root, Bet);
	myFile5 = open_required_file(File5);
	string_len2 = strlen(root) + strlen("M3LOS_Gammabg.txt") + strlen(Bet) + 3;
	File6 = malloc(string_len2);
	sprintf(File6, "%s%sM3LOS_Gammabg.txt", root, Bet);
	myFile6 = open_required_file(File6);
	string_len2 = strlen(root) + strlen("M3LOS_Gammabghe1.txt") + strlen(Bet) + 3;
	File7 = malloc(string_len2);
	sprintf(File7, "%s%sM3LOS_Gammabghe1.txt", root, Bet);
	myFile7 = open_required_file(File7);
	string_len2 = strlen(root) + strlen("M3LOS_Gammabghe2.txt") + strlen(Bet) + 3;
	File8 = malloc(string_len2);
	sprintf(File8, "%s%sM3LOS_Gammabghe2.txt", root, Bet);
	myFile8 = open_required_file(File8);
	string_len2 = strlen(root) + strlen("M3LOS_Heatbg.txt") + strlen(Bet) + 3;
	File9 = malloc(string_len2);
	sprintf(File9, "%s%sM3LOS_Heatbg.txt", root, Bet);
	myFile9 = open_required_file(File9);

	/*
		myFile = fopen("/user1/sindhu/halo_los/halo_loaf1kz5.95M3LOS_nh.txt", "r");
		myFile2 = fopen("/mnt/exports/data/sindhu/halo_los/f1kz5.95M3LOS_nh1.txt", "r");
		myFile3 = fopen("/mnt/exports/data/sindhu/halo_los/f1kz5.95M3LOS_T.txt", "r");
		myFile4 = fopen("/mnt/exports/data/sindhu/halo_los/f1kz5.95M3LOS_vpec.txt", "r");
		myFile5 = fopen("/mnt/exports/data/sindhu/halo_los/f1kz5.95M3LOS_rho.txt", "r");
		myFile6 = fopen("/mnt/exports/data/sindhu/halo_los/f1kz5.95M3LOS_Gammabg.txt","r");
		myFile7 = fopen("/mnt/exports/data/sindhu/halo_los/f1kz5.95M3LOS_Gammabghe1.txt","r");
		myFile8 = fopen("/mnt/exports/data/sindhu/halo_los/f1kz5.95M3LOS_Gammabghe2.txt","r");
		myFile9 = fopen("/mnt/exports/data/sindhu/halo_los/f1kz5.95M3LOS_Heatbg.txt","r");
	*/

	/* dimensions of arrays in input files, all quantities in physical coordinates */
	static double nharr[2048][500];	  // hydrogen number density array in cm^-3, 2048 is the number of cells, 500 is the number of lines of sight
	static double nh1arr[2048][500];	  // neutral hydrogen number density array in cm^-3
	static double Tarr[2048][500];		  //  gas temperature array in K
	static double vpec[2048][500];		  //	peculiar velocity array in km/s
	static double rho[2048][500];		  // gas density array in 2e40h^2kgkpc^-3
	static double gammabg[2048][500];	  // photoionization rate of hydrogen in Myr^-1
	static double gammabghe1[2048][500]; // photoionization rate of helium 1 in Myr^-1
	static double gammabghe2[2048][500]; // photoionization rate of helium 2 in Myr^-1
	static double heatbg[2048][500];	  // background heating rate in Myr^-1

/* comment this out if you are not reading initial input from files */

	for (int i = 0; i < 2048; ++i)
	{
		for (int j = 0; j < 500; ++j)
		{
			fscanf(myFile, "%lf", &nharr[i][j]);
			fscanf(myFile2, "%lf", &nh1arr[i][j]);
			fscanf(myFile3, "%lf", &Tarr[i][j]);
			fscanf(myFile4, "%lf", &vpec[i][j]);
			fscanf(myFile5, "%lf", &rho[i][j]);
			fscanf(myFile6, "%lf", &gammabg[i][j]);
			fscanf(myFile7, "%lf", &gammabghe1[i][j]);
			fscanf(myFile8, "%lf", &gammabghe2[i][j]);
			fscanf(myFile9, "%lf", &heatbg[i][j]);
			// printf("%lf\n",nharr[i][j]);
		}
	}

	fclose(myFile);
	fclose(myFile2);
	fclose(myFile3);
	fclose(myFile4);
	fclose(myFile5);
	fclose(myFile6);
	fclose(myFile7);
	fclose(myFile8);
	fclose(myFile9);


	//#pragma omp parallel for
    //for (int los = 0; los < NLOS; ++los) {
    //    run_simulation_for_LOS(los);
    //}


	/* Initial values */
	for (int i = 0; i < ncells; ++i)
	{
		nh[i] = nharr[i][los] / X;							  // to set it to cosmological mean instead, use getnh(z) and convert to physcial coordinates;
		nh2[i][0] = nharr[i][los] / X - nh1arr[i][los] / X; // ionized hydrogen density in cm^-3

		// rescaling all densities to  zqso; ignore if the quasar is at the same redshift as the simulation (zqso=z)
		nh[i] = nh[i] * pow(1 + z_qso, 3) / pow(1 + z, 3);
		nh2[i][0] = nh2[i][0] * pow(1 + z_qso, 3) / pow(1 + z, 3);

		nhe2[i][0] = x_min_fix * (Y / (4.0 * ((1.0 - Y)))) * nh[i];							// x_min_fix*(Y/(4.0*((1.0-Y))))*nh[i]; // initial He II number density
		nhe3[i][0] = x_min_fix * (Y / (4.0 * ((1.0 - Y)))) * nh[i];							// x_min_fix*(Y/(4.0*((1.0-Y))))*nh[i]; // initial He III number density
		nhe1[i][0] = (Y / (4.0 * ((1.0 - Y)))) * nh[i] - nhe2[i][0] - nhe3[i][0];			// initial He I number density
		T[i][0] = Tarr[i][los];																// initial temperature in K
		Vpec[i] = vpec[i][los];																// peculiar velocity in km/s
		NH[i][0] = cumsum(nh, nh2, i, 0) * cellsize * hubble_inv * kpc * nm2 * ascale(ti);	// neutral hydrogen column density in nm^-2
		NHe1[i][0] = hecumsum(nhe1, i, 0) * cellsize * hubble_inv * kpc * nm2 * ascale(ti); // neutral helium column density in nm^-2
		NHe2[i][0] = hecumsum(nhe2, i, 0) * cellsize * hubble_inv * kpc * nm2 * ascale(ti); // singly ionized helium column density in nm^-2
		dt[i][0] = Dt;																		// initial time step
		Gammabkg[i] = gammabg[i][los];														// background photoionization rate of hydrogen, for reference check Haardt & Madau 2012. Assign to zero if you want to ignore background
		Gammabkg_he1[i] = gammabghe1[i][los];												// background photoionization rate of helium 1
		Gammabkg_he2[i] = gammabghe2[i][los];												// background photoionization rate of helium 2
		Heat_bkg[i] = heatbg[i][los];														// background heating rate
	}

	/*Evolving over global timestep*/
	double tsim;
	int storej = -1;
	double dt_new;
	double RIF, RIFold;
	double NHION;
	double Nphotons[nt];
	double rp[nt];

	int extra = 0; // number of cells to be included redward of the Lyman alpha line. Set to zero if you want the first cell of the sightline to be the location of the Lyman alpha line
	int encells = ncells + extra;

	double *Tau = (double *)malloc(encells * sizeof(double));
	double *NH1f = (double *)malloc(ncells * sizeof(double));
	double *Tf = (double *)malloc(ncells * sizeof(double));
	double dense;
	double UnitMass_in_g = 1.989e43;
	double Hubble_h = 0.678;
	dense = pow(3e21, 3) / (kpc * nm2 * 1.989e43 * pow(Hubble_h, 2)); //*Hubble_h*Hubble_h;

	printf("ton and toff %e %e\n", t_on, t_off);
	printf("LOS : %d\n", los);

	/* Main RT routine */
	/*For each position, solve euler to get nh(xpos,t) and T(xpos,t) where t is tsim= time after each global dt */


	for (int j = 0; j < max_global_steps; j++)
	{

		tsim = ti + Dt;
		tarr[j] = tsim;

		/*evolve all cells for a  global timestep*/
		evolve_eqns(nh, nh2, NH, nhe1, nhe2, nhe3, NHe1, NHe2, T, ncells, ti, tsim, dt, j, Dt); // my equivalent of subroutine rt_solve_cooling of Ramses RT (Rosdahl et al 2013)

		/* quasar lightcurve */
		if (dmod(tarr[j] - t_qso, tcycle) < t_on)
		{
			Ng = 3.15e13 * mag_to_ndot(M1450);
		}
		else
		{
			Ng = 0.0;
		}

		//		if (Ng==0)
		//		{
		//			printf("off\n\n");
		//		}
		//		else{
		//			printf("on\n");
		//		}

		Nphotons[j] = Ng;
		ti = tsim;

		for (int i = 0; i < ncells; ++i)
		{
			Tf[i] = T[i][j + 1];
			NH1f[i] = rho[i][los] * (1.0 - nh2[i][j + 1] / nh[i]) * (UnitMass_in_g / PROTONMASS / Hubble_h); // NH[i][storej];
		}

		//rp[j] = absorb_along_lines_of_sight(Tf, NH1f, Vpec); // lya optical depth and proximity zone size at each timestep

		if (fabs(tf - tsim) <= Dt)
		{
			printf("%lf,%lf\n", tf, tsim);
			storej = j + 1;
			for (int i = 0; i < ncells; ++i)
			{
				Tf[i] = T[i][storej];
				//  printf("%e\n",rho[i][los]);
				NH1f[i] = rho[i][los] * (1.0 - nh2[i][storej] / nh[i]) * (UnitMass_in_g / PROTONMASS / Hubble_h); // NH[i][storej];
			}

			rp[nt - 1] = absorb_along_lines_of_sight(Tf, NH1f, Vpec, Tau); //proximity zone size at the end of the simulation
			break;
		}

		for (int k = 0; k < ncells; ++k)
		{
			dt[k][j + 1] = Dt;
		}
	}

	if (storej < 0)
	{
		fprintf(stderr, "Simulation did not reach final time within %d global steps\n", max_global_steps);
		return 1;
	}

	absorb_along_lines_of_sight(Tf, NH1f, Vpec, Tau); // to get the spectrum at the end of the simulation for the last timestep

	/* Saves distance array in comoving kpc, xHI, gas temperature in K, xHeI, xHeII, xHeIII and Lyman alpha optical depth at each location at the end of the simulation for given quasar lifetime*/
    char *fileSpec;
    fileSpec = outfilename(los, M1450,t_on,fduty, z_qso,tq); // output file name
    FILE *fptr = fopen(fileSpec, "w");
    save_outputs(fptr, ncells, nh, nh2, T, nhe1, nhe2, nhe3, Tau, los, M1450, z_qso,storej);    

	t = clock() - t;
	double time_taken = ((double)t) / CLOCKS_PER_SEC; // in seconds
	printf("fun() took %f seconds to execute \n", time_taken);

	/* Optional- For saving proximity zone sizes at each timestep */

//	for (int ts = 0; ts < nt; ts++)
//	{
//		fprintf(fptr, "%0.12e %e %0.12e\n", truns[ts], Nphotons[ts], rp[ts]);
//	}

/*comment out if you are not reading quasar properties from a file and want to directly assign redshift and magnitude */
	/*FILE *myFile10;
	myFile10 = fopen("/user1/sindhu/example_for_jahan/qsopropsmfp.txt", "r"); // file that includes the list of quasar magnitudes and redshifts
	double **qprops = (double **)malloc(31 * sizeof(int *));

	for (int i = 0; i < 1; i++)
		qprops[i] = (double *)malloc(2 * sizeof(int));

	for (int l = 0; l < 1; l++)
	{
		for (int k = 0; k < 2; k++)
		{
			fscanf(myFile10, "%lf", &qprops[l][k]);
		}
	}
    //double zqso = qprops[m][0]; // redshift of the quasar ; equal to z if the quasar is at the same redshift as the simulation
    */

}
