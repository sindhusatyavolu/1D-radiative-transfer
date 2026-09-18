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

static char *input_filename(const char *input_dir, const char *prefix, const char *suffix)
{
	int length = snprintf(NULL, 0, "%s/%s_%s.txt", input_dir, prefix, suffix) + 1;
	char *path = malloc(length);
	if (path == NULL)
	{
		fprintf(stderr, "Failed to allocate input filename\n");
		exit(1);
	}
	sprintf(path, "%s/%s_%s.txt", input_dir, prefix, suffix);
	return path;
}

static void read_grid_value(FILE *fp, const char *path, double *value)
{
	if (fscanf(fp, "%lf", value) != 1)
	{
		fprintf(stderr, "Error reading numeric value from '%s'\n", path);
		exit(1);
	}
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
        fprintf(stderr, "Usage: %s tq_Myr t_on_Myr fduty LOS_index [M1450 runtime_z z_qso input_dir input_prefix output_dir n_cells n_los output_los legacy_density num_pixels_redward legacy_alphahe2_bug boxsize_override comoving_length_scale initial_only density_multiplier]\n", argv[0]);
        return 1;
    }

    clock_t t;
	t = clock();

	tq = atof(argv[1]);
	t_on = atof(argv[2]);
	fduty = atof(argv[3]);
	int los = atoi(argv[4]); // number of sightline specified in the command line
	const char *input_dir = "example_input_files";
	const char *input_prefix = NULL;
	const char *output_dir = ".";

	if (argc > 5)
		M1450 = atof(argv[5]);
	if (argc > 6)
		z = atof(argv[6]);
	if (argc > 7)
		z_qso = atof(argv[7]);
	if (argc > 8)
		input_dir = argv[8];
	if (argc > 9)
		input_prefix = argv[9];
	if (argc > 10)
		output_dir = argv[10];
	if (argc > 11)
		Ncells = atoi(argv[11]);
	if (argc > 12)
		NLOS = atoi(argv[12]);
	int output_los = los;
	if (argc > 13)
		output_los = atoi(argv[13]);
	int legacy_density = 0;
	if (argc > 14)
		legacy_density = atoi(argv[14]);
	if (argc > 15)
		num_pixels_redward = atoi(argv[15]);
	if (argc > 16)
		legacy_alphahe2_bug = atoi(argv[16]);
	double boxsize_override = 0.0;
	if (argc > 17)
		boxsize_override = atof(argv[17]);
	double comoving_length_scale = 1.0;
	if (argc > 18)
		comoving_length_scale = atof(argv[18]);
	int initial_only = 0;
	if (argc > 19)
		initial_only = atoi(argv[19]);
	double density_multiplier = 1.0;
	if (argc > 20)
		density_multiplier = atof(argv[20]);
	if (initial_only != 0 && initial_only != 1)
	{
		fprintf(stderr, "initial_only must be 0 or 1, got %d\n", initial_only);
		return 1;
	}
	if (!isfinite(comoving_length_scale) || comoving_length_scale <= 0.0)
	{
		fprintf(stderr, "comoving_length_scale must be finite and positive, got %g\n", comoving_length_scale);
		return 1;
	}
	if (!isfinite(density_multiplier) || density_multiplier <= 0.0)
	{
		fprintf(stderr, "density_multiplier must be finite and positive, got %g\n", density_multiplier);
		return 1;
	}
	if (!isfinite(tq) || tq < 0.0 || (!initial_only && tq == 0.0))
	{
		fprintf(stderr, "tq_Myr must be positive (or zero for initial_only), got %g\n", tq);
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
	if (Ncells <= 0 || Ncells > RTE_MAX_CELLS)
	{
		fprintf(stderr, "n_cells must be in [1, %d], got %d\n", RTE_MAX_CELLS, Ncells);
		return 1;
	}
	if (NLOS <= 0 || NLOS > 1000)
	{
		fprintf(stderr, "n_los must be in [1, 1000], got %d\n", NLOS);
		return 1;
	}
	if (num_pixels_redward < 0)
	{
		fprintf(stderr, "num_pixels_redward must be non-negative, got %d\n", num_pixels_redward);
		return 1;
	}
	if (boxsize_override < 0.0)
	{
		fprintf(stderr, "boxsize_override must be non-negative, got %g\n", boxsize_override);
		return 1;
	}

	char default_prefix[64];
	if (input_prefix == NULL)
	{
		snprintf(default_prefix, sizeof(default_prefix), "f1kz%0.2lfM3LOS", z);
		input_prefix = default_prefix;
	}

    printf("Quasar/input properties\n");

	Ng = initial_only ? 0.0 : 3.15e13 * mag_to_ndot(M1450); // number of ionizing photons in Myr^-1
	printf("Mag %lf Ng %e\n", M1450, Ng); 
	printf("z_sim %lf z_qso %lf input %s/%s_*.txt\n", z, z_qso, input_dir, input_prefix);
	printf("n_cells %d n_los %d input_los %d output_los %d legacy_density %d density_multiplier %.17g num_pixels_redward %d\n", Ncells, NLOS, los, output_los, legacy_density, density_multiplier, num_pixels_redward);
	if (initial_only)
		printf("initial_only: synthesize the pre-quasar spectrum without evolving the gas\n");


    // age of the universe at redshift z in Myr
	double tz = (1.0 / H0) * (2.0 / 3.0) * pow(1.0 - omega_m, -0.5) * asinh(pow((1.0 - omega_m) / omega_m, 0.5) * pow(ai(z), 1.5)); 


    /* Set up the spatial grid */
	int ncells = Ncells;
	boxsize = Ncells * cellsize;
	if (boxsize_override > 0.0)
		boxsize = boxsize_override;
	/* Apply to both RT geometry and spectral integration. For a requested
	 * proper-length ratio P between input redshift zs and runtime z, use
	 * comoving_length_scale = P * (1+z)/(1+zs). Density is handled separately. */
	cellsize *= comoving_length_scale;
	Roffset *= comoving_length_scale;
	boxsize *= comoving_length_scale;
	printf("comoving_length_scale %.17g cellsize %.17g Roffset %.17g boxsize %.17g\n",
		comoving_length_scale, cellsize, Roffset, boxsize);
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
	char *File1, *File2, *File3, *File4, *File5, *File6, *File7, *File8, *File9;
	File1 = input_filename(input_dir, input_prefix, "nh");
	myFile = open_required_file(File1);
	File2 = input_filename(input_dir, input_prefix, "nh1");
	myFile2 = open_required_file(File2);

	File3 = input_filename(input_dir, input_prefix, "T");
	myFile3 = open_required_file(File3);

	File4 = input_filename(input_dir, input_prefix, "vpec");
	myFile4 = open_required_file(File4);
	File5 = input_filename(input_dir, input_prefix, "rho");
	myFile5 = open_required_file(File5);
	File6 = input_filename(input_dir, input_prefix, "Gammabg");
	myFile6 = open_required_file(File6);
	File7 = input_filename(input_dir, input_prefix, "Gammabghe1");
	myFile7 = open_required_file(File7);
	File8 = input_filename(input_dir, input_prefix, "Gammabghe2");
	myFile8 = open_required_file(File8);
	File9 = input_filename(input_dir, input_prefix, "Heatbg");
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
	static double nharr[RTE_MAX_CELLS][1000];	  // hydrogen number density array in cm^-3
	static double nh1arr[RTE_MAX_CELLS][1000];	  // neutral hydrogen number density array in cm^-3
	static double Tarr[RTE_MAX_CELLS][1000];		  // gas temperature array in K
	static double vpec[RTE_MAX_CELLS][1000];		  // peculiar velocity array in km/s
	static double rho[RTE_MAX_CELLS][1000];		  // gas density array in 2e40h^2kgkpc^-3
	static double gammabg[RTE_MAX_CELLS][1000];	  // photoionization rate of hydrogen in s^-1
	static double gammabghe1[RTE_MAX_CELLS][1000]; // photoionization rate of helium 1 in s^-1
	static double gammabghe2[RTE_MAX_CELLS][1000]; // photoionization rate of helium 2 in s^-1
	static double heatbg[RTE_MAX_CELLS][1000];	  // background heating rate

/* comment this out if you are not reading initial input from files */

	for (int i = 0; i < Ncells; ++i)
	{
		for (int j = 0; j < NLOS; ++j)
		{
			read_grid_value(myFile, File1, &nharr[i][j]);
			read_grid_value(myFile2, File2, &nh1arr[i][j]);
			read_grid_value(myFile3, File3, &Tarr[i][j]);
			read_grid_value(myFile4, File4, &vpec[i][j]);
			read_grid_value(myFile5, File5, &rho[i][j]);
			read_grid_value(myFile6, File6, &gammabg[i][j]);
			read_grid_value(myFile7, File7, &gammabghe1[i][j]);
			read_grid_value(myFile8, File8, &gammabghe2[i][j]);
			read_grid_value(myFile9, File9, &heatbg[i][j]);
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
		double density_scale = legacy_density ? 1.0 : 1.0 / X;
		double redshift_scale = legacy_density ? 1.0 : pow(1 + z_qso, 3) / pow(1 + z, 3);
		nh[i] = nharr[i][los] * density_scale * density_multiplier;							  // to set it to cosmological mean instead, use getnh(z) and convert to physcial coordinates;
		nh2[i][0] = (nharr[i][los] - nh1arr[i][los]) * density_scale * density_multiplier; // ionized hydrogen density in cm^-3

		// rescaling all densities to  zqso; ignore if the quasar is at the same redshift as the simulation (zqso=z)
		nh[i] = nh[i] * redshift_scale;
		nh2[i][0] = nh2[i][0] * redshift_scale;

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

	int encells = ncells + num_pixels_redward;

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


	if (initial_only)
	{
		storej = 0;
		for (int i = 0; i < ncells; ++i)
		{
			Tf[i] = T[i][0];
			NH1f[i] = rho[i][los] * (1.0 - nh2[i][0] / nh[i]) * (UnitMass_in_g / PROTONMASS / Hubble_h);
		}
	}

	for (int j = 0; !initial_only && j < max_global_steps; j++)
	{

		tsim = ti + Dt;
		tarr[j] = tsim;

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

		/*evolve all cells for a  global timestep*/
		evolve_eqns(nh, nh2, NH, nhe1, nhe2, nhe3, NHe1, NHe2, T, ncells, ti, tsim, dt, j, Dt); // my equivalent of subroutine rt_solve_cooling of Ramses RT (Rosdahl et al 2013)
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
    if (initial_only)
    {
        int length = snprintf(NULL, 0, "%s/noqso_z%0.2lfLOS%d.txt", output_dir, z_qso, output_los) + 1;
        fileSpec = malloc(length);
        if (fileSpec == NULL)
        {
            fprintf(stderr, "Failed to allocate initial spectrum filename\n");
            return 1;
        }
        snprintf(fileSpec, length, "%s/noqso_z%0.2lfLOS%d.txt", output_dir, z_qso, output_los);
    }
    else
        fileSpec = outfilename(output_dir, output_los, M1450,t_on,fduty, z_qso,tq); // output file name
    FILE *fptr = fopen(fileSpec, "w");
    save_outputs(fptr, ncells, nh, nh2, T, nhe1, nhe2, nhe3, Tau, los, M1450, z_qso,storej);    
    fclose(fptr);
    free(fileSpec);

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
