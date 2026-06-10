#include "../include/params.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void load_files(){
/* Loading files containing densities, ionization, temperatures and background photoionization rates */
	FILE *myFile, *myFile2, *myFile3, *myFile4, *myFile5, *myFile6, *myFile7, *myFile8, *myFile9;
	char root[] = "f1kz";
	int length, string_len1, string_len2;
	char *Bet;
	char *File1, *File2, *File3, *File4, *File5, *File6, *File7, *File8, *File9;
	length = snprintf(NULL, 0, "%0.2lf", z) + 1;
	Bet = malloc(length);
	sprintf(Bet, "%0.2lf", z);
	string_len1 = strlen(root) + strlen("M3LOS_nh.txt") + strlen(Bet) + 3;
	File1 = malloc(string_len1);
	sprintf(File1, "%s%sM3LOS_nh.txt", root, Bet);
	myFile = fopen(File1, "r");
	string_len2 = strlen(root) + strlen("M3LOS_nh1.txt") + strlen(Bet) + 3;
	File2 = malloc(string_len2);
	sprintf(File2, "%s%sM3LOS_nh1.txt", root, Bet);
	myFile2 = fopen(File2, "r");

	string_len2 = strlen(root) + strlen("M3LOS_T.txt") + strlen(Bet) + 3;
	File3 = malloc(string_len2);
	sprintf(File3, "%s%sM3LOS_T.txt", root, Bet);
	myFile3 = fopen(File3, "r");

	string_len2 = strlen(root) + strlen("M3LOS_vpec.txt") + strlen(Bet) + 3;
	File4 = malloc(string_len2);
	sprintf(File4, "%s%sM3LOS_vpec.txt", root, Bet);
	myFile4 = fopen(File4, "r");
	string_len2 = strlen(root) + strlen("M3LOS_rho.txt") + strlen(Bet) + 3;
	File5 = malloc(string_len2);
	sprintf(File5, "%s%sM3LOS_rho.txt", root, Bet);
	myFile5 = fopen(File5, "r");
	string_len2 = strlen(root) + strlen("M3LOS_Gammabg.txt") + strlen(Bet) + 3;
	File6 = malloc(string_len2);
	sprintf(File6, "%s%sM3LOS_Gammabg.txt", root, Bet);
	myFile6 = fopen(File6, "r");
	string_len2 = strlen(root) + strlen("M3LOS_Gammabghe1.txt") + strlen(Bet) + 3;
	File7 = malloc(string_len2);
	sprintf(File7, "%s%sM3LOS_Gammabghe1.txt", root, Bet);
	myFile7 = fopen(File7, "r");
	string_len2 = strlen(root) + strlen("M3LOS_Gammabghe2.txt") + strlen(Bet) + 3;
	File8 = malloc(string_len2);
	sprintf(File8, "%s%sM3LOS_Gammabghe2.txt", root, Bet);
	myFile8 = fopen(File8, "r");
	string_len2 = strlen(root) + strlen("M3LOS_Heatbg.txt") + strlen(Bet) + 3;
	File9 = malloc(string_len2);
	sprintf(File9, "%s%sM3LOS_Heatbg.txt", root, Bet);
	myFile9 = fopen(File9, "r");

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
	double nharr[2048][500];	  // hydrogen number density array in cm^-3, 2048 is the number of cells, 500 is the number of lines of sight
	double nh1arr[2048][500];	  // neutral hydrogen number density array in cm^-3
	double Tarr[2048][500];		  //  gas temperature array in K
	double vpec[2048][500];		  //	peculiar velocity array in km/s
	double rho[2048][500];		  // gas density array in 2e40h^2kgkpc^-3
	double gammabg[2048][500];	  // photoionization rate of hydrogen in Myr^-1
	double gammabghe1[2048][500]; // photoionization rate of helium 1 in Myr^-1
	double gammabghe2[2048][500]; // photoionization rate of helium 2 in Myr^-1
	double heatbg[2048][500];	  // background heating rate in Myr^-1


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
}
