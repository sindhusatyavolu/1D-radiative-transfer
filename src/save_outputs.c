#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "../include/params.h"
#include "../include/constants.h"
#include "../include/utils.h"

char *outfilename(int los, float Mag, double t_on,double fduty, double z, double tq){
	char *fileSpec;
	int Length2;
	char *bet2;
	Length2 = snprintf(NULL, 0, "%0.2lf", Mag) + 1;
	bet2 = malloc(Length2);
	sprintf(bet2, "%0.2lf", Mag);

	int Length, String_len1;
	char *bet;
	Length = snprintf(NULL, 0, "%0.2lf", z) + 1;
	bet = malloc(Length);
	sprintf(bet, "%0.2lf", z);

	String_len1 = snprintf(NULL, 0, "z%sLOS%dtq%gmag%s_tep_%gMyr%gfd.txt",
						   bet, los, tq, bet2, t_on, fduty) + 1;
	fileSpec = malloc(String_len1);
	sprintf(fileSpec, "z%sLOS%dtq%gmag%s_tep_%gMyr%gfd.txt",
			bet, los, tq, bet2, t_on, fduty);

	free(bet);
	free(bet2);

    return fileSpec;
}

void save_outputs(FILE *fptr, int ncells, double *nh, double **nh2, double **T, double **nhe1, double **nhe2,
                    double **nhe3, double *Tau, int los, float Mag, double zqso, int storej)
{
    if (fptr == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

        for (int l = 0; l < ncells; ++l)
	{

		fprintf(fptr, "%lf %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e\n", ((double)l * cellsize + Roffset), (nh[l] - nh2[l][storej]) / nh[l], T[l][storej], nhe1[l][storej] / ((Y / (4.0 * ((1.0 - Y))) * nh[l])), nhe2[l][storej] / ((Y / (4.0 * ((1.0 - Y))) * nh[l])), nhe3[l][storej] / ((Y / (4.0 * ((1.0 - Y))) * nh[l])), Tau[l]);
	}

}
