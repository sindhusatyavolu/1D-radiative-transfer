#include "../include/proximity_zone.h"
#include "../include/params.h"
#include "../include/constants.h"
#include <stdio.h>

/* Computing proximity zone size Rp */
double Rp(double *flux_arr, int lenarr)
{
	double rp = 0.0;
	for (int i = 0; i < lenarr; i++)
	{
		if (flux_arr[i] < 0.1)
		{
			rp = (1 + i) * cellsize / (1000 * (1.0 + z));   // in proper Mpc 
			break;
		}
	}
	return rp;
}
