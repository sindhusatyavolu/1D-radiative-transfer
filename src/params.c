#include "../include/params.h"

double epsilon;
double Ng;

/* cosmology */
double X = 0.76, Y = 0.24;
double omega_m = 0.308;
double H0 = 67.8 * 1.05e-6;

/* quasar */
int NLOS = 500, NFREQ = 80;
double z = 5.95;
double M1450 = -27.0;
double alpha = -1.7, alpha2 = -0.61;
double fend = 40.0;
double t_on = 1.0;
double fduty = 1.0;
double tq = 0.0;
double z_qso = 5.95;

/* simulation */
int Ncells = 2048;
double cellsize = 114.88970588235293, boxsize = 235294.118, hubble_inv = 1.0, Roffset = 114.88970588235293;

/* IGM */
double background = 1.0;
double Gammabkg[2048];
double Gammabkg_he1[2048];
double Gammabkg_he2[2048];
double Heat_bkg[2048];
