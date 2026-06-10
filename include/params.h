#ifndef PARAMS_H
#define PARAMS_H

/* Define global variables */

extern double epsilon; // Courant factor; unused in this version
extern double Ng;		// ionising photon rate in Myr^-1

/* cosmology */
extern double X, Y;
extern double omega_m;
extern double H0;

/* quasar */
extern int NLOS, NFREQ;
extern double M1450;
extern double alpha, alpha2;
extern double fend;
extern double z;			   
extern double z_qso; // quasar redshift
extern double t_on;		   // quasar on time in Myr
extern double fduty;	   // quasar duty cycle
extern double tq;	   // quasar life time in Myr

/* simulation */
extern int Ncells;
extern double cellsize, boxsize, hubble_inv, Roffset;

/* IGM */
extern double background;

extern double Gammabkg[];
extern double Gammabkg_he1[];
extern double Gammabkg_he2[];
extern double Heat_bkg[];

void read_params_from_config(const char *filename);

#endif 



