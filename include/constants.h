#ifndef CONSTANTS_H
#define CONSTANTS_H


// physical constants 
#define PI 3.14159265358979323846
#define timefac 3.15e13				// multiply to convert Myr to sec or sec^-1 to Myr^-1
#define kb 1.38e-16 // Boltzmann constant in erg per kelvin
#define HUBBLE 3.2407789e-18  // in s^-1
#define BOLTZMANN 1.3806e-16  // erg K^-1
#define PROTONMASS 1.6726e-24 // g
#define GAMMA_LYA_H1 6.265e8  /* s^-1 */
#define LYMAN_ALPHA 1215.6e-8 /* 1215.6 Angstroem */
#define OSCILLATOR_STRENGTH 0.41615
#define THOMPSON 6.65245e-25 /* cm^2 */
#define Cv 2.9979e10		 /* cm/s */
#define kpc 3.08568025e21	 // in cm
#define nm2 1e-14			 // in cm^2
#define GRAVITY 6.67384e-8	 // cm^3 g^-1 s^-2


// functions
#define max(x, y) (((x) >= (y)) ? (x) : (y))
#define min(x, y) (((x) <= (y)) ? (x) : (y))
// The macro ai(z) is currently unused but retained for potential future use in scale factor calculations.
#define ai(z) 1.0 / (1.0 + z)


// simulation parameters
#define max_global_steps 12000
#define max_subcycling_steps 500000
#define T_MIN 0.1		// minimum temperature in K
#define T2_MIN_FIX 1e-2 // minimum temperature in K
#define X_MIN 1e-40		// minimum ionisation fraction
#define X_FM 1e-6
#define lim 2e-1
#define x_min_fix 1e-20
#define dt_min 1e-4			  // minimum time step in Myr
#define Tini 100.0	 // temperature slab in K

#endif