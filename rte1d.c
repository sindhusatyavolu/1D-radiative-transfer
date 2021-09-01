#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

/*****************************************************************************/

/*
  1D radiative transfer code that outputs Lyman alpha optical depth along a given sightline
  
*/

/*****************************************************************************/
#define count 10000000
#define countl 1000000000
#define max(x,y) (((x) >= (y)) ? (x) : (y))
#define min(x,y) (((x) <= (y)) ? (x) : (y))

/*Physical parameters*/

#define cellsize 0.01  // in *comoving kpc
#define boxsize 6.6 // in *comoving kpc
#define hubble_inv 1.0 // 1/h
#define Roffset 0.01 // in *comoving kpc
#define INITIAL_TEMP 1e4
#define INITIAL_HII 1.2e-3
#define INITIAL_HEII 0.0
#define INITIAL_HEIII 0.0
#define Ng 3.15e13*5e48// timefac*5e53 // in Myr^-1
#define alpha 0.0
#define X  1.0
#define Y 0.0
#define epsilon 5e-1
#define z 0.0 //end redshift
#define background 0.0
#define MONOFREQ 0.0
#define ISOTHERMAL 0.0
#define NOHelium 0.0
#define NFREQ 80
#define fend 40 // freq sampled from nu_0 to fend*nu_0
#define COMOVING 0.0

/* Numerical limits */

#define T_MIN 0.1
#define T2_MIN_FIX 1e-2
#define X_MIN 1e-40
#define X_FM 1e-6
#define lim 2e-1 // max. change allowed bw timesteps
#define x_min_fix 1e-20
#define dt_min 1e-8

/* Physical Constants */

#define PI 3.14159265358979323846
#define timefac 3.15e13 //multiply to convert Myr to sec or sec^-1 to Myr^-1
#define kb 1.38e-16 // erg per kelvin
#define HUBBLE 3.2407789e-18
#define  BOLTZMANN 1.3806e-16
#define  PROTONMASS  1.6726e-24
#define  GAMMA_LYA_H1   6.265e8  /* s^-1 */
#define  LYMAN_ALPHA      1215.6e-8      /* 1215.6 Angstroem */
#define  OSCILLATOR_STRENGTH       0.41615
#define  THOMPSON     6.65245e-25
#define  Cv           2.9979e10
#define kpc 3.08568025e21
#define nm2 1e-14
#define GRAVITY  6.67384e-8

/* Cumilative sum */
double cumsum(double *nH, double **nH2, int i, int j){
  double sum =0;
  for(int l =0;l<i;l++){
  sum = sum + nH[l] -nH2[l][j];
}
return sum;
}

double hecumsum(double **nHe, int i, int j){
    double sum =0;
    for(int l =0;l<i;l++){
    sum = sum + nHe[l][j];
}
return sum;
}

double Hubble(double t){
    double  omega_m = 0.308;
    double omega_l = 1-omega_m;
    double  H0 = 67.8*1.05e-6;
   // double z =  0;//5.945518;
    double a = pow(sqrt(omega_m/omega_l)*sinh(1.5*H0*sqrt(omega_l)*t), 2.0/3.0);//pow(1+z,-3);
    double zz = 1.0/a -1.0;
   // return H0*sqrt(omega_m/pow(a,3) +1-omega_m);
    return sqrt(omega_m)*H0*cosh(1.5*H0*sqrt(1-omega_m)*t)*pow(1+zz,1.5);///sinh(1.5*H0*sqrt(1-omega_m)*t); //sqrt(omega_m/a3 +1-omega_m);
}

double ascale(double t){
    double  omega_m = 0.308;
    double omega_l = 1-omega_m;
    double  H0 = 67.8*1.05e-6;
    if (COMOVING ==0.0){
    return 1.0/(1.0+z);
    }
    else
        return pow(sqrt(omega_m/omega_l)*sinh(1.5*H0*sqrt(omega_l)*t), 2.0/3.0);
}


double getnh(double redshift){
    double omegam = 0.308;
    double omegab = 0.0482;
    double omegal = 0.692;
    double h = 0.678;
    double H0 = 100*h;  //km s^-1 Mpc^-1
    double G = 6.67e-11; //Nm^2/kg^2 == m^3s^-2kg^-1?
    double kmtoMPc = 1.0/3.08e19;
    double rhocr = pow((kmtoMPc),2)*3.0*pow(H0,2)/(8.0*PI*G) ;   //km^2  mpc^-2/m^3kg^-1 = kg/m^3 *(km/Mpc)^2
    double Xfac = 1.0/X ;//#(1.0+4.0*(1.0/X -1.0));
    //omegah^2 = rho,0/rhocr
    //rho,0 = rho,z*(1+z)**3
    double mH = 1.67e-27; //kg
    double rhoz =  rhocr*omegab*pow((1.0+z),3);
    return 1e-6*rhoz/(mH*Xfac); //in cm^-3
}

/*electron collisonal ionisation rates in Myr^-1 cm^3**/
double GE(double Ti){
    double Th1 = 157807;
    double lam = 2*Th1/Ti;
    double T5 =Ti/100000 ;
   // return timefac*5.85e-11*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-157809.1/Ti);
   //return timefac*1.17e-10*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-157809.1/Ti);//bolton
    return  timefac*21.11*pow(Ti,-1.5)*exp(-lam/2)*pow(lam,-1.089)/pow(1.0+pow((lam/0.354),0.874),1.101) ; // HnG 97
}

double betahe1(double Ti){
    double lam  = 2*285335/Ti;
    return timefac*32.38*pow(Ti,-1.5)*exp(-lam/2)*pow(lam,-1.146)/pow(1.0+pow((lam/0.416),0.987),1.056) ; // from HnG
   // double T5 =Ti/100000 ;
   // return timefac*2.38e-11*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-285335.4/Ti); //ramses
   // return timefac*4.76e-11*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-285335.4/Ti);//bolton
}

double betahe2(double Ti){
    double lam = 2*631515/Ti;
    return timefac*19.95*pow(Ti,-1.5)*exp(-lam/2)*pow(lam,-1.089)/pow(1.0+pow((lam/0.553),0.735),1.275) ; //HnG 97
    //double T5 =Ti/100000 ;
    //return timefac*5.68e-12*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-631515.0/Ti); //ramses
    //return timefac*1.14e-11*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-631515.0/Ti); //bolton
}

/*Case A Recombination rates in cm^3 Myr^-1*/

double ralpha(double Ti){
    double Th1 = 157807;
    double lam = 2*Th1/Ti;
    return timefac*1.269e-13*pow(lam,1.503)/pow((1.0+pow(lam/0.522,0.47)),1.923); //  HnG 97
  //  double Tl = kb*Ti/1.6e-12;
  //  return timefac*(exp(-28.6130338-0.72411256*log(Tl)-2.02604473e-2*log(pow(Tl,2))-2.38086188e-3*log(pow(Tl,3))-3.21260521e-4*log(pow(Tl,4))-1.42150291e-5*log(pow(Tl,5))+4.98910892e-6*log(pow(Tl,6))+5.75561414e-7*log(pow(Tl,7))-1.85676704e-8*log(pow(Tl,8))-3.07113524e-9*log(pow(Tl,9))));
}

double alphahe2(double Ti){
    double lam = 570670/Ti;
    double Tl = kb*Ti/1.6e-12;
    return timefac*3e-14*pow(lam,0.654); // HnG 97
//   return timefac*(3.925e-13*pow(Tl,-0.6353)+1.544e-9*pow(Tl,-1.5)*exp(-48.596/Tl)*(0.3+exp(8.1/Tl)));
}

double alphahe3(double Ti){
    double lam = 1263030/Ti;
    return timefac*2.538e-13*pow(lam,1.503)/pow((1.0+pow(lam/0.522,0.47)),1.923); //HnG 97
  //  return 2*ralpha(Ti/4);//bolton
}


/*Photoionisation cross sections in nm^2*/

double sigma(double f){
    double sigma_0=5.475 ; // in nm^2
    double nu_0 = 0.4298 ;
    double y0 = 0;
    double x =  f/nu_0 - y0;
    double yw = 0;
    double P =  2.963;
    double ya = 32.88 ;
    double y1  = 0 ;
    double y = sqrt(x*x+y1*y1);
   return sigma_0*(pow((x-1),2)+pow(yw,2))*pow(y,0.5*P-5.5)/pow(1+sqrt(y/ya),P);
  // return 6.3e-18*1e14*(1.34*pow(f/13.6,-2.99)-0.34*pow(f/13.6,-3.99)); //bolton cross sections
}

double sigmahe1(double f){
    double sigma_0= 9.492e-2 ; // in nm^2
    double nu_0 = 13.61;
    double y0 = 0.4434;
    double x =  f/nu_0 - y0;
    double yw = 2.039;
    double P =  3.188;
    double ya = 1.469 ;
    double y1  = 2.136 ;
    double y = sqrt(x*x+y1*y1);
    return sigma_0*(pow((x-1),2)+pow(yw,2))*pow(y,0.5*P-5.5)/pow(1+sqrt(y/ya),P);
  //return 7.03e-18*1e14*(1.66*pow(f/24.6,-2.05)-0.66*pow(f/24.6,-3.05)); //bolton cross sections
}

double sigmahe2(double f){
    double sigma_0 = 1.369 ; // in nm^2
    double nu_0 = 1.720 ;
    double y0 = 0;
    double x =  f/nu_0 - y0;
    double yw = 0;
    double P =  2.963;
    double ya = 32.88 ;
    double y1  = 0 ;
    double y = sqrt(x*x+y1*y1);
    return sigma_0*(pow((x-1),2)+pow(yw,2))*pow(y,0.5*P-5.5)/pow(1+sqrt(y/ya),P);
  //  return 1.5e-18*1e14*(1.34*pow(f/54.5,-2.99)-0.34*pow(f/54.5,-3.99)); //bolton
}

double PH1(double x,double nh1,double nhe1l,double nhe2l,double t){
    double qH = exp(-sigma(x)*nh1*cellsize*hubble_inv*kpc*nm2*ascale(t));
    double qHe1 = exp(-sigmahe1(x)*nhe1l*cellsize*hubble_inv*kpc*nm2*ascale(t));
    double  qHe2 = exp(-sigmahe2(x)*nhe2l*cellsize*hubble_inv*kpc*nm2*ascale(t)) ;
    double  tau_tot = sigma(x)*nh1*cellsize*hubble_inv*kpc*ascale(t)*nm2+sigmahe1(x)*nhe1l*cellsize*hubble_inv*kpc*nm2*ascale(t)+ sigmahe2(x)*nhe2l*cellsize*hubble_inv*kpc*nm2*ascale(t);
    double  D = (1-qH)*qHe2*qHe1 +(1-qHe1)*qH*qHe2 + (1-qHe2)*qH*qHe1;
    //printf("%0.12e\n",D);
    if (D ==0 || D!= D)
    {
        printf("D %0.12e\n",D);
        printf("%0.12e %0.12e\n",nhe1l,nhe2l );
        exit(0);
    }
  return (1-qH)*qHe2*qHe1*(1-exp(-tau_tot))/D;
}

double PHe1(double x,double nh1,double nhe1l,double nhe2l,double t){
    double qH = exp(-sigma(x)*nh1*cellsize*hubble_inv*kpc*nm2*ascale(t));
    double qHe1 = exp(-sigmahe1(x)*nhe1l*cellsize*hubble_inv*kpc*nm2*ascale(t));
    double  qHe2 = exp(-sigmahe2(x)*nhe2l*cellsize*hubble_inv*kpc*nm2*ascale(t))  ;
    double  tau_tot = sigma(x)*nh1*cellsize*hubble_inv*kpc*nm2*ascale(t)+sigmahe1(x)*nhe1l*cellsize*hubble_inv*kpc*nm2*ascale(t)+ sigmahe2(x)*nhe2l*cellsize*hubble_inv*kpc*nm2*ascale(t);
    double  D = (1-qH)*qHe2*qHe1 +(1-qHe1)*qH*qHe2 + (1-qHe2)*qH*qHe1;
    return (qH)*qHe2*(1-qHe1)*(1-exp(-tau_tot))/D;
}

double PHe2(double x,double nh1,double nhe1l,double nhe2l,double t){
    double qH = exp(-sigma(x)*nh1*cellsize*hubble_inv*kpc*nm2*ascale(t));
    double qHe1 = exp(-sigmahe1(x)*nhe1l*cellsize*hubble_inv*kpc*nm2*ascale(t));
    double  qHe2 = exp(-sigmahe2(x)*nhe2l*cellsize*hubble_inv*kpc*nm2*ascale(t))  ;
    double  tau_tot = sigma(x)*nh1*cellsize*hubble_inv*kpc*nm2*ascale(t)+sigmahe1(x)*nhe1l*cellsize*hubble_inv*kpc*nm2*ascale(t)+ sigmahe2(x)*nhe2l*cellsize*hubble_inv*kpc*nm2*ascale(t);
    double  D = (1-qH)*qHe2*qHe1 +(1-qHe1)*qH*qHe2 + (1-qHe2)*qH*qHe1;
   
    return (qH)*(1-qHe2)*qHe1*(1-exp(-tau_tot))/D;
}



/*Photoionisation rates in energy units eV**/

double Gamma(double NH,double Nhe1, double Nhe2,double nh1,double nhe1l,double nhe2l,double t){
    // SPACING = (b/a)^(1/n) , n= 40, b = 136, a =13.6
// double u[]={13.6,14.2501077683,14.9312920154,15.6450382603,16.392903033,17.1765172688,17.9975898651,18.8579114079,19.7593580771,20.7038957377,21.6935842269,22.7305818467,23.81715007,24.955658473,26.1485899022,27.3985458894,28.708252325,30.0805654024,31.5184778468,33.0251254419,34.6037938684,36.257925869,37.9911287566,39.807182281,41.710046872,43.7038722756,45.7930066044,47.9820058198,50.2756436672,52.6789220867,55.1970821217,57.835615348,60.6002758498,63.4970927684,66.5323834504,69.7127672245,73.0451798366,76.5368885757,80.1955081218,84.0290171524,88.0457757418,92.2545435932,96.6644991413,101.285259569,106.126901779,111.199984373,116.515570673,122.085252852,127.92117721,134.036070669,140.443268521,147.156743511,154.191136314,161.561787456,169.284770772,177.376928458,185.855907805,194.740199677,204.04917884,213.803146214,224.023373144,234.732147789,245.952823726,257.709870885,270.028928902,282.936863046,296.461822797,310.633303239,325.482209382,341.040923556,357.343376036,374.425119031,392.323404218,411.077263981,430.727596533,451.317255102,472.891141391,495.496303493,519.182038499,544.0};  /* 80 bins nu to 40*nu */
    double a = 13.6;
    double b = 13.6*fend;
    double u[NFREQ];
    //double *u = (double*) malloc(NFREQ*sizeof(double));
    for (int i=0; i<NFREQ; i++) {
        u[i] = a*pow(b/a,i/(NFREQ-1.0));
        //printf("%0.12e\n",u[i]);
    }
  
    double fsum =0;
    double x0 = u[0];
   
    int len = sizeof(u)/sizeof(u[0]);
  
 /*Integral for calculating photoionisation rate*/
    for(int i=1;i<len;i++){
   
     double df = (u[i]-u[i-1])*0.5;
   
     fsum = fsum + df*(pow(u[i]/x0,alpha-1)*exp(-sigma(u[i])*NH-sigmahe1(u[i])*Nhe1-sigmahe2(u[i])*Nhe2)*PH1(u[i],nh1,nhe1l,nhe2l,t)+pow(u[i-1]/x0,alpha-1)*exp(-sigma(u[i-1])*NH-sigmahe1(u[i-1])*Nhe1-sigmahe2(u[i-1])*Nhe2)*PH1(u[i-1],nh1,nhe1l,nhe2l,t));
}


#ifdef MONOFREQ
    fsum  = exp(-sigma(u[0])*NH-sigmahe1(u[0])*Nhe1-sigmahe2(u[0])*Nhe2)*PH1(u[0],nh1,nhe1l,nhe2l,t) ;
#endif
 
  /*isothermal case*/
   //  fsum  = df*exp(-sigma(u[0])*NH-sigmahe1(u[0])*Nhe1-sigmahe2(u[0])*Nhe2)*PH1(u[0],nh1,nhe1l,nhe2l) ;
    
    return fsum;
}

double Gammahe1(double NH,double Nhe1, double Nhe2,double nh1,double nhe1l,double nhe2l,double t){
    
   // double u[]={24.6,25.775930228,27.0080723219,28.2991133237,29.6518687214,31.0692885892,32.5544640207,34.1106338703,35.741191816,37.4496937608,39.2398655869,41.1156112815,43.0810214502,45.1403822379,47.298184676,49.5591344765,51.9281622938,54.4104344778,57.0113643405,59.7366239612,62.5921565561,65.5841894395,68.7192476038,72.0041679495,75.4461141949,79.0525924985,82.8314678286,86.7909811153,90.9397672215,95.2868737745,99.8417808967,104.614421879,109.615204846,114.855035449,120.345340653,126.098093656,132.125839999,138.441724924,145.059522044,151.993663379,159.259270827,166.872189147,174.849020506,183.207160691,191.964837042,201.141148205,210.756105777,220.830677952,231.386835248,242.447598416,254.037088647,266.180580175,278.904555392,292.236762604,306.206276543,320.84356177,336.180539118,352.250655298,369.088955842,386.732161534,405.218748481,424.589032029,444.885254682,466.151678218,488.43468022,511.782855215,536.247120647,561.880827918,588.739878735,616.882847021,646.371106654,677.268965306,709.643804688,743.566227495,779.110211375,816.353270258,855.376623398,896.265372495,939.108687286,984.0};
 
    double a = 24.6;
    double b = 24.6*fend;
    
    double u[NFREQ];
    for (int i=0; i<NFREQ; i++) {
        u[i] = a*pow(b/a,i/(NFREQ-1.0));
        //printf("%0.12e\n",u[i]);
    }
    
    double fsum=0 ;
    double x0 = 13.6 ;//u[0];
    
    int len = sizeof(u)/sizeof(u[0]);

 
for(int i=1;i<len;i++){

     double df = (u[i]-u[i-1])*0.5;
  
     fsum = fsum + df*(pow(u[i]/x0,alpha-1)*exp(-sigma(u[i])*NH-sigmahe1(u[i])*Nhe1-sigmahe2(u[i])*Nhe2)*PHe1(u[i],nh1,nhe1l,nhe2l,t)+pow(u[i-1]/x0,alpha-1)*exp(-sigma(u[i-1])*NH-sigmahe1(u[i-1])*Nhe1-sigmahe2(u[i-1])*Nhe2)*PHe1(u[i-1],nh1,nhe1l,nhe2l,t));
}

#ifdef MONOFREQ
    fsum  = exp(-sigma(u[0])*NH-sigmahe1(u[0])*Nhe1-sigmahe2(u[0])*Nhe2)*PHe1(u[0],nh1,nhe1l,nhe2l,t) ;
#endif
 
  /*isothermal case*/
   //  fsum  = df*exp(-sigma(u[0])*NH)*PH1(u[0],nh1,nhe1,nhe2) ;
   
 //   fsum  = df*exp(-sigma(u[0])*NH-sigmahe1(u[0])*Nhe1-sigmahe2(u[0])*Nhe2)*PHe1(u[0],nh1,nhe1,nhe2) ;
    return fsum;
}

double Gammahe2(double NH,double Nhe1, double Nhe2,double nh1,double nhe1l,double nhe2l,double t){
    

   // double u[]={54.5,57.1052112775,59.8349569734,62.6951900871,65.6921481837,68.8323669964,72.12269468,75.5703067451,79.1827217061,82.9678174781,86.9338485564,91.089464018,95.4437263836,100.006131381,104.786628652,109.795643454,115.044099391,120.543442237,126.305664901,132.343333573,138.669615134,145.298305872,152.243861561,159.521428994,167.146879009,175.136841104,183.508739702,192.280832146,201.472248519,211.103033362,221.194189385,231.76772327,242.846693663,254.455261462,266.618742504,279.363662775,292.717816257,306.710325542,321.371705341,336.73392903,352.830498377,369.696516605,387.368764941,405.885782831,425.287951984,445.617584437,466.91901483,489.238697089,512.625305733,537.12984202,562.805745174,589.709008924,617.898303613,647.435104144,678.383824048,710.811955954,744.790218777,780.39271194,817.697076968,856.784666812,897.74072326,940.654562829,985.619771551,1032.73440906,1082.10122244,1133.82787029,1188.02715753,1244.81728136,1304.32208907,1366.67134808,1432.00102897,1500.453602,1572.17834778,1647.33168287,1726.07750081,1808.58752964,1895.04170631,1985.62856915,2080.54566899,2180.0};
    double a =54.5;
    double b = 54.5*fend;
    double u[NFREQ];
    for (int i=0; i<NFREQ; i++) {
        u[i] = a*pow(b/a,i/(NFREQ-1.0));
    }
     
    double fsum=0 ;
    double x0 = 13.6;//u[0];
    
    int len = NFREQ;//sizeof(u)/sizeof(u[0]);
 
 
for(int i=1;i<len;i++){
    
     double df = (u[i]-u[i-1])*0.5;
   
     fsum = fsum + df*(pow(u[i]/x0,alpha-1)*exp(-sigma(u[i])*NH-sigmahe1(u[i])*Nhe1-sigmahe2(u[i])*Nhe2)*PHe2(u[i],nh1,nhe1l,nhe2l,t)+pow(u[i-1]/x0,alpha-1)*exp(-sigma(u[i-1])*NH-sigmahe1(u[i-1])*Nhe1-sigmahe2(u[i-1])*Nhe2)*PHe2(u[i-1],nh1,nhe1l,nhe2l,t));
}

#ifdef MONOFREQ
    fsum  = exp(-sigma(u[0])*NH-sigmahe1(u[0])*Nhe1-sigmahe2(u[0])*Nhe2)*PHe2(u[0],nh1,nhe1l,nhe2l,t) ;
#endif

  /*isothermal case*/
 //    fsum  = df*exp(-sigma(u[0])*NH-sigmahe1(u[0])*Nhe1-sigmahe2(u[0])*Nhe2)*PHe2(u[0],nh1,nhe1,nhe2) ;
   

    return fsum;
}

/*Photoheating rates - ev^2*/
double Heat(double NH,double Nhe1,double Nhe2,double nh1,double nhe1l,double nhe2l,double t){
    
   
    double a = 13.6;
    double b = 13.6*fend;
    double u[NFREQ];
   
    for (int i=0; i<NFREQ; i++) {
        u[i] = a*pow(b/a,i/(NFREQ-1.0));
        //printf("%0.12e\n",u[i]);
    }
    
    double fsum = 0 ;
    double x0 = u[0];
    
    int len = sizeof(u)/sizeof(u[0]);
    double df = 1;//0.5*pow(u[len-1]/x0,1.0/(len-1));
 
    for(int i=1;i<len;i++){
     df = (u[i]-u[i-1])*0.5;
     fsum = fsum + df*((u[i]-x0)*pow(u[i]/x0,alpha-1)*exp(-sigma(u[i])*NH-sigmahe1(u[i])*Nhe1-sigmahe2(u[i])*Nhe2)*PH1(u[i],nh1,nhe1l,nhe2l,t)+(u[i-1]-x0)*pow(u[i-1]/x0,alpha-1)*exp(-sigma(u[i-1])*NH-sigmahe1(u[i-1])*Nhe1-sigmahe2(u[i-1])*Nhe2)*PH1(u[i-1],nh1,nhe1l,nhe2l,t)) ;
    // printf("%0.12e\n",PH1(u[i],nh1,nhe1l,nhe2l));
}

        
    return fsum;
}
double Heat2(double NH,double Nhe1, double Nhe2,double nh1,double nhe1l,double nhe2l,double t){
    
 
    double a = 24.6;
    double b = 24.6*fend;
    
    double u[NFREQ];
    for (int i=0; i<NFREQ; i++) {
        u[i] = a*pow(b/a,i/(NFREQ-1.0));
        //printf("%0.12e\n",u[i]);
    }
    
    double fsum = 0 ;
    double x0 = u[0];
    
    int len = sizeof(u)/sizeof(u[0]);
    double df = 1;//0.5*pow(u[len-1]/x0,1.0/(len-1));

    for(int i=1;i<len;i++){
      df = (u[i]-u[i-1])*0.5;
     fsum = fsum + df*((u[i]-x0)*pow(u[i]/13.6,alpha-1)*exp(-sigma(u[i])*NH-sigmahe1(u[i])*Nhe1-sigmahe2(u[i])*Nhe2)*PHe1(u[i],nh1,nhe1l,nhe2l,t)+(u[i-1]-x0)*pow(u[i-1]/13.6,alpha-1)*exp(-sigma(u[i-1])*NH-sigmahe1(u[i-1])*Nhe1-sigmahe2(u[i-1])*Nhe2)*PHe1(u[i-1],nh1,nhe1l,nhe2l,t)) ;
     
}
 
  
    return fsum;
}

double Heat3(double NH,double Nhe1, double Nhe2,double nh1,double nhe1l,double nhe2l,double t){
  
    double a = 54.5;
    double b = 54.5*fend;
    
    double u[NFREQ];
    for (int i=0; i<NFREQ; i++) {
        u[i] = a*pow(b/a,i/(NFREQ-1.0));
        //printf("%0.12e\n",u[i]);
    }
    
    double fsum=0 ;
    double x0 = u[0];
    
    int len = sizeof(u)/sizeof(u[0]);
    double df = 1;//0.5*pow(u[len-1]/x0,1.0/(len-1));

    for(int i=1;i<len;i++){
      df = (u[i] - u[i-1])*0.5;
     fsum = fsum + df*((u[i]-x0)*pow(u[i]/13.6,alpha-1)*exp(-sigma(u[i])*NH-sigmahe1(u[i])*Nhe1-sigmahe2(u[i])*Nhe2)*PHe2(u[i],nh1,nhe1l,nhe2l,t)+(u[i-1]-x0)*pow(u[i-1]/13.6,alpha-1)*exp(-sigma(u[i-1])*NH-sigmahe1(u[i-1])*Nhe1-sigmahe2(u[i-1])*Nhe2)*PHe2(u[i-1],nh1,nhe1l,nhe2l,t)) ;
     
}
 
    return fsum;
}

double xh2solve(double Ti,double xh2,double xhe1,double xhe2,double xhe3,double dti,double NH1,double Nhe1, double Nhe2,double R,double nhi, double *xrate,double t){
    
    double nhe1l = max((Y/(4.0*((1.0-Y))))*nhi*(x_min_fix),(Y/(4.0*((1.0-Y))))*nhi*(xhe1));
    double nhe2l = max((Y/(4.0*((1.0-Y))))*nhi*xhe2,(Y/(4.0*((1.0-Y))))*nhi*x_min_fix);
    double nhe3l = max((Y/(4.0*((1.0-Y))))*nhi*xhe3,(Y/(4.0*((1.0-Y))))*nhi*x_min_fix);
   
    double ne = xh2*nhi+(Y/(4.0*((1.0-Y))))*nhi*(xhe2+2*xhe3);
    double xe = xh2+(xhe2+2*xhe3);
    double nu = (1.0 + 4.0 * Y) / (1.0 + Y + xe); //bolton
   // double Ti = (2.0/3.0)*Tu*(nu*PROTONMASS)/(kb);
    
    double fac = -alpha*Ng/(pow(3e21,3)*13.6);  // Myr^-1 eV^-1
    /*isothermal case*/
    
#ifdef ISOTHERMAL
     fac = Ng/pow(3e21,3);
#endif
     
    double vol = 1.0/(4*PI*R*R*cellsize*hubble_inv*ascale(t));  // kpc^-3 converted to cm^-3 through fac
    double nh1 = max(nhi*x_min_fix,nhi*(1-xh2));
    double C,D,dbdt,dalphadt,Tnu,dCdx,dDdx,J;
   
    double Gammab;
    if (background==1)
    {
      Gammab = 0.23e-12*timefac; // from haardt madau 2012
    }
    else
      Gammab = 0;//0.23e-12*timefac; // from haardt madau 2012
    
   // Gammab = timefac*1.27e-11*Jb*pow(-alpha+3,-1);
    
    C = GE(Ti)*ne + fac*vol*Gamma(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l,t)/(nh1) + Gammab;
    D = ralpha(Ti)*ne;
    if (C !=C)
    {
        printf("H1 gone wrong\n");
        printf("gamma %0.12e %0.12e %0.12e\n", fac*vol, Gamma(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l,t),1/(nh1));
        printf("elec coll %0.12e %0.12e\n", GE(Ti),ne);
        exit(0);
    }
 
  /*
    dbdt = (1.07732*1e14*exp(-157807/Ti))/(pow((1 + 158648*pow((1/Ti),0.874)),1.101) *pow((1/Ti),1.089) *pow(Ti,(7/2)))
- (1.02402*1e9*exp(-157807/Ti))/(pow((1 + 158648*pow((1/Ti),0.874)),1.101)*pow((1/Ti),1.089)*pow(Ti,(5/2)))
+ (7.43439*1e8*exp(-157807/Ti))/(pow((1 + 158648*pow((1/Ti),0.874)),1.101)*pow((1/Ti),0.089)*pow(Ti,(3/2)))
+ (1.0422*1e14*exp(-157807/Ti)*pow((1/Ti),0.785))/(pow((1 + 158648*pow((1/Ti),0.874)),2.101)*pow(Ti,(3/2)));
    */
    dbdt = -2.91364*exp(-157809.1/Ti)/pow((1+0.00316228*pow(Ti,0.5)),2) +2.90803*exp(-157809.1/Ti)/((1+0.00316228*pow(Ti,0.5))*pow(Ti,1.5)) + 921.375*exp(-157809.1/Ti)/((1+0.00316228*pow(Ti,0.5))*pow(Ti,0.5));
    
    dalphadt = -((1.10726*1e9*pow((1/Ti),2.503))/(1 + 521.553*pow((1/Ti),0.47))*1.923) + (
 3.4727*1e11*pow((1/Ti),2.973))/pow((1 + 521.553*pow((1/Ti),0.47)),2.923);

    nu = (1.0 + 4.0 * Y) / (1.0 + Y + ne); // bolton
    //nu = pow((X*(1+xh2)+Y*(1+xhe2+2*xhe3)/4),(-1));
    Tnu = Ti/nu;
    
    dCdx = nhi*GE(Ti) + ne*Tnu*pow(nu,2)*X*dbdt ;
    dDdx = nhi*ralpha(Ti) + ne*Tnu*pow(nu,2)*X*dalphadt;
    J = dCdx -(C+D) - xh2*(dCdx+dDdx);
    *(xrate) = (double) max(0.0,(C-xh2*(C+D)));
   //return (xh2 + C*dti)/(1+D*dti);
    if ((1-J*dti)==0 || *(xrate)<0)
    {
        printf("x denom is zero\n");
        printf("%0.12e %0.12e %0.12e %0.12e %0.12e negative dervative\n",(C-xh2*(C+D)),xh2,C,D,xh2*(C+D) );
        exit(0);
    }
   return xh2 + dti*(*(xrate))/(1-J*dti);
  //  return  (xh2 + C*dti)/(1+D*dti);
}

double xhe1solve(double Ti,double xh2i,double xhe1,double xhe2,double xhe3,double dti,double NH1,double Nhe1,double Nhe2,double R,double nhi,double t){
        double nhe1l = max((Y/(4.0*((1.0-Y))))*nhi*(x_min_fix),(Y/(4.0*((1.0-Y))))*nhi*(xhe1));
        double nhe2l = max((Y/(4.0*((1.0-Y))))*nhi*xhe2,(Y/(4.0*((1.0-Y))))*nhi*x_min_fix);//(Y/(4*((1.0-Y))))*nhi*xhe2;
        double nhe3l = max((Y/(4.0*((1.0-Y))))*nhi*xhe3,(Y/(4.0*((1.0-Y))))*nhi*x_min_fix);//(Y/(4*((1.0-Y))))*nhi*xhe3;
   // double xhe1 = 1- xhe2 -xhe3;

        double ne = xh2i*nhi+(Y/(4.0*((1.0-Y))))*nhi*(xhe2+2*xhe3);
    double xe = xh2i+(xhe2+2*xhe3);
    double nu = (1.0 + 4.0 * Y) / (1.0 + Y + xe); //bolton
   // double Ti = (2.0/3.0)*Tu*(nu*PROTONMASS)/(kb);
        double fac = -alpha*Ng/(pow(3e21,3)*13.6);

#ifdef ISOTHERMAL
     fac = Ng/pow(3e21,3);
#endif
    
    /*isothermal case */
    //double fac = Ng/pow(3e21,3);
        double vol = 1.0/(4*PI*R*R*cellsize*hubble_inv*ascale(t));
    
        double nh1 = max(nhi*(1-xh2i),nhi*x_min_fix);//nhi*(1-xh2i);
        double C,D,Tnu;
       //nu = pow((X*(1+xh2i)+Y*(1+xhe2+2*xhe3)/4),(-1));
     //  nu = (1.0 + 4.0 * Y ) / (1.0 + Y + ne); //bolton
      // Tnu = Ti/nu;
        double Gammab;
        if (background==1)
         {
           Gammab = 0.154e-12*timefac; // from haardt madau 2012
         }
         else
          Gammab = 0;//0.154e-12*timefac; // from haardt madau 2012


        C = xhe2*alphahe2(Ti)*ne;//+ xhe3*alphahe3(Ti)*ne;
   // D = betahe2(Ti)*ne + alphahe2(Ti)*ne + fac*vol*Gammahe2(NH1,Nhe1,Nhe2,nh1,nhe1,nhe2)/(nhe2) ;
        D = betahe1(Ti)*ne + fac*vol*Gammahe1(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l,t)/nhe1l + Gammab;
        if (D !=D)
    {
        printf("He1 gone wrong\n");
        printf("gamma %0.12e %0.12e %0.12e\n", fac*vol, Gammahe1(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l,t),1/(nhe1l));
        printf("elec coll %0.12e,%0.12e, %0.12e\n",xhe2, betahe1(Ti),ne);
        exit(0);
    }
 
        if ((nhe1l+D*dti)==0)
    {
        printf("xhe1 denom is zero\n");
        exit(0);
    }
   return (xhe1 + C*dti)/(1+D*dti);
   //     return xhe1 + (C-D*xhe1)*dti;
}


double xhe2solve(double Ti,double xh2i,double xhe1,double xhe2,double xhe3,double dti,double NH1,double Nhe1,double Nhe2,double R,double nhi,double *hrate,double t){
    
    double nhe1l = max((Y/(4.0*((1.0-Y))))*nhi*(x_min_fix),(Y/(4.0*((1.0-Y))))*nhi*(xhe1));
    double nhe2l = max((Y/(4.0*((1.0-Y))))*nhi*xhe2,(Y/(4.0*((1.0-Y))))*nhi*x_min_fix);
    double nhe3l = max((Y/(4.0*((1.0-Y))))*nhi*xhe3,(Y/(4.0*((1.0-Y))))*nhi*x_min_fix);//(Y/(4*((1.0-Y))))*nhi*xhe3;
   // double xhe1 = 1- xhe2 -xhe3;

    double ne = xh2i*nhi+(Y/(4.0*((1.0-Y))))*nhi*(xhe2+2*xhe3);
    double xe = xh2i+(xhe2+2*xhe3);
    double nu = (1.0 + 4.0 * Y) / (1.0 + Y + xe); //bolton
  //  double Ti = (2.0/3.0)*Tu*(nu*PROTONMASS)/(kb);
    double fac = -alpha*Ng/(pow(3e21,3)*13.6);
    /*isothermal case*/
#ifdef ISOTHERMAL
     fac = Ng/pow(3e21,3);
#endif
 
   //double fac = Ng/pow(3e21,3);
    double vol = 1.0/(4*PI*R*R*cellsize*hubble_inv*ascale(t));
    
    double nh1 = max(nhi*(1-xh2i),nhi*x_min_fix);
    double C,D,Tnu;
    //nu = pow((X*(1+xh2i)+Y*(1+xhe2+2*xhe3)/4),(-1));
   // nu = (1.0 + 4.0 * Y) / (1.0 + Y + ne); //bolton
   // Tnu = Ti/nu;
    double Gammab1, Gammab2;
    if (background==1)
    {
      Gammab1 =0.154e-12*timefac;
      Gammab2 =0.442e-18*timefac; // from haardt madau 2012
    }
    else{
      Gammab1 =0;
      Gammab2 = 0;
    }

    
    C = xhe1*ne*betahe1(Ti) + fac*vol*Gammahe1(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l,t)/((Y/(4.0*((1.0-Y))))*nhi) + xhe3*alphahe3(Ti)*ne + xhe1*Gammab1;
    D = betahe2(Ti)*ne + alphahe2(Ti)*ne + fac*vol*Gammahe2(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l,t)/(nhe2l) + Gammab2;
   // D = alphahe2(Ti)*ne;

    if (C !=C)
    {
        printf("He2 gone wrong\n");
        printf("gamma %0.12e %0.12e %0.12e\n", fac*vol, Gammahe1(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l,t),1/(nhe2l));
        printf("elec coll %0.12e,%0.12e, %0.12e\n",xhe1, betahe1(Ti),ne);
        exit(0);
    }
     //   printf("gamaa %0.12e\n",Gammahe2(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l,t));
    
    *(hrate) = (C-D*xhe2);
    if ((1+D*dti)==0)
    {
        printf("xhe2 denom is zero\n");
        exit(0);
    }
    return (xhe2 + C*dti)/(1+D*dti);
   // return xhe2 + (C-D*xhe2)*dti;
}

double xhe3solve(double Ti,double xh2i,double xhe1,double xhe2,double xhe3,double dti,double NH1,double Nhe1,double Nhe2,double R,double nhi,double *hrate,double t){
   
    double nhe1l = max((Y/(4.0*((1.0-Y))))*nhi*(x_min_fix),(Y/(4.0*((1.0-Y))))*nhi*(xhe1));
    double nhe2l = max((Y/(4.0*((1.0-Y))))*nhi*xhe2,(Y/(4.0*((1.0-Y))))*nhi*x_min_fix);
    double nhe3l = max((Y/(4.0*((1.0-Y))))*nhi*xhe3,(Y/(4.0*((1.0-Y))))*nhi*x_min_fix);
 
    //double xhe1 = 1- xhe2 -xhe3;

    double ne = xh2i*nhi+(Y/(4.0*((1.0-Y))))*nhi*(xhe2+2*xhe3);
    double xe = xh2i+(xhe2+2*xhe3);
    double nu = (1.0 + 4.0 * Y) / (1.0 + Y + xe); //bolton
  //  double Ti = (2.0/3.0)*Tu*(nu*PROTONMASS)/(kb);
    double fac = -alpha*Ng/(pow(3e21,3)*13.6);

#ifdef ISOTHERMAL
     fac = Ng/pow(3e21,3);
#endif
    /*isothermal case */
   // double fac = Ng/pow(3e21,3);
    double vol = 1.0/(4*PI*R*R*cellsize*hubble_inv*ascale(t));
    //double nh1 = nhi*(1-xh2i);
    double nh1 = max(nhi*(1-xh2i),nhi*x_min_fix);
    double C,D,Tnu;
   // nu = pow((X*(1+xh2i)+Y*(1+xhe2+2*xhe3)/4.0),(-1));
 //   nu = (1.0 + 4.0 * Y) / (1.0 + Y + ne); //bolton
  //  Tnu = Ti/nu;
    double Gammab;
    if (background==1)
    {
      Gammab = 0.442e-18*timefac; // from haardt madau 2012
    }
    else
      Gammab = 0;

    C = xhe2*ne*betahe2(Ti) + fac*vol*Gammahe2(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l,t)/((Y/(4.0*((1.0-Y))))*nhi)+ xhe2*Gammab;
    D = alphahe3(Ti)*ne;

    *(hrate) = (C-D*xhe3);
    if (C !=C)
    {
        printf("He3 gone wrong\n");
        printf("gamma %0.12e %0.12e %0.12e \n", fac*vol, Gammahe2(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l,t),xhe2/(nhe2l));
        printf("elec coll %0.12e %0.12e, %0.12e\n", xhe2,ne,betahe2(Ti));
        exit(0);
    }
    if ((1+D*dti)==0)
    {
        printf("xhe3 denom is zero\n");
        exit(0);
    }
   return (xhe3 + C*dti)/(1+D*dti);
   // return xhe3 + (C-D*xhe3)*dti;
}

double Tsolve(double Ti,double NH1,double Nhe1,double Nhe2,double xh2i,double xhe1,double xhe2,double xhe3,double nhi,double dti,double R,double tsim, double *trate, double *xrate,double t){
    
   
    double fac = -alpha*Ng/(pow(3e21,3)*13.6);
    double fh = fac*1.6e-12; // 1.6e-12 is factor ev to erg
    double vol = 1.0/(4*PI*R*R*cellsize*hubble_inv*ascale(t));
    
    //double K = (2.0/3)*(1.0/(nhi*kb)); // only H # cm^3erg^-1 K
    double nh1 = max(nhi*(1-xh2i),nhi*x_min_fix);
    double nhe1l = max((Y/(4.0*((1.0-Y))))*nhi*(x_min_fix),(Y/(4.0*((1.0-Y))))*nhi*(xhe1));
    double nhe2l = max((Y/(4.0*((1.0-Y))))*nhi*xhe2,(Y/(4.0*((1.0-Y))))*nhi*x_min_fix);
    double nhe3l = max((Y/(4.0*((1.0-Y))))*nhi*xhe3,(Y/(4.0*((1.0-Y))))*nhi*x_min_fix);
  
    double ne = xh2i*nhi+(Y/(4.0*((1.0-Y))))*nhi*(xhe2+2*xhe3);
    double xe = xh2i+(xhe2+2*xhe3);
    double n = nhi + (Y/(4.0*((1.0-Y))))*nhi + ne;
    // double nu = pow((X*(1+xh2i)+Y*(1+xhe2+2*xhe3)/4),(-1));
    double nu = (1.0 + 4.0 * Y) / (1.0 + Y + xe); //bolton

 //   double Ti = (2.0/3.0)*Tu*(nu*PROTONMASS)/(kb);
   // printf("t %lf\n",Ti);
    double K = (2.0/3.0)*(1.0/(n*kb)); // From Chen in units cm^3 erg^-1 K
    
    double T5 = Ti/100000.0;
    
    double f1 = 1.778e-29*timefac;  //#*6.24e11 # in erg cm^3 Myr^-1 K^-1 #4.05   #
    double f2 = 7.5e-19*timefac;  //#*6.24e11 # in erg cm^3 Myr^-1 #1.711e+11
    double f3 = 1.42e-27*timefac; //#*6.24e11# in erg cm^3 K^-1 #489.62
    double gff = 1.5 ;  //# from cen #1.1+ 0.34*(np.exp(-5.5+np.log(T)))**2
    double f4 = 5.65e-36*timefac; //#*6.24e11 #in erg cm^3 Myr^-1 K^-1 #1.289e-6
    double f5 = 1.27e-21*timefac;

    double f11 = 1.0;
    double f12 = 8*1.778e-29*timefac;
    double f21 = timefac*5.54e-17;
   
    double f51 = timefac*9.38e-22;
    double f52 = timefac*4.95e-22;
/*
    f11 = 1
    f12 =8*1.778*1e-29*timefac/(1.38e-16)
    f21 =5.54e-17*timefac/(1.38e-16)
    f22 =9.1e-27*timefac/(1.38e-16)
*/
    double Th1 = 157807;
    //double The = 631515;
    double lamb = 2*Th1/Ti;
    double lamb2 = 1263030.0/Ti;
    double omega_m = 0.308;
    double H0 = 67.8*1.05e-6;
    double rho = 1.0;
    //#fh = 1.04e+8
    double Hf = 1.05e-6;
    double omega_b = 0.0482;
    double rhocb   = 3.0*omega_b*pow(H0/timefac, 2.0)/(8.0*PI*GRAVITY);
    double atime   = ascale(t);
    double DELTA = 0.0;
    double RHOb =   rhocb / pow(atime, 3.0);
    double Xe= ne*PROTONMASS/(RHOb*(1.0+DELTA));
    
    if (K*nu==0 || nu!=nu)
    {
        printf("Knu is zero\n");
        printf("nu %0.12e, %0.12e,%0.12e,%0.12e\n",nu,xh2i,xhe2,xhe3);
        exit(0);
    }

    double zz  = 1.0/ascale(t) - 1.0;
    
    double c1 = (f1)*Ti*pow(lamb,1.965)/pow((1.0+pow((lamb/0.541),0.502)),2.697) ;// recombination cooling coeffecient erg cm^3 Myr^-1 from HnG
    double c2 = (f2)*pow((1.0+pow(T5,1.0/2.0)),-1)*exp(-118355.0/Ti);  //collisional excitation cooling coeffecient
    double c3 = gff*(f3)*pow(Ti,0.5); // Brehmstrallung/ free free emission
    //double c4 = (2.726*(1.0+zz)-Ti)*6.35e-41*timefac*omega_b*0.678*0.678*Xe*pow(1.0+zz,7);
    double c4 = (f4)*(Ti-2.73*(1.0+z))*pow(1.0+z,4.0);  //inverse compton
    double c5 = kb*Th1*GE(Ti);//from HnG in erg cm^3 Myr^-1 //f5*pow(Ti,0.5)*pow((1.0+pow(T5,0.5)),-1)*exp(-157809.1/Ti);// collisional Ionisation cooling

    double c11 = f11*kb*Ti*alphahe2(Ti); // recombination cooling He II erg cm^3 Myr^-1
    double c12 = f12*Ti*pow(lamb2,1.965)/pow((1.0+pow((lamb2/0.541),0.502)),2.697) ; // recombination cooling He III
        
    double c21 = f21*pow(Ti,-0.397)*pow((1.0+pow(T5,1.0/2.0)),-1)*exp(-473638.0/Ti);    //collisional excitation HeII erg cm^3 Myr^-1

    double c51 = kb*285335.0*betahe1(Ti);//f51*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-285335.4/Ti);
    double c52= kb*631515.0*betahe2(Ti);//f52*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-631515.0/Ti);
    //double C =  c1*pow((xh2i*nhi),2) +c2*xh2i*nhi*(nh1)+c3*pow((nhi*xh2i),2) +c4*(nhi*xh2i) + 2*Hubble(tsim)*Ti/(K*nu) ;//+ (Ti/(nhi*xh2i))*(*xrate);
    
    double C;

  //  C = c1*ne*xh2i*nhi + c2*ne*nh1 + c3*(nhi*xh2i+nhe2+4*nhe3)*ne + c4*ne+ 2*Hubble(tsim)*Ti/(K*nu);

   // C = c1*ne*xh2i*nhi + c11*ne*nhe2l + c12*ne*nhe3l +c2*ne*nh1 + c21*ne*nhe2l + c3*(nhi*xh2i+nhe2l+4*nhe3l)*ne + c4*ne + c5*ne*nh1+ c51*ne*nhe1l+ c52*ne*nhe2l +2*Hubble(tsim)*Ti/(K*nu)+ *(xrate)*(Ti/(ne+n))/(K*nu);
    C = c1*ne*xh2i*nhi + c11*ne*nhe2l + c12*ne*nhe3l +c2*ne*nh1 + c21*ne*nhe2l + c3*(nhi*xh2i+nhe2l+4*nhe3l)*ne + c4*ne + c5*ne*nh1+ c51*ne*nhe1l+ c52*ne*nhe2l;// +2*Hubble(tsim)*Ti/(K)+ *(xrate)*(Ti/(n))/(K);
    if (COMOVING == 1.0) {
        C = c1*ne*xh2i*nhi + c11*ne*nhe2l + c12*ne*nhe3l +c2*ne*nh1 + c21*ne*nhe2l + c3*(nhi*xh2i+nhe2l+4*nhe3l)*ne + c4*ne + c5*ne*nh1+ c51*ne*nhe1l+ c52*ne*nhe2l+2*Hubble(tsim)*Ti/(K) ; //+ *(xrate)*(Ti/(n))/(K);
    }
    if (C != C)
    {
        printf("Cooling terms invalid\n");
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
        printf("H %0.12e, %0.12e, %0.12e\n", 2*Hubble(tsim),Ti,(K*nu));
        exit(0);
    }
    //double C =0;
    double bheat1,bheat2,bheat3;
    if (background == 1)
    {
       bheat1 =0.989e-12*timefac*1.6e-12; // background photoheating from haardt madau 2012 erg Myr^-1
       bheat2 =0.112e-11*timefac*1.6e-12;
       bheat3 =0.269e-16*timefac*1.6e-12;
    }
    else{
       bheat1 =0;
       bheat2 =0;
       bheat3 =0;
    }
    
    double L = (fh*vol*Heat(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l,t)+ nh1*bheat1 +fh*vol*Heat2(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l,t)+nhe1l*bheat2+fh*vol*Heat3(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l,t)+nhe3l*bheat3 - C);

    //derivatives of cooling processes
    
    double A1 =  -(6.1712*1e10*f1*pow((1/Ti),1.965))/pow((1 + 784.353*pow((1/Ti),0.502)),2.697) + (
 6.79107*1e13*f1*pow((1/Ti),2.467))/pow((1 + 784.353*pow((1/Ti),0.502)),3.697);
    double A2 = (118355*exp(-118355/Ti)*f2)/((1 + 0.00316228*pow(Ti,0.5))*pow(Ti,2)) - (
 0.00158114*exp(-118355/Ti)*f2)/(pow((1 + 0.00316228*pow(Ti,0.5)),2)*pow(Ti,0.5));
    double A3 = (0.5*f3*gff)/pow(Ti,0.5);
    double A4 =   f4*pow((1+z),4);
    
    double A11 = -4.96095e-13*pow(1/Ti,1.654);
    double A12 = -(9.41467e11*f12*pow(1/Ti,1.965)/(pow((1 + 1573.42*pow((1/Ti),0.502)),2.697)))+ 2.07829e15*f12*pow(1/Ti,2.467)/pow((1 + 1573.42*pow((1/Ti),0.502)),3.697);
    double A21 = 473638.0*exp(-473638/Ti)*f21/((1+pow(T5,1.0/2))*pow(Ti,2.397))- 0.397*exp(-473638/Ti)*f21/((1+pow(T5,1.0/2))*pow(Ti,1.397)) - exp(-473638/Ti)*f21/(200*pow(10,0.5)*pow((1+pow(T5,1.0/2)),2)*pow(Ti,0.897));
    /*Add coll ion derivatives*/
    double dLdt = -(A1+ A2 +A3 +A4 +A11+ A12+ A21);// +2*Hubble(tsim)/(K*nu));//+(1/(nhi*xh2i))*(*xrate));
    //double dLdt =0;
   

  //  *(trate) = L*K;
    *(trate) = L/(nhi*PROTONMASS);
 //  return Ti+ dti*L*K;

    //return Tu + dti*L/(nhi*PROTONMASS) ;
    return Ti+ dti*L*K;

    //return Tu + dti*(L/(nhi*PROTONMASS))/(1-(dLdt/(nhi*PROTONMASS))*dti) ;
   //return (L*K*nu*dti + Ti)/(1+2*Hubble(tsim)*dti) ; // anninos implicit
}


long int xTiter(double dti,double NH1,double Nhe1, double Nhe2,double R,double nhi,double tsim,double ti,double *Tt, double *xh2t,double *xhe1t,double *xhe2t,double *xhe3t,double *dtrec){
 
    int storei=0;
    double frac_change = 0.0;
    //double ne_init;
    double dUU;
    double Trate=0;
    double Xrate=0;
    double he2rate = 0;
    double he3rate = 0;
    double nerate = 0;
    double ne_init,ne_up;
    double mu;
    double tcal;
    tcal = ti;
    for (int i=0;i<countl;i++){
#ifdef ISOTHERMAL
        Tt[i+1] =  Tt[i];
#else
     Tt[i+1] = (double) max(T2_MIN_FIX,Tsolve(Tt[i],NH1,Nhe1,Nhe2,xh2t[i],xhe1t[i],xhe2t[i],xhe3t[i],nhi,dti,R,tsim,&Trate,&nerate,tcal));
#endif

    dUU = (double) fabs(max(T2_MIN_FIX,Tt[i]+Trate*dti)-Tt[i]);
    
    if((double) max(dUU,fabs(Tt[i+1]-Tt[i]))/(Tt[i]+T_MIN)>lim  && dti>dt_min){
        return 100000000000;
    }
   else{
        frac_change = (double) max(frac_change,dUU);
    }

    xh2t[i+1] = (double) max(x_min_fix,xh2solve(Tt[i+1],xh2t[i],xhe1t[i],xhe2t[i],xhe3t[i],dti,NH1,Nhe1,Nhe2,R,nhi,&Xrate,tcal));

/*complete this exactly*/

     if((double) max(fabs(xh2t[i+1]-xh2t[i]),Xrate*dti)/(xh2t[i]+X_FM) >lim && dti>dt_min){
         return 100000000001;
     }
     else{
         frac_change = (double) max(frac_change,fabs(xh2t[i+1]-xh2t[i])/(xh2t[i])); //check
     }

#ifdef NOHelium
        xhe1t[i+1] = xhe1t[i];
        xhe2t[i+1] = xhe2t[i];
        xhe3t[i+1] = xhe3t[i];
   
#else
    xhe1t[i+1] = (double) max(x_min_fix,xhe1solve(Tt[i+1],xh2t[i+1],xhe1t[i],xhe2t[i],xhe3t[i],dti,NH1,Nhe1,Nhe2,R,nhi,tcal));

    xhe2t[i+1] = (double) max(x_min_fix,xhe2solve(Tt[i+1],xh2t[i+1],xhe1t[i],xhe2t[i],xhe3t[i],dti,NH1,Nhe1,Nhe2,R,nhi,&he2rate,tcal));

    xhe3t[i+1] = (double) max(x_min_fix,xhe3solve(Tt[i+1],xh2t[i+1],xhe1t[i],xhe2t[i+1],xhe3t[i],dti,NH1,Nhe1,Nhe2,R,nhi,&he3rate,tcal));

    //xhe1t[i+1] = (double) max(x_min_fix,1-xhe2t[i+1]-xhe3t[i+1]);


        if (xhe1t[i+1]+xhe2t[i+1]+xhe3t[i+1]!=1e0)
        {
            // printf("He fracs don't add up to 1  but %0.12e\n",xhe1t[i+1]+xhe2t[i+1]+xhe3t[i+1]);
            if(max(max(xhe2t[i+1],xhe3t[i+1]),xhe1t[i+1])==xhe2t[i+1] ){
                xhe2t[i+1] = 1-xhe1t[i+1] -xhe3t[i+1] ;
            }
            else if(max(max(xhe2t[i+1],xhe3t[i+1]),xhe1t[i+1])==xhe3t[i+1]){
                xhe3t[i+1] = 1-xhe1t[i+1] -xhe2t[i+1] ;
            }
  
        }

#endif

    ne_init = xh2t[i]*nhi+(Y/(4.0*((1.0-Y))))*nhi*(xhe2t[i]+2*xhe3t[i]);
    ne_up = xh2t[i+1]*nhi+(Y/(4.0*((1.0-Y))))*nhi*(xhe2t[i+1]+2*xhe3t[i+1]);
    nerate =  Xrate*nhi+(Y/(4.0*((1.0-Y))))*nhi*(he2rate+2*he3rate);

    if(fabs(ne_up - ne_init)/(ne_init+X_FM)>lim && dti >dt_min){
    //printf("change %0.12e %0.12e\n", fabs(ne_up - ne_init)/(ne_init+X_FM), lim);
        return 100000000002;
    }
    else {
        frac_change = (double) max(frac_change,fabs(ne_up - ne_init)/(ne_init));
    }

    /* advance timestep */
    tcal = ti + dti;
    ti = tcal ;
    /* set next time step */
    if(frac_change <0.5 && dti>dt_min){  // has value 0.05 in paper but 0.5 here
    *(dtrec)  = 2*dti;
    //    dti = 2*dti;
    }
    else{
    *(dtrec) = dti;
    //  printf("dti %e\n",dti );
    }
   
  
   //     printf("tsim %lf\n", tcal);
    if (tcal >= tsim){
            storei = i+1;
            printf("tsim %lf\n", tcal);
            break ;
    }
   dti = (double) min(*(dtrec),tsim-tcal);

   }
    double nu = (1.0+4.0*Y)/(1.0+Y+xh2t[storei]+xhe2t[storei]+2*xhe3t[storei]);
    printf("T %lf %d\n",Tt[storei],storei);
    return storei;
}


double evolve_eqns(double *nh,double **nh2,double **NH,double **nhe1,double **nhe2,double **nhe3,double **NHe1,double **NHe2,double **T,int ncells,double ti,double tsim,double **dt, int j,double Dt){
  //double dt = Dt;
  
  long int out;
  double *Tt =(double*)malloc(countl*sizeof(double));
  double *xh2t =(double*)malloc(countl*sizeof(double));
  double *xhe1t =(double*)malloc(countl*sizeof(double));
  double *xhe2t =(double*)malloc(countl*sizeof(double));
  double *xhe3t =(double*)malloc(countl*sizeof(double));
  double dt_rec;
  double *Rcell =(double*)malloc(ncells*sizeof(double));
    
    
  for(int k=0; k<ncells;k++){
    printf("dt new %e %d %d\n",dt[k][j],k,j);
    //dt[k][j+1] = dt[k][j];
    printf("cell %d\n",k);
    double nu;
    xh2t[0] = nh2[k][j]/nh[k];
      
    #ifdef NOHelium
        xhe2t[0] = 0.0;
        xhe3t[0] = 0.0;
        xhe1t[0] = 0.0;
#else
    xhe2t[0] = nhe2[k][j]/((Y/(4.0*((1.0-Y))))*nh[k]);
    xhe3t[0] = nhe3[k][j]/((Y/(4.0*((1.0-Y))))*nh[k]);
    xhe1t[0] = 1.0 -xhe2t[0] -xhe3t[0];
 
    #endif
       
    nu = (1.0+4.0*Y)/(1.0+Y+xh2t[0]+xhe2t[0]+2*xhe3t[0]);
  //  Tt[0] = (3.0/2.0)*T[k][j]*kb/(nu*PROTONMASS);
      Tt[0] = T[k][j];
    Rcell[k] = (Roffset + k*cellsize)*hubble_inv*ascale(tsim) ;
  
    out = xTiter(dt[k][j],NH[k][j],NHe1[k][j],NHe2[k][j],Rcell[k],nh[k],tsim,ti,Tt,xh2t,xhe1t,xhe2t,xhe3t,&dt_rec);
    // equivalent of cool_step routine
    //printf("%d\n",out);
  
    while((out==100000000000 || out==100000000001 || out ==100000000002 || out==100000000003)){
    
        dt[k][j] =(double) max(dt_min,dt[k][j]/2);
        printf("%ld %e\n",out-100000000000, dt[k][j]);

        out = xTiter(dt[k][j],NH[k][j],NHe1[k][j],NHe2[k][j],Rcell[k],nh[k],tsim,ti,Tt,xh2t,xhe1t,xhe2t,xhe3t,&dt_rec);

    }
    printf("%ld \n",out);
    
  
    nh2[k][j+1] = pow(ascale(ti),3)*nh[k]*xh2t[out]/pow(ascale(tsim),3);
    nhe2[k][j+1] = (Y/(4.0*((1.0-Y))))*pow(ascale(ti),3)*nh[k]*xhe2t[out]/pow(ascale(tsim),3);
    nhe3[k][j+1] = (Y/(4.0*((1.0-Y))))*pow(ascale(ti),3)*nh[k]*xhe3t[out]/pow(ascale(tsim),3);
    nhe1[k][j+1] = (Y/(4.0*((1.0-Y))))*pow(ascale(ti),3)*nh[k]*xhe1t[out]/pow(ascale(tsim),3);//(Y/(4*((1.0-Y))))*nh[k] - nhe2[k][j+1] - nhe3[k][j+1];
    nu = (1.0+4.0*Y)/(1.0+Y+xh2t[out]+xhe2t[out]+2*xhe3t[out]);
    //T[k][j+1] = (2.0/3.0)*Tt[out]*(nu*PROTONMASS)/(kb);
      T[k][j+1] = Tt[out];
      //dt[k][j+1] = Dt;//(double) min((dt_rec),tsim-ti); // check this criterion
    
}
  for(int k=0;k<ncells;k++){
  NH[k][j+1] = cumsum(nh,nh2,k,j+1)*cellsize*hubble_inv*kpc*nm2*ascale(tsim); // in cm^-2
  NHe1[k][j+1] = hecumsum(nhe1,k,j+1)*cellsize*hubble_inv*kpc*nm2*ascale(tsim) ;
  NHe2[k][j+1] = hecumsum(nhe2,k,j+1)*cellsize*hubble_inv*kpc*nm2*ascale(tsim) ;
  nh[k] = pow(ascale(ti),3)*nh[k]/pow(ascale(tsim),3);
}
  
  return dt_rec;
}




void
absorb_along_lines_of_sight (double *Tau, double *Temp,double *NHI, double *Vpec)
{
  int ncells = boxsize*hubble_inv/(cellsize*hubble_inv);
  double dz, dv, b, fac,aa;
  long long  bin, k;
  double profile_H1, T0, T1, T2;
  double Hubble_h = 0.678;
  double Omega = 0.308;
  double Lambda = 0.692;
  double Time = 1.0 / (1.0 + z);
  double UnitLength_in_cm = 3.085678e21;
  double UnitMass_in_g = 1.989e43;
  double UnitVelocity_in_cm_per_s = 1e5;
  double UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  double Hubble = HUBBLE * UnitTime_in_s;

  double H_a = Hubble * sqrt (Omega / (Time * Time * Time) +
           Lambda + (1 - Omega - Lambda) / (Time * Time));

  double Wmax = Time * H_a * boxsize/Hubble_h;

  printf("Omega=%e\n",Omega);
  printf("Lambda=%e\n",Lambda);
  printf("Time=%e\n",Time);
  printf("z=%e\n",1.0/Time-1.0);
  printf("H_a=%e\n",H_a);
  printf("boxsize=%e\n",boxsize);
  
  printf("T :%0.12e\n",*Temp );
  printf("NH :%0.12e \n", NHI[2046]);
  printf("V : %0.12e\n",*Vpec );
  printf ("Wmax= %g  Time=%g\n", Wmax, Time);

  dz = boxsize / ncells;


    for (bin = 0; bin < ncells; bin++)
    {

      //printf("m=%lld\n", m);
      //printf("bin=%lld\n", bin);
      
      Tau[bin] = 0;

      for (k = 0; k < ncells; k++)
      {
        dv = (k - bin);

        while (dv < -ncells / 2)
          dv += ncells;
        while (dv > ncells / 2)
          dv -= ncells;

        dv = (dv * Wmax / ncells + Vpec[k]) *
          UnitVelocity_in_cm_per_s;
        
        b = sqrt (2 * BOLTZMANN * Temp[k] / PROTONMASS);
      
        aa = GAMMA_LYA_H1*LYMAN_ALPHA/(4.0*PI)/b;

        T0 = dv * dv / b / b;
        T1 = exp(-T0);
        T2 = 1.5/T0;
        
        profile_H1 = (T0 < 1.0e-6)
          ? T1 : T1-aa/sqrt(PI)/T0*(T1*T1*(4.0*T0*T0+7.0*T0+4.0+T2)-T2-1.0);
        
        Tau[bin] += NHI[k]  *profile_H1 / b * dz;
   
      }
  
    }
  

  /* multiply with correct prefactors */

  /*  to get things into cgs units */
  fac = 1.0 / pow (UnitLength_in_cm, 2);
  fac *= Hubble_h * Hubble_h;

  /* Ly-alpha cross section */
  fac *= OSCILLATOR_STRENGTH * PI * LYMAN_ALPHA *
    sqrt (3 * THOMPSON / (8 * PI));

  fac *= Cv / (Time * Time) / sqrt (PI);

  /* Note: For HeII, the oscillator strength is equal to that of HI,
     and the Lyman-alpha wavelength is 4 times shorter */

 
    for (bin = 0; bin < ncells; bin++)
    {
      Tau[bin] *= fac;
     // printf("%le\n",Tau[bin]);
    }

   // Los[m].boxsize = boxsize;
   // Los[m].VelScale = Wmax;
   // Los[m].Time = Time;
  

  printf("Tau=%e\n", Tau[0]);
  
}






int main(int argc,char **argv)
{

  /*For each position, solve euler to get nh(xpos,t) and T(xpos,t) where t is tsim= time after each global dt */
  clock_t t;
      t = clock();
      
      
    double  omega_m = 0.308;
    double  H0 = 67.8*1.05e-6;
    
    double ai = 1.0/(1.0+z);
    double tz = (1.0/H0)*(2.0/3.0)*pow(1.0-omega_m,-0.5)*asinh(pow((1.0-omega_m)/omega_m,0.5)*pow(ai,1.5));
    printf("ti %lf\n",tz );

  
    double ti = tz- atof(argv[1]) ;
    double tf =tz;


  
   double Dt = 3261.6*(epsilon/0.1)*(cellsize*hubble_inv/10.0)*1e-6 ;
    //Dt= 1e-3;
   printf("initial timestep %e\n",Dt);
    
   int ncells = boxsize*hubble_inv/(cellsize*hubble_inv);


/* initialize nh,nh2,NH and T grid */
   double *nh = (double*)malloc(ncells*sizeof(double));
    
   int r  = ncells;
   int c = count ;
   double **nh2 = (double **)malloc(r * sizeof(int *));
      for (int i=0; i<r; i++)
          nh2[i] = (double *)malloc(c * sizeof(int));
   double **nhe2 = (double **)malloc(r * sizeof(int *));
        for (int i=0; i<r; i++)
            nhe2[i] = (double *)malloc(c * sizeof(int));
   double **nhe3 = (double **)malloc(r * sizeof(int *));
        for (int i=0; i<r; i++)
            nhe3[i] = (double *)malloc(c * sizeof(int));
   double **NH = (double **)malloc(r * sizeof(int *));
      for (int i=0; i<r; i++)
    NH[i] = (double *)malloc(c * sizeof(int));
   double **NHe1 = (double **)malloc(r * sizeof(int *));
        for (int i=0; i<r; i++)
        NHe1[i] = (double *)malloc(c * sizeof(int));
   double **NHe2 = (double **)malloc(r * sizeof(int *));
        for (int i=0; i<r; i++)
        NHe2[i] = (double *)malloc(c * sizeof(int));
   double **T = (double **)malloc(r * sizeof(int *));
      for (int i=0; i<r; i++)
          T[i] = (double *)malloc(c * sizeof(int));
   double **dt = (double **)malloc(r * sizeof(int *));
        for (int i=0; i<r; i++)
            dt[i] = (double *)malloc(c * sizeof(int));
// intermediate
   double **nhe1 = (double **)malloc(r * sizeof(int *));
        for (int i=0; i<r; i++)
            nhe1[i] = (double *)malloc(c * sizeof(int));


   double *Tau = (double*)malloc(ncells*sizeof(double));
   double *Vpec = (double*)malloc(ncells*sizeof(double));
   // double *rho = (double*)malloc(ncells*sizeof(double));
   // double **rho = (double **)malloc(ncells * sizeof(int *));
   //     for (int i=0; i<r; i++)
   //         rho[i] = (double *)malloc(5000 * sizeof(int));
/*
FILE *myFile,*myFile2, *myFile3,*myFile4,*myFile5;
    myFile = fopen("LOS23_nh.txt", "r");
    myFile2 = fopen("LOS23_nh1.txt", "r");
    myFile3 = fopen("LOS23_T.txt", "r");
    myFile4 = fopen("LOS23_Vpec.txt", "r");
    myFile5 = fopen("LOS23_rho.txt", "r");
    double nharr[2048];
    double nh1arr[2048];
    double Tarr[2048];
    double vpec[2048];
    double rho[2048];
    
    for (int i = 0; i < 2048; i++)
    {
        fscanf(myFile, "%lf", &nharr[i]);
        fscanf(myFile2, "%lf", &nh1arr[i]);
        fscanf(myFile3, "%lf", &Tarr[i]);
        fscanf(myFile4, "%lf", &vpec[i]);
        fscanf(myFile5, "%lf", &rho[i]);
    }
 */
    /*
    for (int i = 0; i < 2048; ++i)
    {
      for (int j = 0; j < 5000; ++j)
      {
        fscanf(myFile5,"%lf",&rho[i][j]);
      }
      
    }
    */
/*
    for (int i = 0; i < 2048; i++)
    {
        printf("Number is: %0.12e\n", numberArray[i]);
    }
*/

   printf("nh at z %lf %0.12e\n",z,getnh(z));

  /* Initial values*/
  for (int i = 0; i < ncells; ++i){
      nh[i] =1e-3;//getnh(z);
  nh2[i][0] = max(x_min_fix*nh[i],INITIAL_HII*nh[i]);//nharr[i]- 1e-4*nharr[i];//nh[i]- nh1arr[i];//x_min_fix*nh[i];
  nhe2[i][0] = max(x_min_fix*(Y/(4.0*((1.0-Y))))*nh[i],INITIAL_HEII*(Y/(4.0*((1.0-Y))))*nh[i]);
  nhe3[i][0] = max(x_min_fix*(Y/(4.0*((1.0-Y))))*nh[i],INITIAL_HEIII*(Y/(4.0*((1.0-Y))))*nh[i]);
  nhe1[i][0] = (Y/(4.0*((1.0-Y))))*nh[i] - nhe2[i][0] - nhe3[i][0];
  T[i][0] = INITIAL_TEMP;//Tarr[i];
  NH[i][0] = cumsum(nh,nh2,i,0)*cellsize*hubble_inv*kpc*nm2*ascale(ti); // in nm^-2
  NHe1[i][0] = hecumsum(nhe1,i,0)*cellsize*hubble_inv*kpc*nm2*ascale(ti) ;
  NHe2[i][0] = hecumsum(nhe2,i,0)*cellsize*hubble_inv*kpc*nm2*ascale(ti) ;
  dt[i][0] = Dt;
  Vpec[i] =  0.0;
}

  
/* printing 2d array
  for (int i =0;i<ncells;++i){
  printf("[");
  for (int j =0;j<count;j++){
  printf("%lf,",nh2[i][j]);
}
printf("]\n");
}
*/

/*Evolving over global timestep*/
double tsim ;
int storej;
double dt_new;
double RIF,RIFold;
double NHION;
int cks;

for (int ck = 0; ck < ncells; ++ck)
    {
        if (nh2[ck][0]/nh[ck]<=0.5)
        {
            RIF = (Roffset+(ck)*cellsize)*hubble_inv*ascale(ti);
            cks =ck;
            break;
        }
    }
  

    printf("RIF %0.12e\n",RIF );

    NHION = NH[cks+1][0];//nm^-2//*(pow(3e21,2))/1e-14 ; // in kpc^-2
    printf("%0.12e %d\n",NHION,cks );
    // Dt = 4*PI*pow(RIF*3e28,2)*NHION*pow(Ng,-1); //in Myr
    printf("Dt %0.12e\n",Dt );
/*
    for (int i = 0; i < ncells; ++i)
{
  dt[i][0] =Dt;
}
*/

for (int j=0;j<count;j++){
  printf("Dt %0.12e\n",Dt );
  tsim = ti + Dt;
 
  /*evolve all cells for a  global timestep*/
  evolve_eqns(nh,nh2,NH,nhe1,nhe2,nhe3,NHe1,NHe2,T,ncells,ti,tsim,dt,j,Dt); // my equivalent of subroutine rt_solve_cooling

    for (int ck = 0; ck < ncells; ++ck)
    {
        if (nh2[ck][j]/nh[ck]<=0.5)
        {
            RIFold = RIF;
            RIF = (Roffset+(ck)*cellsize)*hubble_inv*ascale(tsim);
            cks =ck;
            break;
        }
    }
    printf("RIF %0.12e\n",RIF );
    
    ti = tsim;
    printf("%lf\n",tsim);
    if(fabs(tf-tsim)<=Dt){
      printf("%lf,%lf\n",tf,tsim);
      storej = j+1;
      break;
}
  
   //if (cks==ncells-1)

   //if(fabs(RIFold - RIF)<epsilon)
  /*   if(Dt*j>1e-1)
   {
        Dt = 1e-2;//3261.6*(epsilon/0.1)*(cellsize*hubble_inv/10.0)*1e-6;//nm^-2//*(pow(3e21,2))/1e-14 ; // in kpc^-2
   }
 
   */
 // else{
   NHION = NH[cks+1][j];//nm^-2//*(pow(3e21,2))/1e-14 ; // in kpc^-2
   printf("%0.12e %d\n",NHION,cks );
   Dt = 4*PI*pow(RIF*3e28,2)*NHION*pow(Ng,-1); //in Myr
   printf("Dt %0.12e\n",Dt );
  // }
  // */
   for (int k= 0; k < ncells; ++k)
   {
     dt[k][j+1]= Dt;
   }
   
}

double *Tf = (double*)malloc(ncells*sizeof(double));
double *NH1f = (double*)malloc(ncells*sizeof(double));
double dense ;
double UnitMass_in_g = 1.989e43;
double Hubble_h = 0.678;
dense = pow(3e21,3)/(kpc*nm2*1.989e43*pow(Hubble_h,2));//*Hubble_h*Hubble_h;

for (int i = 0; i < ncells; ++i)
{
  Tf[i] = T[i][storej];
  NH1f[i] = 1*(1-nh2[i][storej]/nh[i])* (UnitMass_in_g / PROTONMASS / Hubble_h) ;//NH[i][storej];
}
absorb_along_lines_of_sight(Tau,Tf,NH1f,Vpec);

printf("[");
//printf("%d\n",storej);

for (int i=0;i<ncells;i++){

printf("%e,",1-nh2[i][storej]/nh[i]);
}
printf("]\n");




for (int i=0;i<ncells;i++){

printf("%e,",T[i][storej]);
}
printf("]\n");

for (int i=0;i<ncells;i++){

printf("%e,",Tau[i]);
}
printf("]");

t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
  
    printf("fun() took %f seconds to execute \n", time_taken);


//FILE *fptr = fopen("./boltA_w_scalefac.txt","w");
   // FILE *fptr = fopen("./molweighdav_w_scalefac_old.txt","w");
    FILE *fptr = fopen("./isothermaltestnewcode.txt","w");
//    FILE *fptr = fopen("./bolt_w_scalefac.txt","w");
      for (int l = 0; l < ncells; ++l)
    {
#ifdef NOHelium
        fprintf(fptr,"%lf %0.12e %0.12e\n", ((double)l*cellsize+Roffset),(nh[l]-nh2[l][storej])/nh[l],T[l][storej]);
#else
    fprintf(fptr,"%lf %0.12e %0.12e %0.12e %0.12e %0.12e\n", ((double)l*cellsize+Roffset),(nh[l]-nh2[l][storej])/nh[l],T[l][storej],nhe1[l][storej]/((Y/(4.0*((1.0-Y)))*nh[l])),nhe2[l][storej]/((Y/(4.0*((1.0-Y)))*nh[l])),nhe3[l][storej]/((Y/(4.0*((1.0-Y)))*nh[l])));
#endif
    }

}



