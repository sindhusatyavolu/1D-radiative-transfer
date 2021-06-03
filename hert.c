#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h> 
#include<time.h>

#ifndef PI
#define PI 3.1415
#endif
#define cellsize 2
#define boxsize 200
#define count 100000
#define countl 1000000000
#define Ng 3.15e13*5e56
#define alpha -1.5
#define kb 1.38e-16
#define max(x,y) (((x) >= (y)) ? (x) : (y))
#define min(x,y) (((x) <= (y)) ? (x) : (y))
#define T_MIN 0.1
#define T2_MIN_FIX 1e-2
#define X_MIN  1e-20
#define X_FM 1e-6
#define X 0.76
#define Y 0.24
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



double GE(double Ti){
    double Th1 = 157807; 
    double lam = 2*Th1/Ti;
    double T5 =Ti/100000 ;
    return 3.15e13*5.85e-11*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-157809.1/Ti);
    //return 3.15e13*21.11*pow(Ti,-3/2)*exp(-lam/2)*pow(lam,-1.089)/pow(1+pow((lam/0.354),0.874),1.101) ;
}

double betahe1(double Ti){
    double T5 =Ti/100000 ;
    return 3.15e13*2.38e-11*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-285335.4/Ti);
}

double betahe2(double Ti){
    double T5 =Ti/100000 ;
    return 3.15e13*5.68e-12*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-631515.0/Ti);
}


double ralpha(double Ti){
    double Th1 = 157807; 
    double lam = 2*Th1/Ti;

    return 3.15e13*1.269e-13*pow(lam,1.503)/pow((1+pow(lam/0.522,0.47)),1.923);
}

double alphahe2(double Ti){

    double lam = 570670/Ti;

    return 3.15e13*3e-14*pow(lam,0.654);

}


double alphahe3(double Ti){
  
    double lam = 1263030/Ti;

    return 3.15e13*2.538e-13*pow(lam,1.503)/pow((1+pow(lam/0.522,0.47)),1.923);
}


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
}

double sigmahe1(double f){
    double sigma_0= 9.492e-2 ; // in nm^2
    double nu_0 = 13.61 ;  
    double y0 = 0.4434;
    double x =  f/nu_0 - y0;
    double yw = 2.039;
    double P =  3.188;
    double ya = 1.469 ;
    double y1  = 2.136 ;
    double y = sqrt(x*x+y1*y1);
    return sigma_0*(pow((x-1),2)+pow(yw,2))*pow(y,0.5*P-5.5)/pow(1+sqrt(y/ya),P);
}

double sigmahe2(double f){
    double sigma_0=1.369 ; // in nm^2
    double nu_0 = 1.720 ;  
    double y0 = 0;
    double x =  f/nu_0 - y0;
    double yw = 0;
    double P =  2.963;
    double ya = 32.88 ;
    double y1  = 0 ;
    double y = sqrt(x*x+y1*y1);
    return sigma_0*(pow((x-1),2)+pow(yw,2))*pow(y,0.5*P-5.5)/pow(1+sqrt(y/ya),P);
}

double PH1(double x,double nh1,double nhe1l,double nhe2l){
	double qH = exp(-sigma(x)*nh1*cellsize*3e7);
	double qHe1 = exp(-sigmahe1(x)*nhe1l*cellsize*3e7);
    double 	qHe2 = exp(-sigmahe2(x)*nhe2l*cellsize*3e7)	;
    double	tau_tot = sigma(x)*nh1*cellsize*3e7+sigmahe1(x)*nhe1l*cellsize*3e7+ sigmahe2(x)*nhe2l*cellsize*3e7;
    double	D = (1-qH)*qHe2*qHe1 +(1-qHe1)*qH*qHe2 + (1-qHe2)*qH*qHe1;

	return (1-qH)*qHe2*qHe1*(1-exp(-tau_tot))/D;
}

double PHe1(double x,double nh1,double nhe1l,double nhe2l){
    double qH = exp(-sigma(x)*nh1*cellsize*3e7);
    double qHe1 = exp(-sigmahe1(x)*nhe1l*cellsize*3e7);
    double  qHe2 = exp(-sigmahe2(x)*nhe2l*cellsize*3e7)  ;
    double  tau_tot = sigma(x)*nh1*cellsize*3e7+sigmahe1(x)*nhe1l*cellsize*3e7+ sigmahe2(x)*nhe2l*cellsize*3e7;
    double  D = (1-qH)*qHe2*qHe1 +(1-qHe1)*qH*qHe2 + (1-qHe2)*qH*qHe1;

    return (qH)*qHe2*(1-qHe1)*(1-exp(-tau_tot))/D;
}

double PHe2(double x,double nh1,double nhe1l,double nhe2l){
    double qH = exp(-sigma(x)*nh1*cellsize*3e7);
    double qHe1 = exp(-sigmahe1(x)*nhe1l*cellsize*3e7);
    double  qHe2 = exp(-sigmahe2(x)*nhe2l*cellsize*3e7)  ;
    double  tau_tot = sigma(x)*nh1*cellsize*3e7+sigmahe1(x)*nhe1l*cellsize*3e7+ sigmahe2(x)*nhe2l*cellsize*3e7;
    double  D = (1-qH)*qHe2*qHe1 +(1-qHe1)*qH*qHe2 + (1-qHe2)*qH*qHe1;

    return (qH)*(1-qHe2)*qHe1*(1-exp(-tau_tot))/D;
}

double Gamma(double NH,double Nhe1, double Nhe2,double nh1,double nhe1l,double nhe2l){
    
    double u[] ={13.6, 14.4059, 15.2595, 16.1636, 17.1214, 18.1359, 19.2105, 20.3488,21.5545, 22.8317, 24.1846, 25.6176, 27.1356, 28.7435, 30.4466, 32.2507, 34.1617, 36.1859, 38.33, 40.6012, 43.007, 45.5553, 48.2546, 51.1139, 54.1426,57.3507, 60.749, 64.3486, 68.1615, 72.2003, 76.4784, 81.0101, 85.8102, 90.8948, 96.2806, 101.986, 108.029, 114.43, 121.21, 128.392, 136};
    // SPACING = (b/a)^(1/n) , n= 40, b = 136, a =13.6 
    double fsum=0 ;
    double x0 = u[0];
    
    int len = sizeof(u)/sizeof(u[0]);
    double df = 1;//0.5*pow(u[len-1]/x0,1.0/(len-1));
   // fsum = integral(&u[0],len,NH,nh1,nhe1,nhe2);
    //for(int i=0;i<sizeof(u)/sizeof(u[0]);i++){
 

    for(int i=0;i<40;i++){
     fsum = fsum + df*pow(u[i]/x0,alpha-1)*exp(-sigma(u[i])*NH-sigmahe1(u[i])*Nhe1-sigmahe2(u[i])*Nhe2)*PH1(u[i],nh1,nhe1l,nhe2l) ;	
     
}
 
  /*isothermal case*/
   //  fsum  = df*exp(-sigma(u[0])*NH)*PH1(u[0],nh1,nhe1,nhe2) ;  
  //  fsum  = df*exp(-sigma(u[0])*NH-sigmahe1(u[0])*Nhe1-sigmahe2(u[0])*Nhe2)*PH1(u[0],nh1,nhe1,nhe2) ;
    return fsum;
}

double Gammahe1(double NH,double Nhe1, double Nhe2,double nh1,double nhe1l,double nhe2l){
    
    //double u[] ={13.6, 14.4059, 15.2595, 16.1636, 17.1214, 18.1359, 19.2105, 20.3488,21.5545, 22.8317, 24.1846, 25.6176, 27.1356, 28.7435, 30.4466, 32.2507, 34.1617, 36.1859, 38.33, 40.6012, 43.007, 45.5553, 48.2546, 51.1139, 54.1426,57.3507, 60.749, 64.3486, 68.1615, 72.2003, 76.4784, 81.0101, 85.8102, 90.8948, 96.2806, 101.986, 108.029, 114.43, 121.21, 128.392, 136};
    double u[] = {54.5,57.7293, 61.15, 64.7734, 68.6114, 72.6769, 76.9833, 81.5448, 86.3767, 91.4948, 96.9162, 102.659, 108.742, 115.185, 122.01, 129.24, 136.898, 145.01, 153.602, 162.703, 172.344, 182.556, 193.373, 204.831, 216.968};
    // SPACING = (b/a)^(1/n) , n= 40, b = 136, a =13.6 
    double fsum=0 ;
    double x0 = u[0];
    
    int len = sizeof(u)/sizeof(u[0]);
    double df = 1;//0.5*pow(u[len-1]/x0,1.0/(len-1));
   // fsum = integral(&u[0],len,NH,nh1,nhe1,nhe2);
    //for(int i=0;i<sizeof(u)/sizeof(u[0]);i++){
 

    for(int i=0;i<len;i++){
     fsum = fsum + df*pow(u[i]/x0,alpha-1)*exp(-sigma(u[i])*NH-sigmahe1(u[i])*Nhe1-sigmahe2(u[i])*Nhe2)*PHe1(u[i],nh1,nhe1l,nhe2l) ;   
     
}
 
  /*isothermal case*/
   //  fsum  = df*exp(-sigma(u[0])*NH)*PH1(u[0],nh1,nhe1,nhe2) ;  
   
 //   fsum  = df*exp(-sigma(u[0])*NH-sigmahe1(u[0])*Nhe1-sigmahe2(u[0])*Nhe2)*PHe1(u[0],nh1,nhe1,nhe2) ;
    return fsum;
}

double Gammahe2(double NH,double Nhe1, double Nhe2,double nh1,double nhe1l,double nhe2l){
    
    //double u[] ={13.6, 14.4059, 15.2595, 16.1636, 17.1214, 18.1359, 19.2105, 20.3488,21.5545, 22.8317, 24.1846, 25.6176, 27.1356, 28.7435, 30.4466, 32.2507, 34.1617, 36.1859, 38.33, 40.6012, 43.007, 45.5553, 48.2546, 51.1139, 54.1426,57.3507, 60.749, 64.3486, 68.1615, 72.2003, 76.4784, 81.0101, 85.8102, 90.8948, 96.2806, 101.986, 108.029, 114.43, 121.21, 128.392, 136};
    double u[] = {54.5,57.7293, 61.15, 64.7734, 68.6114, 72.6769, 76.9833, 81.5448, 86.3767, 91.4948, 96.9162, 102.659, 108.742, 115.185, 122.01, 129.24, 136.898, 145.01, 153.602, 162.703, 172.344, 182.556, 193.373, 204.831, 216.968};
    // SPACING = (b/a)^(1/n) , n= 40, b = 136, a =13.6 
    double fsum=0 ;
    double x0 = u[0];
    
    int len = sizeof(u)/sizeof(u[0]);
    double df = 1;//0.5*pow(u[len-1]/x0,1.0/(len-1));
   // fsum = integral(&u[0],len,NH,nh1,nhe1,nhe2);
    //for(int i=0;i<sizeof(u)/sizeof(u[0]);i++){
 

    for(int i=0;i<len;i++){
     fsum = fsum + df*pow(u[i]/x0,alpha-1)*exp(-sigma(u[i])*NH-sigmahe1(u[i])*Nhe1-sigmahe2(u[i])*Nhe2)*PHe2(u[i],nh1,nhe1l,nhe2l) ;   
     
}
 
  /*isothermal case*/
 //    fsum  = df*exp(-sigma(u[0])*NH-sigmahe1(u[0])*Nhe1-sigmahe2(u[0])*Nhe2)*PHe2(u[0],nh1,nhe1,nhe2) ;  
   

    return fsum;
}

double Heat(double NH,double Nhe1,double Nhe2,double nh1,double nhe1l,double nhe2l){
    
    double u[] ={13.6, 14.4059, 15.2595, 16.1636, 17.1214, 18.1359, 19.2105, 20.3488,21.5545, 22.8317, 24.1846, 25.6176, 27.1356, 28.7435, 30.4466, 32.2507, 34.1617, 36.1859, 38.33, 40.6012, 43.007, 45.5553, 48.2546, 51.1139, 54.1426,57.3507, 60.749, 64.3486, 68.1615, 72.2003, 76.4784, 81.0101, 85.8102, 90.8948, 96.2806, 101.986, 108.029, 114.43, 121.21, 128.392, 136};
    // SPACING = (b/a)^(1/n) , n= 40, b = 136, a =13.6 
    double fsum=0 ;
    double x0 = u[0];
    
    int len = sizeof(u)/sizeof(u[0]);   
    double df = 1;//0.5*pow(u[len-1]/x0,1.0/(len-1));
   // for(int i=0;i<sizeof(u)/sizeof(u[0]);i++){
 //   fsum  = hintegral(&u[0],len,NH,nh1,nhe1,nhe2);

    for(int i=0;i<40;i++){	
     fsum = fsum + df*(u[i]-x0)*pow(u[i]/x0,alpha-1)*exp(-sigma(u[i])*NH-sigmahe1(u[i])*Nhe1-sigmahe2(u[i])*Nhe2)*PH1(u[i],nh1,nhe1l,nhe2l) ;	
}
    
    return fsum;
}
double Heat2(double NH,double Nhe1, double Nhe2,double nh1,double nhe1l,double nhe2l){
    
    //double u[] ={13.6, 14.4059, 15.2595, 16.1636, 17.1214, 18.1359, 19.2105, 20.3488,21.5545, 22.8317, 24.1846, 25.6176, 27.1356, 28.7435, 30.4466, 32.2507, 34.1617, 36.1859, 38.33, 40.6012, 43.007, 45.5553, 48.2546, 51.1139, 54.1426,57.3507, 60.749, 64.3486, 68.1615, 72.2003, 76.4784, 81.0101, 85.8102, 90.8948, 96.2806, 101.986, 108.029, 114.43, 121.21, 128.392, 136};
    double u[] = {54.5,57.7293, 61.15, 64.7734, 68.6114, 72.6769, 76.9833, 81.5448, 86.3767, 91.4948, 96.9162, 102.659, 108.742, 115.185, 122.01, 129.24, 136.898, 145.01, 153.602, 162.703, 172.344, 182.556, 193.373, 204.831, 216.968};
    // SPACING = (b/a)^(1/n) , n= 40, b = 136, a =13.6 
    double fsum=0 ;
    double x0 = u[0];
    
    int len = sizeof(u)/sizeof(u[0]);
    double df = 1;//0.5*pow(u[len-1]/x0,1.0/(len-1));

    for(int i=0;i<len;i++){
     fsum = fsum + df*(u[i]-x0)*pow(u[i]/x0,alpha-1)*exp(-sigma(u[i])*NH-sigmahe1(u[i])*Nhe1-sigmahe2(u[i])*Nhe2)*PHe1(u[i],nh1,nhe1l,nhe2l) ;   
     
}
 
  /*isothermal case*/
   //  fsum  = df*exp(-sigma(u[0])*NH)*PH1(u[0],nh1,nhe1,nhe2) ;  
   
 //   fsum  = df*exp(-sigma(u[0])*NH-sigmahe1(u[0])*Nhe1-sigmahe2(u[0])*Nhe2)*PHe1(u[0],nh1,nhe1,nhe2) ;
    return fsum;
}

double Heat3(double NH,double Nhe1, double Nhe2,double nh1,double nhe1l,double nhe2l){
    
    //double u[] ={13.6, 14.4059, 15.2595, 16.1636, 17.1214, 18.1359, 19.2105, 20.3488,21.5545, 22.8317, 24.1846, 25.6176, 27.1356, 28.7435, 30.4466, 32.2507, 34.1617, 36.1859, 38.33, 40.6012, 43.007, 45.5553, 48.2546, 51.1139, 54.1426,57.3507, 60.749, 64.3486, 68.1615, 72.2003, 76.4784, 81.0101, 85.8102, 90.8948, 96.2806, 101.986, 108.029, 114.43, 121.21, 128.392, 136};
    double u[] = {54.5,57.7293, 61.15, 64.7734, 68.6114, 72.6769, 76.9833, 81.5448, 86.3767, 91.4948, 96.9162, 102.659, 108.742, 115.185, 122.01, 129.24, 136.898, 145.01, 153.602, 162.703, 172.344, 182.556, 193.373, 204.831, 216.968};
    // SPACING = (b/a)^(1/n) , n= 40, b = 136, a =13.6 
    double fsum=0 ;
    double x0 = u[0];
    
    int len = sizeof(u)/sizeof(u[0]);
    double df = 1;//0.5*pow(u[len-1]/x0,1.0/(len-1));

    for(int i=0;i<len;i++){
     fsum = fsum + df*(u[i]-x0)*pow(u[i]/x0,alpha-1)*exp(-sigma(u[i])*NH-sigmahe1(u[i])*Nhe1-sigmahe2(u[i])*Nhe2)*PHe2(u[i],nh1,nhe1l,nhe2l) ;   
     
}
 
  /*isothermal case*/
   //  fsum  = df*exp(-sigma(u[0])*NH)*PH1(u[0],nh1,nhe1,nhe2) ;  
   
 //   fsum  = df*exp(-sigma(u[0])*NH-sigmahe1(u[0])*Nhe1-sigmahe2(u[0])*Nhe2)*PHe1(u[0],nh1,nhe1,nhe2) ;
    return fsum;
}

double Hubble(double t){
    double  omega_m = 0.32;
    double  H0 = 67.32*1.05e-6;
    double z =  7;//5.945518;
    return sqrt(1-omega_m)*H0*cosh(1.5*H0*sqrt(1-omega_m)*t)*pow(1+z,1.5);///sinh(1.5*H0*sqrt(1-omega_m)*t);
}

double xh2solve(double Ti,double xh2,double xhe1,double xhe2,double xhe3,double dti,double NH1,double Nhe1, double Nhe2,double R,double nhi, double *xrate){
    double nhe1l = (Y/(4*(1-Y)))*nhi*(xhe1);
    double nhe2l = (Y/(4*(1-Y)))*nhi*xhe2;
    double nhe3l = (Y/(4*(1-Y)))*nhi*xhe3;


 //   double nhe2 =0;
 //   double nhe3 =0;
 //   double nhe1 =0;
 //   double xhe2,xhe3;
   
    double ne = xh2*nhi+(Y/(4*(1-Y)))*nhi*(xhe2+2*xhe3);	
//    xhe2 =0;
//    xhe3= 0;
    double fac = -alpha*Ng/pow(3e21,3);
    /*isothermal case */
   // double fac = Ng/pow(3e21,3);
    double vol = 1.0/(4*PI*R*R*cellsize);
    double nh1 = nhi*(1-xh2);
    double C,D,dbdt,dalphadt,nu,Tnu,dCdx,dDdx,J;
    
    C = GE(Ti)*ne + fac*vol*Gamma(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l)/(nh1);
    D = ralpha(Ti)*ne;
    //printf("%lf,%lf\n",C,D);
 //   double X = 1;
 //   double Y = 0;
  /* 
    dbdt = (1.07732*1e14*exp(-157807/Ti))/(pow((1 + 158648*pow((1/Ti),0.874)),1.101) *pow((1/Ti),1.089) *pow(Ti,(7/2)))
- (1.02402*1e9*exp(-157807/Ti))/(pow((1 + 158648*pow((1/Ti),0.874)),1.101)*pow((1/Ti),1.089)*pow(Ti,(5/2))) 
+ (7.43439*1e8*exp(-157807/Ti))/(pow((1 + 158648*pow((1/Ti),0.874)),1.101)*pow((1/Ti),0.089)*pow(Ti,(3/2))) 
+ (1.0422*1e14*exp(-157807/Ti)*pow((1/Ti),0.785))/(pow((1 + 158648*pow((1/Ti),0.874)),2.101)*pow(Ti,(3/2)));
    */
    dbdt = -2.91364*exp(-157809.1/Ti)/pow((1+0.00316228*pow(Ti,0.5)),2) +2.90803*exp(-157809.1/Ti)/((1+0.00316228*pow(Ti,0.5))*pow(Ti,1.5)) + 921.375*exp(-157809.1/Ti)/((1+0.00316228*pow(Ti,0.5))*pow(Ti,0.5));
    
    dalphadt = -((1.10726*1e9*pow((1/Ti),2.503))/(1 + 521.553*pow((1/Ti),0.47))*1.923) + (
 3.4727*1e11*pow((1/Ti),2.973))/pow((1 + 521.553*pow((1/Ti),0.47)),2.923);

    nu = pow((X*(1+xh2)+Y*(1+xhe2+2*xhe3)/4),(-1));
    Tnu = Ti/nu;
    
    dCdx = nhi*GE(Ti) + ne*Tnu*pow(nu,2)*X*dbdt ;
    dDdx = nhi*ralpha(Ti) + ne*Tnu*pow(nu,2)*X*dalphadt;
    J = dCdx -(C+D) - xh2*(dCdx+dDdx);
    *(xrate) = (C-xh2*(C+D));
    
    return xh2 + dti*(C-xh2*(C+D))/(1-J*dti);
}

double xhe1solve(double Ti,double xh2i,double xhe1,double xhe2,double xhe3,double dti,double NH1,double Nhe1,double Nhe2,double R,double nhi, double *hrate){
        double nhe1l = (Y/(4*(1-Y)))*nhi*(xhe1);
        double nhe2l = (Y/(4*(1-Y)))*nhi*xhe2;
        double nhe3l = (Y/(4*(1-Y)))*nhi*xhe3;
   // double xhe1 = 1- xhe2 -xhe3;

        double ne = xh2i*nhi+(Y/(4*(1-Y)))*nhi*(xhe2+2*xhe3); 

        double fac = -alpha*Ng/pow(3e21,3);
    /*isothermal case */
   // double fac = Ng/pow(3e21,3);
        double vol = 1.0/(4*PI*R*R*cellsize);
    
        double nh1 = nhi*(1-xh2i);
        double C,D,nu,Tnu;
        nu = pow((X*(1+xh2i)+Y*(1+xhe2+2*xhe3)/4),(-1));
        Tnu = Ti/nu;
    
        C = xhe2*alphahe2(Ti)*ne;//+ xhe3*alphahe3(Ti)*ne;
   // D = betahe2(Ti)*ne + alphahe2(Ti)*ne + fac*vol*Gammahe2(NH1,Nhe1,Nhe2,nh1,nhe1,nhe2)/(nhe2) ;
        D = betahe1(Ti)*ne + fac*vol*Gammahe1(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l)/(nhe1l);
  // printf("C %e D %e\n",C,D );
        *(hrate) = (C-D*xhe1);
   return (xhe1 + C*dti)/(1+D*dti);
   //     return xhe1 + (C-D*xhe1)*dti;
}


double xhe2solve(double Ti,double xh2i,double xhe1,double xhe2,double xhe3,double dti,double NH1,double Nhe1,double Nhe2,double R,double nhi,double *hrate){
    double nhe1l = (Y/(4*(1-Y)))*nhi*(xhe1);
    double nhe2l = (Y/(4*(1-Y)))*nhi*xhe2;
    double nhe3l = (Y/(4*(1-Y)))*nhi*xhe3;
   // double xhe1 = 1- xhe2 -xhe3;

    double ne = xh2i*nhi+(Y/(4*(1-Y)))*nhi*(xhe2+2*xhe3); 

    double fac = -alpha*Ng/pow(3e21,3);
    /*isothermal case */
   // double fac = Ng/pow(3e21,3);
    double vol = 1.0/(4*PI*R*R*cellsize);
    
    double nh1 = nhi*(1-xh2i);
    double C,D,nu,Tnu;
    nu = pow((X*(1+xh2i)+Y*(1+xhe2+2*xhe3)/4),(-1));
    Tnu = Ti/nu;
    
    C = xhe1*ne*betahe1(Ti) + xhe1*fac*vol*Gammahe1(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l)/(nhe1l) + xhe3*alphahe3(Ti)*ne;
    D = betahe2(Ti)*ne + alphahe2(Ti)*ne + fac*vol*Gammahe2(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l)/(nhe2l) ;
   // D = alphahe2(Ti)*ne;
  // printf("C %e D %e\n",C,D );
    *(hrate) = (C-D*xhe2);
    return (xhe2 + C*dti)/(1+D*dti);
   // return xhe2 + (C-D*xhe2)*dti;
}

double xhe3solve(double Ti,double xh2i,double xhe1,double xhe2,double xhe3,double dti,double NH1,double Nhe1,double Nhe2,double R,double nhi,double *hrate){
    double nhe1l = (Y/(4*(1-Y)))*nhi*(xhe1);
    double nhe2l = (Y/(4*(1-Y)))*nhi*xhe2;
    double nhe3l = (Y/(4*(1-Y)))*nhi*xhe3;
    //double xhe1 = 1- xhe2 -xhe3;

    double ne = xh2i*nhi+(Y/(4*(1-Y)))*nhi*(xhe2+2*xhe3);    
    double fac = -alpha*Ng/pow(3e21,3);
    /*isothermal case */
   // double fac = Ng/pow(3e21,3);
    double vol = 1.0/(4*PI*R*R*cellsize);
    double nh1 = nhi*(1-xh2i);
    double C,D,nu,Tnu;
    nu = pow((X*(1+xh2i)+Y*(1+xhe2+2*xhe3)/4),(-1));
    Tnu = Ti/nu;
    
    C = xhe2*ne*betahe2(Ti) + fac*vol*Gammahe2(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l)/((Y/(4*(1-Y)))*nhi);
    D = alphahe3(Ti)*ne;
  //  printf("C %e D %e\n",C,D );
    *(hrate) = (C-D*xhe3);
   return (xhe3 + C*dti)/(1+D*dti);
   // return xhe3 + (C-D*xhe3)*dti;
}

double Tsolve(double Ti,double NH1,double Nhe1,double Nhe2,double xh2i,double xhe1,double xhe2,double xhe3,double nhi,double dti,double R,double tsim, double *trate){
    
    
    double fac = -alpha*Ng/pow(3e21,3);
    double fh = fac*1.6e-12; // 1.6e-12 is factor ev to erg
    double vol = 1.0/(4*PI*R*R*cellsize);
    
    double K = (2.0/3)*(1.0/(nhi*kb)); // only H # cm^3erg^-1 K    
    
    double nh1 = nhi*(1-xh2i);
    double nhe1l = (Y/(4*(1-Y)))*nhi*(xhe1);
    double nhe2l = (Y/(4*(1-Y)))*nhi*xhe2;
    double nhe3l = (Y/(4*(1-Y)))*nhi*xhe3;
    double ne = xh2i*nhi+(Y/(4*(1-Y)))*nhi*(xhe2+2*xhe3);

    double T5 = Ti/100000.0;
    
    double f1 = 1.778e-29*3.15e13;	//#*6.24e11 # in erg cm^3 Myr^-1 K^-1 #4.05   #
    double f2 = 7.5e-19*3.15e13;	//#*6.24e11 # in erg cm^3 Myr^-1 #1.711e+11
    double f3 = 1.42e-27*3.15e13;	//#*6.24e11# in erg cm^3 K^-1 #489.62
    double gff = 1.5 ;	//# from cen #1.1+ 0.34*(np.exp(-5.5+np.log(T)))**2
    double f4 = 5.65e-36*3.15e13;	//#*6.24e11 #in erg cm^3 Myr^-1 K^-1 #1.289e-6

    double f11 = 1.0;
    double f12 = 8*1.778e-29*3.15e13;
    double f21 = 3.15e13*5.54e-17;
/*
    f11 = 1
    f12 =8*1.778*1e-29*3.15e13/(1.38e-16)
    f21 =5.54e-17*3.15e13/(1.38e-16)
    f22 =9.1e-27*3.15e13/(1.38e-16)
*/
    double Th1 = 157807;
    //double The = 631515;
    double lamb = 2*Th1/Ti;
    double lamb2 = 1263030.0/Ti;
    double omega_m = 0.32;
    double H0 = 67.32*1.05e-6;
    double rho = 1.0;
    //#fh = 1.04e+8
    double Hf = 1.05e-6;
    
    double z = 7;//5.945518;
    
  //  double X = 1;
  //  double Y = 0;
  //  double xhe2 = 0; 
  //  double xhe3 = 0;
    double nu = pow((X*(1+xh2i)+Y*(1+xhe2+2*xhe3)/4),(-1));



    
    double c1 = (f1)*Ti*pow(lamb,1.965)/pow((1.0+pow((lamb/0.541),0.502)),2.697) ;// recombination cooling coeffecient erg cm^3 Myr^-1
    double c2 = (f2)*pow((1+pow(T5,1.0/2)),-1)*exp(-118355.0/Ti);  //collisional excitation cooling coeffecient
    double c3 = gff*(f3)*pow(Ti,0.5); // Brehmstrallung/ free free emission
    double c4 = (f4)*(Ti-2.73*(1+z))*pow(1+z,4);
    
    double c11 = f11*kb*Ti*alphahe2(Ti); // recombination cooling He II erg cm^3 Myr^-1
    double c12 = f12*Ti*pow(lamb2,1.965)/pow((1.0+pow((lamb2/0.541),0.502)),2.697) ; // recombination cooling He III
    
    double c21 = f21*pow(Ti,-0.397)*pow((1+pow(T5,1.0/2)),-1)*exp(-473638.0/Ti);    //collisional excitation HeII erg cm^3 Myr^-1

    //double C =  c1*pow((xh2i*nhi),2) +c2*xh2i*nhi*(nh1)+c3*pow((nhi*xh2i),2) +c4*(nhi*xh2i) + 2*Hubble(tsim)*Ti/(K*nu) ;//+ (Ti/(nhi*xh2i))*(*xrate);
    
    double C; 

  //  C = c1*ne*xh2i*nhi + c2*ne*nh1 + c3*(nhi*xh2i+nhe2+4*nhe3)*ne + c4*ne+ 2*Hubble(tsim)*Ti/(K*nu);

    C = c1*ne*xh2i*nhi + c11*ne*nhe2l + c12*ne*nhe3l +c2*ne*nh1 + c21*ne*nhe2l + c3*(nhi*xh2i+nhe2l+4*nhe3l)*ne + c4*ne + 2*Hubble(tsim)*Ti/(K*nu);
    //double C =0;
       
    double L = (fh*vol*Heat(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l)+fh*vol*Heat2(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l)+fh*vol*Heat3(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l) - C);
    
     //printf("%e, %e\n",c1, c2);
    //derivatives of cooling processes 
    
    double A1 =  -(6.1712*1e10*f1*pow((1/Ti),1.965))/pow((1 + 784.353*pow((1/Ti),0.502)),2.697) + (
 6.79107*1e13*f1*pow((1/Ti),2.467))/pow((1 + 784.353*pow((1/Ti),0.502)),3.697); 
    double A2 = (118355*exp(-118355/Ti)*f2)/((1 + 0.00316228*pow(Ti,0.5))*pow(Ti,2)) - (
 0.00158114*exp(-118355/Ti)*f2)/(pow((1 + 0.00316228*pow(Ti,0.5)),2)*pow(Ti,0.5)); 
    double A3 = (0.5*f3*gff)/pow(Ti,0.5); 
    double A4 =   f4*pow((1+z),4);
    
    double A11 = -4.96095e-13*pow(1/Ti,1.654);
    double A12 = -(9.41467e11*f12*pow(1/Ti,1.965)/(pow((1 + 1573.42*pow((1/Ti),0.502)),2.697)))+ 2.07829e15*f12*pow(1/Ti,2.467)/pow((1 + 1573.42*pow((1/Ti),0.502)),3.697);
    double A21 = 473638*exp(-473638/Ti)*f21/((1+pow(T5,1.0/2))*pow(Ti,2.397))- 0.397*exp(-473638/Ti)*f21/((1+pow(T5,1.0/2))*pow(Ti,1.397)) - exp(-473638/Ti)*f21/(200*pow(10,0.5)*pow((1+pow(T5,1.0/2)),2)*pow(Ti,0.897));
    
    double dLdt = -(A1+ A2 +A3 +A4 +A11+ A12+ A21 +2*Hubble(tsim)/(K*nu));//+(1/(nhi*xh2i))*(*xrate));
    //double dLdt =0;
   // printf("%e\n", K)
	//printf("%e %e\n",fh*vol*Heat(NH1,nh1,nhe1,nhe2)-C, dti);
    *(trate) = L*K*nu;
   return Ti+ dti*L*K*nu;
   //return Ti + nu*dti*L*K/(1-nu*dLdt*K*dti) ;
   //return (L*K*nu*dti + Ti)/(1+2*Hubble(tsim)*Ti*dti) ; // anninos implicit
}


int xTiter(double dti,double NH1,double Nhe1, double Nhe2,double R,double nhi,double tsim,double ti,double *Tt, double *xh2t,double *xhe1t,double *xhe2t,double *xhe3t,double *dtrec){
  //  double *Tt =(double*)malloc(countl*sizeof(double));
  //  Tt[0] = T[k][j];
  //  double *xh2t =(double*)malloc(countl*sizeof(double)); 	
   // xh2t[0] = nh2[k][j]/nhi;
   // double *NHt =(double*)malloc(countl*sizeof(double)); 	
   // NHt[0] = NH[k][j];
   // printf("%d %lf\n",k,Tt[0]);
    int storei;
    double frac_change = 0.0;
    //double ne_init;
    double dUU;
    double Trate=0;
    double Xrate=0;
    double he1rate = 0;
    double he2rate = 0;
    double he3rate = 0;
    for (int i=0;i<countl;i++){
        

        //Tt[i+1] = (double) max(T2_MIN_FIX,Tsolve(Tt[i],NH1,Nhe1,Nhe2,xh2t[i],xhe1t[i],xhe2t[i],xhe3t[i],nhi,dti,R,tsim,&Trate));
        Tt[i+1] = Tsolve(Tt[i],NH1,Nhe1,Nhe2,xh2t[i],xhe1t[i],xhe2t[i],xhe3t[i],nhi,dti,R,tsim,&Trate);
      //  printf("T %e\n",Tt[i+1] );
        dUU = (double) fabs(max(T2_MIN_FIX,Tt[i]+Trate*dti)-Tt[i]);
     //   printf("heat %e\n",Heat(NH1,nhi,0,0)*(-alpha*Ng/pow(3e21,3))* 1.6e-12*1.0/(4*PI*R*R*cellsize*kb));
 //   printf("%e\n",fabs(Tt[i+1]-Tt[i])/(Tt[i]+T_MIN));
 //   printf("%e\n",Tt[i+1]);

//if(dti>1e-7){
      //  printf("i %d\n", i);
   //     printf("change %e\n",(double) max(dUU,fabs(Tt[i+1]-Tt[i]))/(Tt[i]+T_MIN) );


	if((double) max(dUU,fabs(Tt[i+1]-Tt[i]))/(Tt[i]+T_MIN)>0.1){
	//printf("hii");
	return 0;
}
	else{
	//	printf("hi");
	//frac_change = (double) max(frac_change,fabs(Tt[i+1]-Tt[i])/Tt[i]);
        frac_change = (double) max(frac_change,dUU);

}

//}

xh2t[i+1] = xh2solve(Tt[i+1],xh2t[i],xhe1t[i],xhe2t[i],xhe3t[i],dti,NH1,Nhe1,Nhe2,R,nhi,&Xrate);
//printf("xh2 %e\n",xh2t[i+1] );
/*complete this exactly*/

     if((double) max(fabs(xh2t[i+1]-xh2t[i]),Xrate*dti)/(xh2t[i]+X_FM) >0.1){
	return 1;
}
else{
	frac_change = (double) max(frac_change,fabs(xh2t[i+1]-xh2t[i])/xh2t[i]);
}

//xhe1t[i+1] = xhe1solve(Tt[i+1],xh2t[i],xhe2t[i],xhe3t[i],dti,NH1,Nhe1,Nhe2,R,nhi);
xhe1t[i+1] = xhe1solve(Tt[i+1],xh2t[i+1],xhe1t[i],xhe2t[i],xhe3t[i],dti,NH1,Nhe1,Nhe2,R,nhi,&he1rate);
/*
if((double) max(fabs(xhe1t[i+1]-xhe1t[i]),he1rate*dti)/(xhe1t[i]+X_FM) >0.1){
    return 2;
}
else{
    frac_change = (double) max(frac_change,fabs(xhe1t[i+1]-xhe1t[i])/xhe1t[i]);
}
*/
xhe2t[i+1] = xhe2solve(Tt[i+1],xh2t[i+1],xhe1t[i+1],xhe2t[i],xhe3t[i],dti,NH1,Nhe1,Nhe2,R,nhi,&he2rate);
/*
if((double) max(fabs(xhe2t[i+1]-xhe2t[i]),he2rate*dti)/(xhe2t[i]+X_FM) >0.1){
    return 3;
}
else{
    frac_change = (double) max(frac_change,fabs(xhe2t[i+1]-xhe2t[i])/xhe2t[i]);
}
*/
xhe3t[i+1] = xhe3solve(Tt[i+1],xh2t[i+1],xhe1t[i+1],xhe2t[i+1],xhe3t[i],dti,NH1,Nhe1,Nhe2,R,nhi,&he3rate);
/*
if((double) max(fabs(xhe3t[i+1]-xhe3t[i]),he3rate*dti)/(xhe3t[i]+X_FM) >0.1){
    return 4;
}
else{
    frac_change = (double) max(frac_change,fabs(xhe3t[i+1]-xhe3t[i])/xhe3t[i]);
}
//}
*/
/*
	ne_init = xh2t[i]*nhi;
	if(fabs(xh2t[i+1]*nhi-ne_init)/ne_init>0.1){
	return 2;
}
else {
	frac_change = (double) max(frac_change,fabs(xh2t[i+1]*nhi-ne_init)/(ne_init+X_FM));
}
*/	
	double tcal = ti + dti;
    ti = tcal ;

//if (dti>1e-7){
  

    if(frac_change <0.05){  // has value 0.05 in paper but 0.5 here
    *(dtrec)  = 2*dti;
    //    dti = 2*dti;
///}
/*
else{
    *(dtrec) = dti;
    printf("dti %e\n",dti );
}
*/
}
else{
    *(dtrec) = dti;
  //  dti = dti;
  //  printf("dti %e\n",dti );
}
   
   // printf("dti %e\n", dti);
	
    
        if (tcal >= tsim){
            storei = i+1;
            printf("tsim %lf\n", tcal);
            break ;
        }
 dti = (double) min(*(dtrec),tsim-tcal);
  //printf("i %d dti %e\n", i, dti);

   }

	printf("%lf %d\n",Tt[storei],storei);
//	T[k][j+1] = Tt[storei];
//	nh2[k][j+1] = xh2t[storei]*nhi;

    return storei;//Tt[storei];
}


double evolve_eqns(double *nh,double **nh2,double **NH,double **nhe1,double **nhe2,double **nhe3,double **NHe1,double **NHe2,double **T,int ncells,double ti,double tsim,double **dt, int j,double Dt){
	//double dt = Dt;
	//double X =1 ;
	//double Y =0;
	//double xhe2 =0;
	//double xhe3 =0;
	int out;
	double *Tt =(double*)malloc(countl*sizeof(double));
	double *xh2t =(double*)malloc(countl*sizeof(double));
    double *xhe1t =(double*)malloc(countl*sizeof(double));
    double *xhe2t =(double*)malloc(countl*sizeof(double));
    double *xhe3t =(double*)malloc(countl*sizeof(double));
	double dt_rec;
	for(int k=0; k<ncells;k++){
    printf("dt new %e %d %d\n",dt[k][j],k,j);
    //dt[k][j+1] = dt[k][j];
    printf("cell %d\n",k);
	//mu = pow((X*(1+nh2[k][j]/nh[k])+Y*(1+nh2[k][j]/nh[k]+xhe2+2*xhe3)/4),(-1));
    Tt[0] = T[k][j];
    xh2t[0] = nh2[k][j]/nh[k];
    xhe2t[0] = nhe2[k][j]/((Y/(4*(1-Y)))*nh[k]);
    xhe3t[0] = nhe3[k][j]/((Y/(4*(1-Y)))*nh[k]);
    xhe1t[0] = 1-xhe2t[0] -xhe3t[0];
	out = xTiter(dt[k][j],NH[k][j],NHe1[k][j],NHe2[k][j],2+(k)*cellsize,nh[k],tsim,ti,Tt,xh2t,xhe1t,xhe2t,xhe3t,&dt_rec);  // equivalent of cool_step routine
	//printf("%d\n",out);
	

    while((out==0 || out==1 || out ==2 || out ==3 || out ==4)){
	dt[k][j] =dt[k][j]/2;
	printf("%d %e\n",out, dt[k][j]);
    /*if (dt[k][j]<1e-8)
    {
        out = xTiter(dt[k][j],NH[k][j],(1+k)*cellsize,nh[k],tsim,ti,Tt,xh2t,&dt_rec);
        break;
    }*/
	out = xTiter(dt[k][j],NH[k][j],NHe1[k][j],NHe2[k][j],2+(k)*cellsize,nh[k],tsim,ti,Tt,xh2t,xhe1t,xhe2t,xhe3t,&dt_rec);

	//printf("%d\n",out);
} 
	printf("%d\n",out);
	T[k][j+1] = Tt[out];
	nh2[k][j+1] = nh[k]*xh2t[out];
    nhe2[k][j+1] = (Y/(4*(1-Y)))*nh[k]*xhe2t[out];
    nhe3[k][j+1] = (Y/(4*(1-Y)))*nh[k]*xhe3t[out];
    nhe1[k][j+1] = (Y/(4*(1-Y)))*nh[k]*xhe1t[out];//(Y/(4*(1-Y)))*nh[k] - nhe2[k][j+1] - nhe3[k][j+1];
 	//nh2[k][j+1] = nh[k]*xh2iter(T[k][j],nh2[k][j],nh2[k][j]/nh[k],dt,NH[k][j],(1+k)*cellsize,nh[k],tsim,ti);	
	//printf("%e\n",nh2[k][j+1]);	
	//*(*(nh2+k)+j+1) = 0.1*(*(*(nh2+k)+j));
	//T[k][j+1] = T[k][j];
//	printf("dtrec %e\n",dt_rec );

	dt[k][j+1] = Dt;//(double) min((dt_rec),tsim-ti); // check this criterion
    

    //dt[k][j+1] =(double) max(dt_rec,1e-8);
//    printf("%e\n",(double) min((dt_rec),tsim-ti-dt[k][j]));
  //  printf("dts %e %e %d\n",(dt_rec), tsim-ti-dt[0][j],j);
  //  printf("dtrec %e %d %d\n", dt[k][j+1],k,j+1);
}
	for(int k=0;k<ncells;k++){
	NH[k][j+1] = cumsum(nh,nh2,k,j+1)*cellsize*3e7; // in cm^-2
    NHe1[k][j+1] = hecumsum(nhe1,k,j+1)*cellsize*3e7 ;
    NHe2[k][j+1] = hecumsum(nhe2,k,j+1)*cellsize*3e7 ;
}
		return dt_rec;
}


int main(int argc,char **argv)
{

	/*For each position, solve euler to get nh(xpos,t) and T(xpos,t) where t is tsim= time after each global dt */	
	clock_t t;
    	t = clock();
    	
    	//printf("%lf\n",(double) max(3e-3,4.4e-3));
	
	double ti = 935.7;
	double tf = ti + atof(argv[1]);
//	printf("%lf\n",tf);
	double epsilon = 0.05;	
	double Dt = 3261.6*(epsilon/0.1)*(cellsize/10.0)*1e-6 ;
	//double Dt = 4500e-6;	
	printf("%e\n",Dt);
	int ncells = boxsize/cellsize;	
/* setting up grid */
	double *x = (double*)malloc(ncells*sizeof(double));
	for (int i = 0; i < ncells; ++i){
	x[i] = i*cellsize + cellsize;
	//printf("%lf\n",x[i]);
}

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

/* Initial values*/
	for (int i = 0; i < ncells; ++i){
	//printf("%d\n",i);
	nh[i] = 7e-2;
	nh2[i][0] = 0.05*nh[i];
    nhe2[i][0] = 1e-9*(Y/(4*(1-Y)))*nh[i];
    nhe3[i][0] = 0;
    nhe1[i][0] = (Y/(4*(1-Y)))*nh[i] - nhe2[i][0] - nhe3[i][0];
	T[i][0] = 10;
	NH[i][0] = cumsum(nh,nh2,i,0)*cellsize*3e7; // in cm^-2
    NHe1[i][0] = hecumsum(nhe1,i,0)*cellsize*3e7 ;
    NHe2[i][0] = hecumsum(nhe2,i,0)*cellsize*3e7 ;
    dt[i][0] = Dt;	
	//printf("%lf\n",NHe2[i][0]);
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
//double dt;


int storej;
double dt_new;
for (int j=0;j<count;j++){
	tsim = ti + Dt;
	//printf("%d\n",j);	
 //   printf("dt new out %e \n",dt[0][j]);
    /*evolve all cells for a  global timestep*/
    evolve_eqns(nh,nh2,NH,nhe1,nhe2,nhe3,NHe1,NHe2,T,ncells,ti,tsim,dt,j,Dt);	// my equivalent of subroutine rt_solve_cooling
//	printf("dt new out %e \n",dt[0][j+1]);

    ti = tsim;
	printf("%lf\n",tsim);
	if(fabs(tf-tsim)<=Dt){
	printf("%lf,%lf\n",tf,tsim);
	storej = j+1;
	break;
} 
    
}

printf("[");
//printf("%d\n",storej);

for (int i=0;i<ncells;i++){

printf("%e,",1-nh2[i][storej]/nh[i]);
}
printf("]\n");

for (int i=0;i<ncells;i++){

printf("%e,",T[i][storej]);
}
printf("]");

t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
  
    printf("fun() took %f seconds to execute \n", time_taken);	    


 FILE *fptr = fopen("xh1_1myr.txt","w");
 FILE *fdptr = fopen("xhe2_1myr.txt","w");
FILE *fd1ptr = fopen("xhe3_1myr.txt","w");
FILE *fd2ptr = fopen("T_1myr.txt","w");
FILE *fd3ptr = fopen("xhe1_1myr.txt","w");
/*
FILE *fptr = fopen("intblum_alphacrct_xh1_001myr.txt","w");
FILE *fdptr = fopen("intblum_alphacrct__T_001myr.txt","w");    
*/
     	for (int l = 0; l < ncells; ++l)
    {
    fprintf(fptr,"%lf %0.12e \n", (double)l*cellsize+cellsize,(nh[l]-nh2[l][storej])/nh[l]);
	fprintf(fd2ptr,"%0.12e \n",T[l][storej]);
    fprintf(fd3ptr,"%0.12e \n", nhe1[l][storej]/((Y/(4*(1-Y))*nh[l])));
    fprintf(fdptr,"%0.12e \n",nhe2[l][storej]/((Y/(4*(1-Y))*nh[l])));
    fprintf(fd1ptr,"%0.12e \n",nhe3[l][storej]/((Y/(4*(1-Y))*nh[l])));
    } 
	

}

