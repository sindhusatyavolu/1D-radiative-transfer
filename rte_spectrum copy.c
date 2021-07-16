#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h> 
#include<time.h>


#ifndef PI
#define PI 3.1415
#endif
#define cellsize 55.46875//78.125
#define boxsize 113600.0//160000.0
#define count 10000000
#define countl 1000000000
#define Ng 3.15e13*1e57
#define alpha -1.5
#define kb 1.38e-16
#define max(x,y) (((x) >= (y)) ? (x) : (y))
#define min(x,y) (((x) <= (y)) ? (x) : (y))
#define T_MIN 0.1
#define T2_MIN_FIX 1e-2
#define X_MIN  1e-20
#define X_FM 1e-6
#define X  0.742
#define Y 0.258
#define lim 2e-1
#define epsilon 5e-3
#define x_min_fix 1e-20
#define dt_min 1e-8
#define z 5.945518
#define background 0.0
#define HUBBLE 3.2407789e-18
#define  BOLTZMANN 1.3806e-16
#define  PROTONMASS  1.6726e-24
#define  GAMMA_LYA_H1   6.265e8  /* s^-1 */
#define  LYMAN_ALPHA      1215.6e-8      /* 1215.6 Angstroem */
#define  OSCILLATOR_STRENGTH       0.41615
#define  THOMPSON     6.65245e-25
#define  Cv           2.9979e10


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
   // return 3.15e13*5.85e-11*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-157809.1/Ti);
    //return 3.15e13*1.17e-10*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-157809.1/Ti);//bolton
   return  3.15e13*21.11*pow(Ti,-1.5)*exp(-lam/2)*pow(lam,-1.089)/pow(1+pow((lam/0.354),0.874),1.101) ;
}

double betahe1(double Ti){
    double lam  = 2*285335/Ti;
    return 3.15e13*32.38*pow(Ti,-1.5)*exp(-lam/2)*pow(lam,-1.146)/pow(1+pow((lam/0.416),0.987),1.056) ; // from HnG
   // double T5 =Ti/100000 ;
    
   // return 3.15e13*2.38e-11*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-285335.4/Ti); //ramses
   // return 3.15e13*4.76e-11*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-285335.4/Ti);//bolton
}

double betahe2(double Ti){
    double lam = 2*631515/Ti;
    return 3.15e13*19.95*pow(Ti,-1.5)*exp(-lam/2)*pow(lam,-1.089)/pow(1+pow((lam/0.553),0.735),1.275) ;
    //double T5 =Ti/100000 ;
    //return 3.15e13*5.68e-12*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-631515.0/Ti); //ramses
    //return 3.15e13*1.14e-11*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-631515.0/Ti); //bolton
}


double ralpha(double Ti){
    double Th1 = 157807; 
    double lam = 2*Th1/Ti;

   return 3.15e13*1.269e-13*pow(lam,1.503)/pow((1+pow(lam/0.522,0.47)),1.923);
  //  double Tl = kb*Ti/1.6e-12;
  //  return 3.15e13*(exp(-28.6130338-0.72411256*log(Tl)-2.02604473e-2*log(pow(Tl,2))-2.38086188e-3*log(pow(Tl,3))-3.21260521e-4*log(pow(Tl,4))-1.42150291e-5*log(pow(Tl,5))+4.98910892e-6*log(pow(Tl,6))+5.75561414e-7*log(pow(Tl,7))-1.85676704e-8*log(pow(Tl,8))-3.07113524e-9*log(pow(Tl,9))));

}

double alphahe2(double Ti){

    double lam = 570670/Ti;
    double Tl = kb*Ti/1.6e-12;
   return 3.15e13*3e-14*pow(lam,0.654);
//   return 3.15e13*(3.925e-13*pow(Tl,-0.6353)+1.544e-9*pow(Tl,-1.5)*exp(-48.596/Tl)*(0.3+exp(8.1/Tl)));
}


double alphahe3(double Ti){
  
    double lam = 1263030/Ti;

    return 3.15e13*2.538e-13*pow(lam,1.503)/pow((1+pow(lam/0.522,0.47)),1.923);
  //  return 2*ralpha(Ti/4);//bolton
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

double PH1(double x,double nh1,double nhe1l,double nhe2l){
	
    double qH = exp(-sigma(x)*nh1*cellsize*3e7);
	double qHe1 = exp(-sigmahe1(x)*nhe1l*cellsize*3e7);
    double 	qHe2 = exp(-sigmahe2(x)*nhe2l*cellsize*3e7)	;
    double	tau_tot = sigma(x)*nh1*cellsize*3e7+sigmahe1(x)*nhe1l*cellsize*3e7+ sigmahe2(x)*nhe2l*cellsize*3e7;
    double	D = (1-qH)*qHe2*qHe1 +(1-qHe1)*qH*qHe2 + (1-qHe2)*qH*qHe1;
    //printf("%0.12e\n",D);
    if (D ==0 || D!= D)
    {
        printf("D %0.12e\n",D);
        printf("%0.12e %0.12e\n",nhe1l,nhe2l );
        exit(0);
    }
    
	return (1-qH)*qHe2*qHe1*(1-exp(-tau_tot))/D;
}

double PHe1(double x,double nh1,double nhe1l,double nhe2l){
    double qH = exp(-sigma(x)*nh1*cellsize*3e7);
    double qHe1 = exp(-sigmahe1(x)*nhe1l*cellsize*3e7);
    double  qHe2 = exp(-sigmahe2(x)*nhe2l*cellsize*3e7)  ;
    double  tau_tot = sigma(x)*nh1*cellsize*3e7+sigmahe1(x)*nhe1l*cellsize*3e7+ sigmahe2(x)*nhe2l*cellsize*3e7;
    double  D = (1-qH)*qHe2*qHe1 +(1-qHe1)*qH*qHe2 + (1-qHe2)*qH*qHe1;
  // if ((qH)*qHe2*(1-qHe1)*(1-exp(-tau_tot))/D==0)
  //  {
  //      printf("%0.12e,%0.12e,%0.12e,%0.12e\n",qH,qHe1,qHe2, 1-exp(-tau_tot));
   // }
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
    
    //double u[] ={13.6, 14.4059, 15.2595, 16.1636, 17.1214, 18.1359, 19.2105, 20.3488,21.5545, 22.8317, 24.1846, 25.6176, 27.1356, 28.7435, 30.4466, 32.2507, 34.1617, 36.1859, 38.33, 40.6012, 43.007, 45.5553, 48.2546, 51.1139, 54.1426,57.3507, 60.749, 64.3486, 68.1615, 72.2003, 76.4784, 81.0101, 85.8102, 90.8948, 96.2806, 101.986, 108.029, 114.43, 121.21, 128.392, 136};
    // SPACING = (b/a)^(1/n) , n= 40, b = 136, a =13.6 
    //double u[]={13.6,15.2594509785,17.1213856004,19.2105106069,21.5545474175,24.1845999765,27.1355674836,30.4466074845,34.1616554685,38.3300078652,43.0069761783,48.2546209358,54.1425751953,60.7489685325,68.1614637733,76.4784202259,85.8101988493,96.2806266762,108.028639923,121.210127586,136.0};
  double u[]={13.6,15.3521953269,17.3301398056,19.5629184807,22.083363653,24.9285376673,28.1402779032,31.7658119796,35.8584522227,40.4783796062,45.6935286935,51.580586594,58.2261206226,65.7278512446,74.1960890239,83.755356705,94.5462202801,106.727355968,120.477883496,136.0};
    double fsum =0;
    double x0 = u[0];
    //printf("U %0.12e\n",u[0] );
    int len = sizeof(u)/sizeof(u[0]);
  //  double df = 1;//0.5*pow(u[len-1]/x0,1.0/(len));
   // fsum = integral(&u[0],len,NH,nh1,nhe1,nhe2);
    //for(int i=0;i<sizeof(u)/sizeof(u[0]);i++){
    //printf("%0.12e df\n",df );
  // printf("len %d\n",len );
    double lent =len;
    //fsum = 0.5*(u[1]-u[0])*pow(u[0]/x0,alpha-1)*exp(-sigma(u[0])*NH-sigmahe1(u[0])*Nhe1-sigmahe2(u[0])*Nhe2)*PH1(u[0],nh1,nhe1l,nhe2l)+0.5*(u[19]-u[18])*pow(u[19]/x0,alpha-1)*exp(-sigma(u[19])*NH-sigmahe1(u[19])*Nhe1-sigmahe2(u[19])*Nhe2)*PH1(u[19],nh1,nhe1l,nhe2l);
   
    for(int i=0;i<len-1;i++){
    // double df = 0.5*(u[i+1]-u[i-1]);//0.5*x0*pow(u[len-1]/x0,i/(lent))*(u[len-1]/x0-1);
    // double x1 = u[i] ;//(u[i+1]+u[i])/2.0;
    // double x2 = u[i+1];
     double df = (u[i+1]-u[i]);
    // fsum = fsum + 0.5*df*(pow(x1/x0,alpha-1)*exp(-sigma(x1)*NH-sigmahe1(x1)*Nhe1-sigmahe2(x1)*Nhe2)*PH1(x1,nh1,nhe1l,nhe2l)+pow(x2/x0,alpha-1)*exp(-sigma(x2)*NH-sigmahe1(x2)*Nhe1-sigmahe2(x2)*Nhe2)*PH1(x2,nh1,nhe1l,nhe2l)) ;	
     fsum = fsum + df*pow(u[i]/x0,alpha-1)*exp(-sigma(u[i])*NH-sigmahe1(u[i])*Nhe1-sigmahe2(u[i])*Nhe2)*PH1(u[i],nh1,nhe1l,nhe2l);
}

// fsum  = fsum  + (u[len-1]-u[len-2])*pow(u[len-1]/x0,alpha-1)*exp(-sigma(u[len-1])*NH-sigmahe1(u[len-1])*Nhe1-sigmahe2(u[len-1])*Nhe2)*PH1(u[len-1],nh1,nhe1l,nhe2l) ;  
     
  /*isothermal case*/
   //  fsum  = df*exp(-sigma(u[0])*NH)*PH1(u[0],nh1,nhe1,nhe2) ;  
   //  fsum  = df*exp(-sigma(u[0])*NH-sigmahe1(u[0])*Nhe1-sigmahe2(u[0])*Nhe2)*PH1(u[0],nh1,nhe1l,nhe2l) ;
    //fsum  = df*exp(-sigma(u[0])*NH)*PH1(u[0],nh1,nhe1l,nhe2l) ;
    //printf("%0.12e\n",nhe2l );
    return fsum;
}

double Gammahe1(double NH,double Nhe1, double Nhe2,double nh1,double nhe1l,double nhe2l){
    
    //double u[] ={13.6, 14.4059, 15.2595, 16.1636, 17.1214, 18.1359, 19.2105, 20.3488,21.5545, 22.8317, 24.1846, 25.6176, 27.1356, 28.7435, 30.4466, 32.2507, 34.1617, 36.1859, 38.33, 40.6012, 43.007, 45.5553, 48.2546, 51.1139, 54.1426,57.3507, 60.749, 64.3486, 68.1615, 72.2003, 76.4784, 81.0101, 85.8102, 90.8948, 96.2806, 101.986, 108.029, 114.43, 121.21, 128.392, 136};
    //double u[] = {54.5,57.7293, 61.15, 64.7734, 68.6114, 72.6769, 76.9833, 81.5448, 86.3767, 91.4948, 96.9162, 102.659, 108.742, 115.185, 122.01, 129.24, 136.898, 145.01, 153.602, 162.703, 172.344, 182.556, 193.373, 204.831, 216.968};
    // SPACING = (b/a)^(1/n) , n= 40, b = 136, a =13.6 
    //double u[]={54.5,61.1500057595,68.6114349428,76.9832961819,86.3766789891,96.9162278471,108.741796166,122.010302052,136.897810517,153.601869754,172.344132479,193.373297132,216.968407952,243.442557722,273.147042327,306.476022229,343.871752742,385.830452489,432.908887925,485.731761283,545.0};
    double u[]={24.6,27.6016539758,30.9695651301,34.7484235977,38.9883725345,43.745673487,49.0834529482,55.0725400088,61.7924062151,69.3322201091,77.7920304401,87.2840937515,97.9343639562,109.884163669,123.292059472,138.335965997,155.215506742,174.154662958,195.404745742,219.247730781,246.0};
    double fsum=0 ;
    double x0 = 13.6 ;//u[0];
    
    int len = sizeof(u)/sizeof(u[0]);
    //double df = 1;//0.5*pow(u[len-1]/x0,1.0/(len-1));
   // fsum = integral(&u[0],len,NH,nh1,nhe1,nhe2);
    //for(int i=0;i<sizeof(u)/sizeof(u[0]);i++){
 
for(int i=0;i<len-1;i++){
    // double df = 0.5*(u[i+1]-u[i-1]);//0.5*x0*pow(u[len-1]/x0,i/(lent))*(u[len-1]/x0-1);
     double x1 = u[i] ;//(u[i+1]+u[i])/2.0;
     double x2 = u[i+1];
     double df = (u[i+1]-u[i]);
    // fsum = fsum + 0.5*df*(pow(x1/x0,alpha-1)*exp(-sigma(x1)*NH-sigmahe1(x1)*Nhe1-sigmahe2(x1)*Nhe2)*PHe1(x1,nh1,nhe1l,nhe2l)+pow(x2/x0,alpha-1)*exp(-sigma(x2)*NH-sigmahe1(x2)*Nhe1-sigmahe2(x2)*Nhe2)*PH1(x2,nh1,nhe1l,nhe2l)) ; 
     fsum = fsum + df*pow(u[i]/x0,alpha-1)*exp(-sigma(u[i])*NH-sigmahe1(u[i])*Nhe1-sigmahe2(u[i])*Nhe2)*PHe1(u[i],nh1,nhe1l,nhe2l);
}

/*
    for(int i=0;i<len;i++){
     fsum = fsum + df*pow(u[i]/x0,alpha-1)*exp(-sigma(u[i])*NH-sigmahe1(u[i])*Nhe1-sigmahe2(u[i])*Nhe2)*PHe1(u[i],nh1,nhe1l,nhe2l) ;   
     
}
*/
// if (fsum==0)
// {
//     printf("fsum zero %0.12e,%0.12e,%0.12e\n", pow(u[19]/x0,alpha-1),exp(-sigma(u[0])*NH-sigmahe1(u[0])*Nhe1-sigmahe2(u[0])*Nhe2),PHe1(u[19],nh1,nhe1l,nhe2l));
// }
  /*isothermal case*/
   //  fsum  = df*exp(-sigma(u[0])*NH)*PH1(u[0],nh1,nhe1,nhe2) ;  
   
 //   fsum  = df*exp(-sigma(u[0])*NH-sigmahe1(u[0])*Nhe1-sigmahe2(u[0])*Nhe2)*PHe1(u[0],nh1,nhe1,nhe2) ;
    return fsum;
}

double Gammahe2(double NH,double Nhe1, double Nhe2,double nh1,double nhe1l,double nhe2l){
    
    //double u[] ={13.6, 14.4059, 15.2595, 16.1636, 17.1214, 18.1359, 19.2105, 20.3488,21.5545, 22.8317, 24.1846, 25.6176, 27.1356, 28.7435, 30.4466, 32.2507, 34.1617, 36.1859, 38.33, 40.6012, 43.007, 45.5553, 48.2546, 51.1139, 54.1426,57.3507, 60.749, 64.3486, 68.1615, 72.2003, 76.4784, 81.0101, 85.8102, 90.8948, 96.2806, 101.986, 108.029, 114.43, 121.21, 128.392, 136};
    //double u[] = {54.5,57.7293, 61.15, 64.7734, 68.6114, 72.6769, 76.9833, 81.5448, 86.3767, 91.4948, 96.9162, 102.659, 108.742, 115.185, 122.01, 129.24, 136.898, 145.01, 153.602, 162.703, 172.344, 182.556, 193.373, 204.831, 216.968};
    // SPACING = (b/a)^(1/n) , n= 40, b = 136, a =13.6 
    double u[]={54.5,61.1500057595,68.6114349428,76.9832961819,86.3766789891,96.9162278471,108.741796166,122.010302052,136.897810517,153.601869754,172.344132479,193.373297132,216.968407952,243.442557722,273.147042327,306.476022229,343.871752742,385.830452489,432.908887925,485.731761283,545.0};

    double fsum=0 ;
    double x0 = 13.6;//u[0];
    
    int len = sizeof(u)/sizeof(u[0]);
    //double df = 1;//0.5*pow(u[len-1]/x0,1.0/(len-1));
   // fsum = integral(&u[0],len,NH,nh1,nhe1,nhe2);
    //for(int i=0;i<sizeof(u)/sizeof(u[0]);i++){
 
for(int i=0;i<len-1;i++){
    // double df = 0.5*(u[i+1]-u[i-1]);//0.5*x0*pow(u[len-1]/x0,i/(lent))*(u[len-1]/x0-1);
     double x1 = u[i] ;//(u[i+1]+u[i])/2.0;
     double x2 = u[i+1];
     double df = (u[i+1]-u[i]);
   //  fsum = fsum + 0.5*df*(pow(x1/x0,alpha-1)*exp(-sigma(x1)*NH-sigmahe1(x1)*Nhe1-sigmahe2(x1)*Nhe2)*PHe2(x1,nh1,nhe1l,nhe2l)+pow(x2/x0,alpha-1)*exp(-sigma(x2)*NH-sigmahe1(x2)*Nhe1-sigmahe2(x2)*Nhe2)*PH1(x2,nh1,nhe1l,nhe2l)) ; 
     fsum = fsum + df*pow(u[i]/x0,alpha-1)*exp(-sigma(u[i])*NH-sigmahe1(u[i])*Nhe1-sigmahe2(u[i])*Nhe2)*PHe2(u[i],nh1,nhe1l,nhe2l);
}

/*
    for(int i=0;i<len;i++){
     fsum = fsum + df*pow(u[i]/x0,alpha-1)*exp(-sigma(u[i])*NH-sigmahe1(u[i])*Nhe1-sigmahe2(u[i])*Nhe2)*PHe2(u[i],nh1,nhe1l,nhe2l) ;   
     
}
 */
  /*isothermal case*/
 //    fsum  = df*exp(-sigma(u[0])*NH-sigmahe1(u[0])*Nhe1-sigmahe2(u[0])*Nhe2)*PHe2(u[0],nh1,nhe1,nhe2) ;  
   

    return fsum;
}

double Heat(double NH,double Nhe1,double Nhe2,double nh1,double nhe1l,double nhe2l){
    
    //double u[] ={13.6, 14.4059, 15.2595, 16.1636, 17.1214, 18.1359, 19.2105, 20.3488,21.5545, 22.8317, 24.1846, 25.6176, 27.1356, 28.7435, 30.4466, 32.2507, 34.1617, 36.1859, 38.33, 40.6012, 43.007, 45.5553, 48.2546, 51.1139, 54.1426,57.3507, 60.749, 64.3486, 68.1615, 72.2003, 76.4784, 81.0101, 85.8102, 90.8948, 96.2806, 101.986, 108.029, 114.43, 121.21, 128.392, 136};
    // SPACING = (b/a)^(1/n) , n= 40, b = 136, a =13.6 
    //double u[]={13.6,15.2594509785,17.1213856004,19.2105106069,21.5545474175,24.1845999765,27.1355674836,30.4466074845,34.1616554685,38.3300078652,43.0069761783,48.2546209358,54.1425751953,60.7489685325,68.1614637733,76.4784202259,85.8101988493,96.2806266762,108.028639923,121.210127586,136.0};
    double u[]={13.6,15.3521953269,17.3301398056,19.5629184807,22.083363653,24.9285376673,28.1402779032,31.7658119796,35.8584522227,40.4783796062,45.6935286935,51.580586594,58.2261206226,65.7278512446,74.1960890239,83.755356705,94.5462202801,106.727355968,120.477883496,136.0};
    double fsum = 0 ;
    double x0 = u[0];
    
    int len = sizeof(u)/sizeof(u[0]);   
    double df = 1;//0.5*pow(u[len-1]/x0,1.0/(len-1));
   // for(int i=0;i<sizeof(u)/sizeof(u[0]);i++){
 //   fsum  = hintegral(&u[0],len,NH,nh1,nhe1,nhe2);

    for(int i=0;i<len-1;i++){	
     df = u[i+1]-u[i]; 
     fsum = fsum + df*(u[i]-x0)*pow(u[i]/x0,alpha-1)*exp(-sigma(u[i])*NH-sigmahe1(u[i])*Nhe1-sigmahe2(u[i])*Nhe2)*PH1(u[i],nh1,nhe1l,nhe2l) ;	
    // printf("%0.12e\n",PH1(u[i],nh1,nhe1l,nhe2l));
}

        
    return fsum;
}
double Heat2(double NH,double Nhe1, double Nhe2,double nh1,double nhe1l,double nhe2l){
    
    //double u[] ={13.6, 14.4059, 15.2595, 16.1636, 17.1214, 18.1359, 19.2105, 20.3488,21.5545, 22.8317, 24.1846, 25.6176, 27.1356, 28.7435, 30.4466, 32.2507, 34.1617, 36.1859, 38.33, 40.6012, 43.007, 45.5553, 48.2546, 51.1139, 54.1426,57.3507, 60.749, 64.3486, 68.1615, 72.2003, 76.4784, 81.0101, 85.8102, 90.8948, 96.2806, 101.986, 108.029, 114.43, 121.21, 128.392, 136};
    //double u[] = {54.5,57.7293, 61.15, 64.7734, 68.6114, 72.6769, 76.9833, 81.5448, 86.3767, 91.4948, 96.9162, 102.659, 108.742, 115.185, 122.01, 129.24, 136.898, 145.01, 153.602, 162.703, 172.344, 182.556, 193.373, 204.831, 216.968};
    // SPACING = (b/a)^(1/n) , n= 40, b = 136, a =13.6 
    //double u[]={54.5,61.1500057595,68.6114349428,76.9832961819,86.3766789891,96.9162278471,108.741796166,122.010302052,136.897810517,153.601869754,172.344132479,193.373297132,216.968407952,243.442557722,273.147042327,306.476022229,343.871752742,385.830452489,432.908887925,485.731761283,545.0};
    double u[]={24.6,27.6016539758,30.9695651301,34.7484235977,38.9883725345,43.745673487,49.0834529482,55.0725400088,61.7924062151,69.3322201091,77.7920304401,87.2840937515,97.9343639562,109.884163669,123.292059472,138.335965997,155.215506742,174.154662958,195.404745742,219.247730781,246.0};
    double fsum = 0 ;
    double x0 = u[0];
    
    int len = sizeof(u)/sizeof(u[0]);
    double df = 1;//0.5*pow(u[len-1]/x0,1.0/(len-1));

    for(int i=0;i<len-1;i++){
      df = u[i+1]-u[i];
     fsum = fsum + df*(u[i]-x0)*pow(u[i]/13.6,alpha-1)*exp(-sigma(u[i])*NH-sigmahe1(u[i])*Nhe1-sigmahe2(u[i])*Nhe2)*PHe1(u[i],nh1,nhe1l,nhe2l) ;   
     
}
 
  /*isothermal case*/
   //  fsum  = df*exp(-sigma(u[0])*NH)*PH1(u[0],nh1,nhe1,nhe2) ;  
   
 //   fsum  = df*exp(-sigma(u[0])*NH-sigmahe1(u[0])*Nhe1-sigmahe2(u[0])*Nhe2)*PHe1(u[0],nh1,nhe1,nhe2) ;
    return fsum;
}

double Heat3(double NH,double Nhe1, double Nhe2,double nh1,double nhe1l,double nhe2l){
    
    //double u[] ={13.6, 14.4059, 15.2595, 16.1636, 17.1214, 18.1359, 19.2105, 20.3488,21.5545, 22.8317, 24.1846, 25.6176, 27.1356, 28.7435, 30.4466, 32.2507, 34.1617, 36.1859, 38.33, 40.6012, 43.007, 45.5553, 48.2546, 51.1139, 54.1426,57.3507, 60.749, 64.3486, 68.1615, 72.2003, 76.4784, 81.0101, 85.8102, 90.8948, 96.2806, 101.986, 108.029, 114.43, 121.21, 128.392, 136};
    //double u[] = {54.5,57.7293, 61.15, 64.7734, 68.6114, 72.6769, 76.9833, 81.5448, 86.3767, 91.4948, 96.9162, 102.659, 108.742, 115.185, 122.01, 129.24, 136.898, 145.01, 153.602, 162.703, 172.344, 182.556, 193.373, 204.831, 216.968};
    // SPACING = (b/a)^(1/n) , n= 40, b = 136, a =13.6 
    double u[]={54.5,61.1500057595,68.6114349428,76.9832961819,86.3766789891,96.9162278471,108.741796166,122.010302052,136.897810517,153.601869754,172.344132479,193.373297132,216.968407952,243.442557722,273.147042327,306.476022229,343.871752742,385.830452489,432.908887925,485.731761283,545.0};

    double fsum=0 ;
    double x0 = u[0];
    
    int len = sizeof(u)/sizeof(u[0]);
    double df = 1;//0.5*pow(u[len-1]/x0,1.0/(len-1));

    for(int i=0;i<len-1;i++){
      df = u[i+1] - u[i];
     fsum = fsum + df*(u[i]-x0)*pow(u[i]/13.6,alpha-1)*exp(-sigma(u[i])*NH-sigmahe1(u[i])*Nhe1-sigmahe2(u[i])*Nhe2)*PHe2(u[i],nh1,nhe1l,nhe2l) ;   
     
}
 
  /*isothermal case*/
   //  fsum  = df*exp(-sigma(u[0])*NH)*PH1(u[0],nh1,nhe1,nhe2) ;  
   
 //   fsum  = df*exp(-sigma(u[0])*NH-sigmahe1(u[0])*Nhe1-sigmahe2(u[0])*Nhe2)*PHe1(u[0],nh1,nhe1,nhe2) ;
    return fsum;
}

double Hubble(double t){
    double  omega_m = 0.32;
    double  H0 = 67.32*1.05e-6;
   // double z =  0;//5.945518;
    return sqrt(omega_m)*H0*cosh(1.5*H0*sqrt(1-omega_m)*t)*pow(1+z,1.5);///sinh(1.5*H0*sqrt(1-omega_m)*t);
}

double xh2solve(double Ti,double xh2,double xhe1,double xhe2,double xhe3,double dti,double NH1,double Nhe1, double Nhe2,double R,double nhi, double *xrate){
    //double nhe1l = (Y/(4*(1-Y)))*nhi*(xhe1);
    double nhe1l = max((Y/(4*(1-Y)))*nhi*(x_min_fix),(Y/(4*(1-Y)))*nhi*(xhe1));
     double nhe2l = max((Y/(4*(1-Y)))*nhi*xhe2,(Y/(4*(1-Y)))*nhi*x_min_fix);
    double nhe3l = max((Y/(4*(1-Y)))*nhi*xhe3,(Y/(4*(1-Y)))*nhi*x_min_fix);
    //double nhe2l = (Y/(4*(1-Y)))*nhi*xhe2;
    //double nhe3l = (Y/(4*(1-Y)))*nhi*xhe3;

    //printf("%0.12e %0.12e \n",Nhe1,Nhe2 );
 //   double nhe2 =0;
 //   double nhe3 =0;
 //   double nhe1 =0;
 //   double xhe2,xhe3;
   
    double ne = xh2*nhi+(Y/(4*(1-Y)))*nhi*(xhe2+2*xhe3);	
//    xhe2 =0;
//    xhe3= 0;
    double fac = -alpha*Ng/(pow(3e21,3)*13.6);
    /*isothermal case */
   // double fac = Ng/pow(3e21,3);
    double vol = 1.0/(4*PI*R*R*cellsize);
    double nh1 = max(nhi*x_min_fix,nhi*(1-xh2));
    double C,D,dbdt,dalphadt,nu,Tnu,dCdx,dDdx,J;
    
    double Gammab; 
    if (background==1)
    {
      Gammab = 0.23e-12*3.15e13; // from haardt madau 2012
    }
    else
      Gammab = 0;//0.23e-12*3.15e13; // from haardt madau 2012
    
   // Gammab = 3.15e13*1.27e-11*Jb*pow(-alpha+3,-1);
    
    C = GE(Ti)*ne + fac*vol*Gamma(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l)/(nh1) + Gammab;
    D = ralpha(Ti)*ne;
    if (C !=C)
    {
        printf("H1 gone wrong\n");
        printf("gamma %0.12e %0.12e %0.12e\n", fac*vol, Gamma(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l),1/(nh1));
        printf("elec coll %0.12e %0.12e\n", GE(Ti),ne);
        exit(0);
    }
   // printf("%lf,%lf\n",C,Ti);
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
    *(xrate) = (double) max(0.0,(C-xh2*(C+D)));
   //return (xh2 + C*dti)/(1+D*dti);
    if ((1-J*dti)==0 || *(xrate)<0)
    {
        printf("x denom is zero\n");
        printf("%0.12e %0.12e %0.12e %0.12e %0.12e negative dervative\n",(C-xh2*(C+D)),xh2,C,D,xh2*(C+D) );
        exit(0);
    }
   return xh2 + dti*(*(xrate))/(1-J*dti);
}

double xhe1solve(double Ti,double xh2i,double xhe1,double xhe2,double xhe3,double dti,double NH1,double Nhe1,double Nhe2,double R,double nhi){
        double nhe1l = max((Y/(4*(1-Y)))*nhi*(x_min_fix),(Y/(4*(1-Y)))*nhi*(xhe1));
        double nhe2l = max((Y/(4*(1-Y)))*nhi*xhe2,(Y/(4*(1-Y)))*nhi*x_min_fix);//(Y/(4*(1-Y)))*nhi*xhe2;
        double nhe3l = max((Y/(4*(1-Y)))*nhi*xhe3,(Y/(4*(1-Y)))*nhi*x_min_fix);//(Y/(4*(1-Y)))*nhi*xhe3;
   // double xhe1 = 1- xhe2 -xhe3;

        double ne = xh2i*nhi+(Y/(4*(1-Y)))*nhi*(xhe2+2*xhe3); 

        double fac = -alpha*Ng/(pow(3e21,3)*13.6);
    /*isothermal case */
    //double fac = Ng/pow(3e21,3);
        double vol = 1.0/(4*PI*R*R*cellsize);
    
        double nh1 = max(nhi*(1-xh2i),nhi*x_min_fix);//nhi*(1-xh2i);
        double C,D,nu,Tnu;
        nu = pow((X*(1+xh2i)+Y*(1+xhe2+2*xhe3)/4),(-1));
        Tnu = Ti/nu;
        double Gammab;
        if (background==1)
         {
           Gammab = 0.154e-12*3.15e13; // from haardt madau 2012
         } 
         else
          Gammab = 0;//0.154e-12*3.15e13; // from haardt madau 2012


        C = xhe2*alphahe2(Ti)*ne;//+ xhe3*alphahe3(Ti)*ne;
   // D = betahe2(Ti)*ne + alphahe2(Ti)*ne + fac*vol*Gammahe2(NH1,Nhe1,Nhe2,nh1,nhe1,nhe2)/(nhe2) ;
        D = betahe1(Ti)*ne*nhe1l + fac*vol*Gammahe1(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l) + nhe1l*Gammab;
        if (D !=D)
    {
        printf("He1 gone wrong\n");
        printf("gamma %0.12e %0.12e %0.12e\n", fac*vol, Gammahe1(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l),1/(nhe1l));
        printf("elec coll %0.12e,%0.12e, %0.12e\n",xhe2, betahe1(Ti),ne);
        exit(0);
    }
  // printf("C %e D %e\n",C,D );
        if ((nhe1l+D*dti)==0)
    {
        printf("xhe1 denom is zero\n");
        exit(0);
    }
   return nhe1l*(xhe1 + C*dti)/(nhe1l+D*dti);
   //     return xhe1 + (C-D*xhe1)*dti;
}


double xhe2solve(double Ti,double xh2i,double xhe1,double xhe2,double xhe3,double dti,double NH1,double Nhe1,double Nhe2,double R,double nhi,double *hrate){
    //double nhe1l = (Y/(4*(1-Y)))*nhi*(xhe1);
    double nhe1l = max((Y/(4*(1-Y)))*nhi*(x_min_fix),(Y/(4*(1-Y)))*nhi*(xhe1));
    double nhe2l = max((Y/(4*(1-Y)))*nhi*xhe2,(Y/(4*(1-Y)))*nhi*x_min_fix);
    double nhe3l = max((Y/(4*(1-Y)))*nhi*xhe3,(Y/(4*(1-Y)))*nhi*x_min_fix);//(Y/(4*(1-Y)))*nhi*xhe3;
   // double xhe1 = 1- xhe2 -xhe3;

    double ne = xh2i*nhi+(Y/(4*(1-Y)))*nhi*(xhe2+2*xhe3); 

    double fac = -alpha*Ng/(pow(3e21,3)*13.6);
    /*isothermal case */
   //double fac = Ng/pow(3e21,3);
    double vol = 1.0/(4*PI*R*R*cellsize);
    
    double nh1 = max(nhi*(1-xh2i),nhi*x_min_fix);
    double C,D,nu,Tnu;
    nu = pow((X*(1+xh2i)+Y*(1+xhe2+2*xhe3)/4),(-1));
    Tnu = Ti/nu;
    double Gammab1, Gammab2;
    if (background==1)
    {
      Gammab1 =0.154e-12*3.15e13; 
      Gammab2 =0.442e-18*3.15e13; // from haardt madau 2012
    }
    else{
      Gammab1 =0;
      Gammab2 = 0;
    }

    
    C = xhe1*ne*betahe1(Ti) + fac*vol*Gammahe1(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l)/((Y/(4*(1-Y)))*nhi) + xhe3*alphahe3(Ti)*ne + xhe1*Gammab1;
    D = betahe2(Ti)*ne + alphahe2(Ti)*ne + fac*vol*Gammahe2(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l)/(nhe2l) + Gammab2;
   // D = alphahe2(Ti)*ne;
  // printf("C %e D %e\n",C,D );
    if (C !=C)
    {
        printf("He2 gone wrong\n");
        printf("gamma %0.12e %0.12e %0.12e\n", fac*vol, Gammahe1(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l),1/(nhe2l));
        printf("elec coll %0.12e,%0.12e, %0.12e\n",xhe1, betahe1(Ti),ne);
        exit(0);
    }
    *(hrate) = (C-D*xhe2);
    if ((1+D*dti)==0)
    {
        printf("xhe2 denom is zero\n");
        exit(0);
    }
    return (xhe2 + C*dti)/(1+D*dti);
   // return xhe2 + (C-D*xhe2)*dti;
}

double xhe3solve(double Ti,double xh2i,double xhe1,double xhe2,double xhe3,double dti,double NH1,double Nhe1,double Nhe2,double R,double nhi,double *hrate){
    //double nhe1l = (Y/(4*(1-Y)))*nhi*(xhe1);
    double nhe1l = max((Y/(4*(1-Y)))*nhi*(x_min_fix),(Y/(4*(1-Y)))*nhi*(xhe1));
    double nhe2l = max((Y/(4*(1-Y)))*nhi*xhe2,(Y/(4*(1-Y)))*nhi*x_min_fix);
    double nhe3l = max((Y/(4*(1-Y)))*nhi*xhe3,(Y/(4*(1-Y)))*nhi*x_min_fix);
    //double nhe2l = (Y/(4*(1-Y)))*nhi*xhe2;
    //double nhe3l = (Y/(4*(1-Y)))*nhi*xhe3;
    //double xhe1 = 1- xhe2 -xhe3;

    double ne = xh2i*nhi+(Y/(4*(1-Y)))*nhi*(xhe2+2*xhe3);    
    double fac = -alpha*Ng/(pow(3e21,3)*13.6);
    /*isothermal case */
   // double fac = Ng/pow(3e21,3);
    double vol = 1.0/(4*PI*R*R*cellsize);
    //double nh1 = nhi*(1-xh2i);
    double nh1 = max(nhi*(1-xh2i),nhi*x_min_fix);
    double C,D,nu,Tnu;
    nu = pow((X*(1+xh2i)+Y*(1+xhe2+2*xhe3)/4.0),(-1));
    Tnu = Ti/nu;
    double Gammab; 
    if (background==1)
    {
      Gammab = 0.442e-18*3.15e13; // from haardt madau 2012
    }
    else
      Gammab = 0;

    C = xhe2*ne*betahe2(Ti) + fac*vol*Gammahe2(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l)/((Y/(4*(1-Y)))*nhi)+ xhe2*Gammab;
    D = alphahe3(Ti)*ne;
  //  printf("C %e D %e\n",C,D );
    *(hrate) = (C-D*xhe3);
    if (C !=C)
    {
        printf("He3 gone wrong\n");
        printf("gamma %0.12e %0.12e %0.12e\n", fac*vol, Gammahe2(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l),xhe2/(nhe2l));
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

double Tsolve(double Ti,double NH1,double Nhe1,double Nhe2,double xh2i,double xhe1,double xhe2,double xhe3,double nhi,double dti,double R,double tsim, double *trate, double *xrate){
    
    
    double fac = -alpha*Ng/(pow(3e21,3)*13.6);
    //double fac = Ng/(pow(3e21,3));
    double fh = fac*1.6e-12; // 1.6e-12 is factor ev to erg
    double vol = 1.0/(4*PI*R*R*cellsize);
    
    //double K = (2.0/3)*(1.0/(nhi*kb)); // only H # cm^3erg^-1 K    
     double nh1 = max(nhi*(1-xh2i),nhi*x_min_fix);
    //double nh1 = nhi*(1.0-xh2i);
    //double nhe1l = (Y/(4*(1-Y)))*nhi*(xhe1);
    double nhe1l = max((Y/(4*(1-Y)))*nhi*(x_min_fix),(Y/(4*(1-Y)))*nhi*(xhe1));
    double nhe2l = max((Y/(4*(1-Y)))*nhi*xhe2,(Y/(4*(1-Y)))*nhi*x_min_fix);
    double nhe3l = max((Y/(4*(1-Y)))*nhi*xhe3,(Y/(4*(1-Y)))*nhi*x_min_fix);
   // double nhe2l = (Y/(4*(1-Y)))*nhi*xhe2;
   // double nhe3l = (Y/(4*(1-Y)))*nhi*xhe3;
    double ne = xh2i*nhi+(Y/(4*(1-Y)))*nhi*(xhe2+2*xhe3);
    double n = nhi + (Y/(4*(1-Y)))*nhi + ne;

    double K = (2.0/3)*(1.0/(n*kb)); // From Chen
    
    double T5 = Ti/100000.0;
    
    double f1 = 1.778e-29*3.15e13;	//#*6.24e11 # in erg cm^3 Myr^-1 K^-1 #4.05   #
    double f2 = 7.5e-19*3.15e13;	//#*6.24e11 # in erg cm^3 Myr^-1 #1.711e+11
    double f3 = 1.42e-27*3.15e13;	//#*6.24e11# in erg cm^3 K^-1 #489.62
    double gff = 1.5 ;	//# from cen #1.1+ 0.34*(np.exp(-5.5+np.log(T)))**2
    double f4 = 5.65e-36*3.15e13;	//#*6.24e11 #in erg cm^3 Myr^-1 K^-1 #1.289e-6
    double f5 = 1.27e-21*3.15e13;

    double f11 = 1.0;
    double f12 = 8*1.778e-29*3.15e13;
    double f21 = 3.15e13*5.54e-17;
   
    double f51 = 3.15e13*9.38e-22;
    double f52 = 3.15e13*4.95e-22;
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
    
    //double z = 0;//5.945518;
    
  //  double X = 1;
  //  double Y = 0;
  //  double xhe2 = 0; 
  //  double xhe3 = 0;
    double nu = pow((X*(1+xh2i)+Y*(1+xhe2+2*xhe3)/4),(-1));
    if (K*nu==0 || nu!=nu)
    {
        printf("Knu is zero\n");
        printf("nu %0.12e, %0.12e,%0.12e,%0.12e\n",nu,xh2i,xhe2,xhe3);
        exit(0);
    }


    
    double c1 = (f1)*Ti*pow(lamb,1.965)/pow((1.0+pow((lamb/0.541),0.502)),2.697) ;// recombination cooling coeffecient erg cm^3 Myr^-1
    double c2 = (f2)*pow((1+pow(T5,1.0/2)),-1)*exp(-118355.0/Ti);  //collisional excitation cooling coeffecient
    double c3 = gff*(f3)*pow(Ti,0.5); // Brehmstrallung/ free free emission
    double c4 = (f4)*(Ti-2.73*(1+z))*pow(1+z,4);  //inverse compton
    double c5 = f5*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-157809.1/Ti);// collisional Ionisation cooling

    double c11 = f11*kb*Ti*alphahe2(Ti); // recombination cooling He II erg cm^3 Myr^-1
    double c12 = f12*Ti*pow(lamb2,1.965)/pow((1.0+pow((lamb2/0.541),0.502)),2.697) ; // recombination cooling He III
        
    double c21 = f21*pow(Ti,-0.397)*pow((1+pow(T5,1.0/2)),-1)*exp(-473638.0/Ti);    //collisional excitation HeII erg cm^3 Myr^-1

    double c51 = f51*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-285335.4/Ti);
    double c52= f52*pow(Ti,0.5)*pow((1+pow(T5,0.5)),-1)*exp(-631515.0/Ti);
    //double C =  c1*pow((xh2i*nhi),2) +c2*xh2i*nhi*(nh1)+c3*pow((nhi*xh2i),2) +c4*(nhi*xh2i) + 2*Hubble(tsim)*Ti/(K*nu) ;//+ (Ti/(nhi*xh2i))*(*xrate);
    
    double C; 

  //  C = c1*ne*xh2i*nhi + c2*ne*nh1 + c3*(nhi*xh2i+nhe2+4*nhe3)*ne + c4*ne+ 2*Hubble(tsim)*Ti/(K*nu);

    C = c1*ne*xh2i*nhi + c11*ne*nhe2l + c12*ne*nhe3l +c2*ne*nh1 + c21*ne*nhe2l + c3*(nhi*xh2i+nhe2l+4*nhe3l)*ne + c4*ne + c5*ne*nh1+ c51*ne*nhe1l+ c52*ne*nhe2l +2*Hubble(tsim)*Ti/(K*nu)+ *(xrate)*(Ti/n)/(K*nu);
    if (C != C)
    {
        printf("Cooling is gone bonkers\n");
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
       bheat1 =0.989e-12*3.15e13*1.6e-12; // background photoheating from haardt madau 2012 erg Myr^-1
       bheat2 =0.112e-11*3.15e13*1.6e-12; 
       bheat3 =0.269e-16*3.15e13*1.6e-12;    
    }
    else{
       bheat1 =0;
       bheat2 =0; 
       bheat3 =0; 
    }
    //printf("background %0.12e\n", bheat3);
    double L = (fh*vol*Heat(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l)+ nh1*bheat1 +fh*vol*Heat2(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l)+nhe1l*bheat2+fh*vol*Heat3(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l)+nhe3l*bheat3 - C);
    
   //  printf("%e, %e\n",Heat(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l), C);
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
    double dLdt = -(A1+ A2 +A3 +A4 +A11+ A12+ A21 +2*Hubble(tsim)/(K*nu));//+(1/(nhi*xh2i))*(*xrate));
    //double dLdt =0;
   // printf("%e\n", K)
	//printf("%e %e\n",fh*vol*Heat(NH1,nh1,nhe1,nhe2)-C, dti);
    *(trate) = L*K*nu;
   return Ti+ dti*L*K*nu;
   //return Ti + nu*dti*L*K/(1-nu*dLdt*K*dti) ;
   //return (L*K*nu*dti + Ti)/(1+2*Hubble(tsim)*dti) ; // anninos implicit
}


long int xTiter(double dti,double NH1,double Nhe1, double Nhe2,double R,double nhi,double tsim,double ti,double *Tt, double *xh2t,double *xhe1t,double *xhe2t,double *xhe3t,double *dtrec){
  //  double *Tt =(double*)malloc(countl*sizeof(double));
  //  Tt[0] = T[k][j];
  //  double *xh2t =(double*)malloc(countl*sizeof(double)); 	
   // xh2t[0] = nh2[k][j]/nhi;
   // double *NHt =(double*)malloc(countl*sizeof(double)); 	
   // NHt[0] = NH[k][j];
   // printf("%d %lf\n",k,Tt[0]);
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
    for (int i=0;i<countl;i++){
        

        //Tt[i+1] = (double) max(T2_MIN_FIX,Tsolve(Tt[i],NH1,Nhe1,Nhe2,xh2t[i],xhe1t[i],xhe2t[i],xhe3t[i],nhi,dti,R,tsim,&Trate));
        Tt[i+1] = (double) max(T2_MIN_FIX,Tsolve(Tt[i],NH1,Nhe1,Nhe2,xh2t[i],xhe1t[i],xhe2t[i],xhe3t[i],nhi,dti,R,tsim,&Trate,&nerate));
      //  printf("T %e\n",Tt[i+1] );
        dUU = (double) fabs(max(T2_MIN_FIX,Tt[i]+Trate*dti)-Tt[i]);
     //   printf("heat %e\n",Heat(NH1,nhi,0,0)*(-alpha*Ng/pow(3e21,3))* 1.6e-12*1.0/(4*PI*R*R*cellsize*kb));
 //   printf("%e\n",fabs(Tt[i+1]-Tt[i])/(Tt[i]+T_MIN));
 //   printf("%e\n",Tt[i+1]);

//if(dti>1e-7){
      //  printf("i %d\n", i);
      //  printf("change %e\n",(double) max(dUU,fabs(Tt[i+1]-Tt[i]))/(Tt[i]+T_MIN) );

//printf("%0.12e %0.12e\n",dti,dt_min );
	if((double) max(dUU,fabs(Tt[i+1]-Tt[i]))/(Tt[i]+T_MIN)>lim  && dti>dt_min){
	//printf("hii");
  //  printf("change %e\n",(double) max(dUU,fabs(Tt[i+1]-Tt[i]))/(Tt[i]+T_MIN) );
	return 100000000000;
}
	else{
	//	printf("hi");
	//frac_change = (double) max(frac_change,fabs(Tt[i+1]-Tt[i])/Tt[i]);
        frac_change = (double) max(frac_change,dUU);

}

//}

xh2t[i+1] = (double) max(x_min_fix,xh2solve(Tt[i+1],xh2t[i],xhe1t[i],xhe2t[i],xhe3t[i],dti,NH1,Nhe1,Nhe2,R,nhi,&Xrate));
//printf("xh2 %e\n",xh2t[i+1] );
/*complete this exactly*/

     if((double) max(fabs(xh2t[i+1]-xh2t[i]),Xrate*dti)/(xh2t[i]+X_FM) >lim && dti>dt_min){
	return 100000000001;
}
else{
	frac_change = (double) max(frac_change,fabs(xh2t[i+1]-xh2t[i])/(xh2t[i])); //check
}

//xhe1t[i+1] = xhe1solve(Tt[i+1],xh2t[i],xhe2t[i],xhe3t[i],dti,NH1,Nhe1,Nhe2,R,nhi);
xhe1t[i+1] = (double) max(x_min_fix,xhe1solve(Tt[i+1],xh2t[i+1],xhe1t[i],xhe2t[i],xhe3t[i],dti,NH1,Nhe1,Nhe2,R,nhi));

xhe2t[i+1] = (double) max(x_min_fix,xhe2solve(Tt[i+1],xh2t[i+1],xhe1t[i],xhe2t[i],xhe3t[i],dti,NH1,Nhe1,Nhe2,R,nhi,&he2rate));

xhe3t[i+1] = (double) max(x_min_fix,xhe3solve(Tt[i+1],xh2t[i+1],xhe1t[i],xhe2t[i+1],xhe3t[i],dti,NH1,Nhe1,Nhe2,R,nhi,&he3rate));

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
  //  printf("He fracs don't add up to 1  but %0.12e\n",xhe1t[i+1]+xhe2t[i+1]+xhe3t[i+1]);
}


ne_init = xh2t[i]*nhi+(Y/(4*(1-Y)))*nhi*(xhe2t[i]+2*xhe3t[i]);
    ne_up = xh2t[i+1]*nhi+(Y/(4*(1-Y)))*nhi*(xhe2t[i+1]+2*xhe3t[i+1]);
    nerate =  Xrate*nhi+(Y/(4*(1-Y)))*nhi*(he2rate+2*he3rate);

    if(fabs(ne_up - ne_init)/(ne_init+X_FM)>lim && dti >dt_min){
    //printf("change %0.12e %0.12e\n", fabs(ne_up - ne_init)/(ne_init+X_FM), lim);
    return 100000000002;
}
else {
    frac_change = (double) max(frac_change,fabs(ne_up - ne_init)/(ne_init));
}

/*
	ne_init = xh2t[i]*nhi;
	if(fabs(xh2t[i+1]*nhi-ne_init)/ne_init>0.1 & dti>dt_min){
	return 100000000003;
}

else {
	frac_change = (double) max(frac_change,fabs(xh2t[i+1]*nhi-ne_init)/(ne_init+X_FM));
}
*/
	double tcal = ti + dti;
    ti = tcal ;

//if (dti>1e-7){
  

    if(frac_change <0.5 && dti>dt_min){  // has value 0.05 in paper but 0.5 here
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
	
   //     printf("tsim %lf\n", tcal);
        if (tcal >= tsim){
            storei = i+1;
            printf("tsim %lf\n", tcal);
            break ;
        }
 dti = (double) min(*(dtrec),tsim-tcal);
 // printf("i %d dti %e\n", i, dti);

   }

	printf("T %lf %d\n",Tt[storei],storei);
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
	long int out;
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
    xhe2t[0] = nhe2[k][j]/((Y/(4.0*(1-Y)))*nh[k]);
    xhe3t[0] = nhe3[k][j]/((Y/(4.0*(1-Y)))*nh[k]);
    xhe1t[0] = 1.0 -xhe2t[0] -xhe3t[0];
	out = xTiter(dt[k][j],NH[k][j],NHe1[k][j],NHe2[k][j],(1+k)*cellsize,nh[k],tsim,ti,Tt,xh2t,xhe1t,xhe2t,xhe3t,&dt_rec);  // equivalent of cool_step routine
	//printf("%d\n",out);
	

    while((out==100000000000 || out==100000000001 || out ==100000000002 || out==100000000003)){
	dt[k][j] =(double) max(dt_min,dt[k][j]/2);
	printf("%ld %e\n",out-100000000000, dt[k][j]);
    /*if (dt[k][j]<1e-8)
    {
        out = xTiter(dt[k][j],NH[k][j],(1+k)*cellsize,nh[k],tsim,ti,Tt,xh2t,&dt_rec);
        break;
    }*/
	out = xTiter(dt[k][j],NH[k][j],NHe1[k][j],NHe2[k][j],(1+k)*cellsize,nh[k],tsim,ti,Tt,xh2t,xhe1t,xhe2t,xhe3t,&dt_rec);

	//printf("%d\n",out);
} 
	printf("%ld \n",out);
    //int out = outi;
    //printf("%d\n",out);
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

	//dt[k][j+1] = Dt;//(double) min((dt_rec),tsim-ti); // check this criterion
    

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




void
absorb_along_lines_of_sight (double *Tau, double *Temp,double *NHI, double *Vpec)
{
  int ncells = boxsize/cellsize;  
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
    	
    	//printf("%lf\n",(double) max(3e-3,4.4e-3));
    double  omega_m = 0.32;
    double  H0 = 67.32*1.05e-6;
  //  double z =  0.0;    
    double a = 1/(1+z);
    double ti = (1.0/H0)*(2.0/3)*pow(1-omega_m,-0.5)*asinh(pow((1-omega_m)/omega_m,0.5)*pow(a,1.5));
	printf("ti %lf\n",ti );
  //  double ti = 935.7;
	double tf = ti + atof(argv[1]);
//	printf("%lf\n",tf);
	//double epsilon = 0.05;	
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


    double *Tau = (double*)malloc(ncells*sizeof(double));
    double *Vpec = (double*)malloc(ncells*sizeof(double));
   // double *rho = (double*)malloc(ncells*sizeof(double));
    double **rho = (double **)malloc(ncells * sizeof(int *));
        for (int i=0; i<r; i++)
            rho[i] = (double *)malloc(5000 * sizeof(int));

FILE *myFile,*myFile2, *myFile3,*myFile4,*myFile5;
    myFile = fopen("LOS9_nh.txt", "r");
    myFile2 = fopen("LOS9_nh1.txt", "r");
    myFile3 = fopen("LOS9_T.txt", "r");
    myFile4 = fopen("LOS9_Vpec.txt", "r");
    myFile5 = fopen("LOS_rho.txt", "r");
    double nharr[2048];
    double nh1arr[2048];
    double Tarr[2048];
    double vpec[2048];
    //double rho[5000][5000];
    
    for (int i = 0; i < 2048; i++)
    { 
        fscanf(myFile, "%lf", &nharr[i]);
        fscanf(myFile2, "%lf", &nh1arr[i]);
        fscanf(myFile3, "%lf", &Tarr[i]);
        fscanf(myFile4, "%lf", &vpec[i]);
        //fscanf(myFile5, "%lf", &rho[i]);
    }
    for (int i = 0; i < 2048; ++i)
    {
      for (int j = 0; j < 5000; ++j)
      {
        fscanf(myFile5,"%lf",&rho[i][j]);
      }
      
    }
/*
    for (int i = 0; i < 2048; i++)
    {
        printf("Number is: %0.12e\n", numberArray[i]);
    }          
*/

/* Initial values*/
	for (int i = 0; i < ncells; ++i){
	//printf("%d\n",i);
	nh[i] = nharr[i];//1e-3;
  Vpec[i] = vpec[i];
	nh2[i][0] = nharr[i]- 1e-4*nharr[i];//nh[i]- nh1arr[i];//x_min_fix*nh[i];
    nhe2[i][0] = x_min_fix*(Y/(4*(1-Y)))*nh[i];
    nhe3[i][0] = x_min_fix*(Y/(4*(1-Y)))*nh[i];
    nhe1[i][0] = (Y/(4*(1-Y)))*nh[i] - nhe2[i][0] - nhe3[i][0];
	T[i][0] = Tarr[i];
	NH[i][0] = cumsum(nh,nh2,i,0)*cellsize*3e7; // in nm^-2
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
double RIF,RIFold;
double NHION;
int cks;

for (int ck = 0; ck < ncells; ++ck)
    {
        if (nh2[ck][0]/nh[ck]<=0.5)
        {
            RIF = (1+ck)*cellsize;
            cks =ck;
            break; 
        }
    }
    printf("RIF %0.12e\n",RIF );
/*
NHION = NH[cks+1][0];//nm^-2//*(pow(3e21,2))/1e-14 ; // in kpc^-2 
   printf("%0.12e %d\n",NHION,cks );
   Dt = 4*PI*pow(RIF*3e28,2)*NHION*pow(Ng,-1); //in Myr
   printf("Dt %0.12e\n",Dt );
for (int i = 0; i < ncells; ++i)
{
  dt[i][0] =Dt;
}
*/

for (int j=0;j<count;j++){
	printf("Dt %0.12e\n",Dt );
  tsim = ti + Dt;
	//printf("%d\n",j);	
 //   printf("dt new out %e \n",dt[0][j]);
    /*evolve all cells for a  global timestep*/
    evolve_eqns(nh,nh2,NH,nhe1,nhe2,nhe3,NHe1,NHe2,T,ncells,ti,tsim,dt,j,Dt);	// my equivalent of subroutine rt_solve_cooling
//	printf("dt new out %e \n",dt[0][j+1]);
    for (int ck = 0; ck < ncells; ++ck)
    {
        if (nh2[ck][j]/nh[ck]<=0.5)
        {
            RIFold = RIF;
            RIF = (1+ck)*cellsize;
            cks =ck;
            break; 
        }
    }
    printf("RIF %0.12e\n",RIF );
    ti = tsim;
	printf("%lf\n",tsim);
	if(fabs(tf-tsim)<=Dt){
	printf("%lf,%lf\n",tf,tsim);
  //printf("%0.12e %0.12e %0.12e\n",tsim-ti, tsim, ti);
	storej = j+1;
	break;
} 
  
   //if (cks==ncells-1)

   if(fabs(RIFold - RIF)<epsilon)
   {
        Dt = 1e0;//3261.6*(epsilon/0.1)*(cellsize/10.0)*1e-6;//nm^-2//*(pow(3e21,2))/1e-14 ; // in kpc^-2 
   }
 
   /*
   else{
   NHION = NH[cks+1][j];//nm^-2//*(pow(3e21,2))/1e-14 ; // in kpc^-2 
   printf("%0.12e %d\n",NHION,cks );
   Dt = 4*PI*pow(RIF*3e28,2)*NHION*pow(Ng,-1); //in Myr
   printf("Dt %0.12e\n",Dt );    
   }
   */
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
dense = pow(3e21,3)/(3e7*1.989e43*pow(Hubble_h,2));//*Hubble_h*Hubble_h;

for (int i = 0; i < ncells; ++i)
{
  Tf[i] = T[i][storej];
  NH1f[i] = rho[i][9]*(1-nh2[i][storej]/nh[i])* (UnitMass_in_g / PROTONMASS / Hubble_h) ;//NH[i][storej]; 
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


 FILE *fptr = fopen("./speceiler_9_10myr.txt","w");
// FILE *fdptr = fopen("./c2raytest/dgc2xhe2_10myr_b1.txt","w");
//FILE *fd1ptr = fopen("./c2raytest/dgc2xhe3_10myr_b1.txt","w");
//FILE *fd2ptr = fopen("./c2raytest/dgc2T_10myr_b1.txt","w");
//FILE *fd3ptr = fopen("./c2raytest/dgc2xhe1_10myr_b1.txt","w");
/*
FILE *fptr = fopen("intblum_alphacrct_xh1_001myr.txt","w");
FILE *fdptr = fopen("intblum_alphacrct__T_001myr.txt","w");    
*/
     	for (int l = 0; l < ncells; ++l)
    {
    fprintf(fptr,"%lf %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e\n", (double)l*cellsize+cellsize,(nh[l]-nh2[l][storej])/nh[l],T[l][storej],nhe1[l][storej]/((Y/(4*(1-Y))*nh[l])),nhe2[l][storej]/((Y/(4*(1-Y))*nh[l])),nhe3[l][storej]/((Y/(4*(1-Y))*nh[l])),Tau[l]);
	// fprintf(fd2ptr,"%0.12e \n",T[l][storej]);
  //  fprintf(fd3ptr,"%0.12e \n", nhe1[l][storej]/((Y/(4*(1-Y))*nh[l])));
  //  fprintf(fdptr,"%0.12e \n",nhe2[l][storej]/((Y/(4*(1-Y))*nh[l])));
  //  fprintf(fd1ptr,"%0.12e \n",nhe3[l][storej]/((Y/(4*(1-Y))*nh[l])));
    } 
	

}


