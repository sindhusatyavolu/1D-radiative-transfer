#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h> 
#include<time.h>

#ifndef PI
#define PI 3.1415
#endif
#define cellsize 50
#define boxsize 5000
#define count 100000
#define countl 10000000000
#define Ng 3.15e13*1e57
#define alpha -1.5
#define kb 1.38e-16
#define max(x,y) (((x) >= (y)) ? (x) : (y))
#define min(x,y) (((x) <= (y)) ? (x) : (y))
#define T_MIN 1e-1
#define T2_MIN_FIX 1e-2
#define X_MIN  1e-20
#define X_FM 1e-6
#define X 0.76
#define Y 0.24
#define lim 1e-1
#define epsilon 5e-2
#define dtmin 1e-15
#define x_min_fix 1e-20
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
    //printf("%0.12e\n",D );
    if (D ==0 || D!=D)
    {
        printf(" D %0.12e\n",D);
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
    double u[] = {13.6,14.2417926431,14.9138718889,15.6176669814,16.3546746116,17.1264621003,17.9346707311,18.7810192407,19.6673074745,20.5954202132,21.5673311819,22.5851072468,23.6509128112,24.7670144174,25.9357855677,27.1597117714,28.4413958305,29.7835633748,31.1890686586,32.6609006299,34.2021892874,35.8162123361,37.5064021582,39.2763531122,41.1298291765,43.0707719543,45.1033090552,47.2317628736,49.4606597805,51.7947397488,54.2389664342,56.7985377303,59.4788968225,62.2857437637,65.2250475959,68.3030590438,71.5263238082,74.9016964852,78.436355144,82.1378165913,86.0139523566,90.073005432,94.3236078017,98.7747987986,103.436044327,108.317256995,113.428817191,118.78159516,124.386974124,130.256874484,136.403779172,142.840760198,149.581506451,156.640352802,164.032310599,171.773099583,179.879181317,188.367794201,197.256990121,206.565672846,216.313638225,226.521616284,237.211315316,248.405468039,260.127879943,272.403479914,285.25837325,298.719897171,312.816678962,327.578696844,343.037343731,359.22549399,376.177573347,393.929632102,412.519421792,431.986475471,452.372191785,473.719923008,496.075067234,519.485164923,544.0};
    // SPACING = (b/a)^(1/n) , n= 80, b = 13.6*40, a =13.6 
    double fsum=0 ;
    double x0 = u[0];
    
    int len = sizeof(u)/sizeof(u[0]);
    double df = 1;//0.5*pow(u[len-1]/x0,1.0/(len-1));
   // fsum = integral(&u[0],len,NH,nh1,nhe1,nhe2);
    //for(int i=0;i<sizeof(u)/sizeof(u[0]);i++){
 

    for(int i=0;i<len;i++){
     fsum = fsum + df*pow(u[i]/x0,alpha-1)*exp(-sigma(u[i])*NH-sigmahe1(u[i])*Nhe1-sigmahe2(u[i])*Nhe2)*PH1(u[i],nh1,nhe1l,nhe2l) ;	
     
}
 
  /*isothermal case*/
   //  fsum  = df*exp(-sigma(u[0])*NH)*PH1(u[0],nh1,nhe1,nhe2) ;  
  //  fsum  = df*exp(-sigma(u[0])*NH-sigmahe1(u[0])*Nhe1-sigmahe2(u[0])*Nhe2)*PH1(u[0],nh1,nhe1,nhe2) ;
    return fsum;
}

double Gammahe1(double NH,double Nhe1, double Nhe2,double nh1,double nhe1l,double nhe2l){
    
    //double u[] ={13.6, 14.4059, 15.2595, 16.1636, 17.1214, 18.1359, 19.2105, 20.3488,21.5545, 22.8317, 24.1846, 25.6176, 27.1356, 28.7435, 30.4466, 32.2507, 34.1617, 36.1859, 38.33, 40.6012, 43.007, 45.5553, 48.2546, 51.1139, 54.1426,57.3507, 60.749, 64.3486, 68.1615, 72.2003, 76.4784, 81.0101, 85.8102, 90.8948, 96.2806, 101.986, 108.029, 114.43, 121.21, 128.392, 136};
   // double u[] = {54.5,57.7293, 61.15, 64.7734, 68.6114, 72.6769, 76.9833, 81.5448, 86.3767, 91.4948, 96.9162, 102.659, 108.742, 115.185, 122.01, 129.24, 136.898, 145.01, 153.602, 162.703, 172.344, 182.556, 193.373, 204.831, 216.968};
   // double u[]= {54.5,57.0718896359,59.7651483783,62.5855037123,65.5389534069,68.6317782696,71.8705555032,75.2621726926,78.8138424528,82.533117766,86.427908045,90.5064959524,94.7775550153,99.2501680696,103.933846576,108.838550849,113.974711232,119.353250289,124.985606022,130.883756201,137.060243835,143.528203847,150.301391002,157.394209163,164.821741921,172.59978467,180.744878199,189.274343869,198.206320444,207.55980267,217.354681667,227.611787228,238.35293212,249.600958465,261.379786322,273.714464551,286.631224084,300.157533709,314.322158482,329.155220899,344.688264958,360.954323239,377.987987147,395.82548048,414.504736459,434.065478399,454.549304183,475.999774723,498.462506601,521.985269071,546.618085651,572.413340501,599.42588982,627.713178509,657.335362328,688.355435828,720.839366309,754.856234114,790.47837953,827.781556626,866.845094356,907.752065258,950.589462113,995.448382951,1042.42422477,1091.61688642,1143.13098104,1197.07605852,1253.56683849,1312.72345426,1374.67170833,1439.54333988,1507.47630496,1578.61506982,1653.11091821,1731.12227303,1812.81503326,1898.36292676,1987.94787972,2081.76040355,2180.0};
    double u[]={24.6,25.7608896338,26.9765623873,28.2496035105,29.5827202534,30.9787476226,32.4406544106,33.971549509,35.57468852,37.2534806797,39.0114961084,40.8524734024,42.7803275849,44.7991584314,46.9132591886,49.1271257042,51.4454659875,53.8732102222,56.4155212501,59.0778055512,61.8657247404,64.7852076079,67.8424627273,71.0439916588,74.3966027752,77.9074257408,81.5839266734,85.4339240214,89.4656051911,93.6875439574,98.1087186972,102.738531483,107.586828076,112.663918867,117.980600798,123.548180329,129.378497477,135.483950995,141.877524746,148.572815305,155.58406088,162.92617159,170.614761171,178.66618018,187.097550769,195.926803094,205.172713448,214.854944187,224.994085548,235.611699434,246.730365267,258.373728006,270.566548433,283.334755804,296.705502996,310.707224245,325.369695618,340.724098334,356.803085072,373.640849413,391.273198553,409.737629456,429.073408587,449.321655424,470.525429897,492.729823963,515.982057496,540.331578707,565.830169299,592.532054585,620.494018808,649.775525893,680.438845906,712.549187478,746.174836476,781.38730122,818.261464553,856.875743088,897.312253967,939.656989493,984.0};
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
   // double u[] = {54.5,57.7293, 61.15, 64.7734, 68.6114, 72.6769, 76.9833, 81.5448, 86.3767, 91.4948, 96.9162, 102.659, 108.742, 115.185, 122.01, 129.24, 136.898, 145.01, 153.602, 162.703, 172.344, 182.556, 193.373, 204.831, 216.968};
    // SPACING = (b/a)^(1/n) , n= 40, b = 136, a =13.6 
    double u[] = {54.5,57.0718896359,59.7651483783,62.5855037123,65.5389534069,68.6317782696,71.8705555032,75.2621726926,78.8138424528,82.533117766,86.427908045,90.5064959524,94.7775550153,99.2501680696,103.933846576,108.838550849,113.974711232,119.353250289,124.985606022,130.883756201,137.060243835,143.528203847,150.301391002,157.394209163,164.821741921,172.59978467,180.744878199,189.274343869,198.206320444,207.55980267,217.354681667,227.611787228,238.35293212,249.600958465,261.379786322,273.714464551,286.631224084,300.157533709,314.322158482,329.155220899,344.688264958,360.954323239,377.987987147,395.82548048,414.504736459,434.065478399,454.549304183,475.999774723,498.462506601,521.985269071,546.618085651,572.413340501,599.42588982,627.713178509,657.335362328,688.355435828,720.839366309,754.856234114,790.47837953,827.781556626,866.845094356,907.752065258,950.589462113,995.448382951,1042.42422477,1091.61688642,1143.13098104,1197.07605852,1253.56683849,1312.72345426,1374.67170833,1439.54333988,1507.47630496,1578.61506982,1653.11091821,1731.12227303,1812.81503326,1898.36292676,1987.94787972,2081.76040355,2180.0};
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
    
   // double u[] ={13.6, 14.4059, 15.2595, 16.1636, 17.1214, 18.1359, 19.2105, 20.3488,21.5545, 22.8317, 24.1846, 25.6176, 27.1356, 28.7435, 30.4466, 32.2507, 34.1617, 36.1859, 38.33, 40.6012, 43.007, 45.5553, 48.2546, 51.1139, 54.1426,57.3507, 60.749, 64.3486, 68.1615, 72.2003, 76.4784, 81.0101, 85.8102, 90.8948, 96.2806, 101.986, 108.029, 114.43, 121.21, 128.392, 136};
    // SPACING = (b/a)^(1/n) , n= 40, b = 136, a =13.6 
    double fsum=0 ;
    double u[] = {13.6,14.2417926431,14.9138718889,15.6176669814,16.3546746116,17.1264621003,17.9346707311,18.7810192407,19.6673074745,20.5954202132,21.5673311819,22.5851072468,23.6509128112,24.7670144174,25.9357855677,27.1597117714,28.4413958305,29.7835633748,31.1890686586,32.6609006299,34.2021892874,35.8162123361,37.5064021582,39.2763531122,41.1298291765,43.0707719543,45.1033090552,47.2317628736,49.4606597805,51.7947397488,54.2389664342,56.7985377303,59.4788968225,62.2857437637,65.2250475959,68.3030590438,71.5263238082,74.9016964852,78.436355144,82.1378165913,86.0139523566,90.073005432,94.3236078017,98.7747987986,103.436044327,108.317256995,113.428817191,118.78159516,124.386974124,130.256874484,136.403779172,142.840760198,149.581506451,156.640352802,164.032310599,171.773099583,179.879181317,188.367794201,197.256990121,206.565672846,216.313638225,226.521616284,237.211315316,248.405468039,260.127879943,272.403479914,285.25837325,298.719897171,312.816678962,327.578696844,343.037343731,359.22549399,376.177573347,393.929632102,412.519421792,431.986475471,452.372191785,473.719923008,496.075067234,519.485164923,544.0};
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
    //double u[] = {54.5,57.7293, 61.15, 64.7734, 68.6114, 72.6769, 76.9833, 81.5448, 86.3767, 91.4948, 96.9162, 102.659, 108.742, 115.185, 122.01, 129.24, 136.898, 145.01, 153.602, 162.703, 172.344, 182.556, 193.373, 204.831, 216.968};
    // SPACING = (b/a)^(1/n) , n= 40, b = 136, a =13.6 
   // double u[]={54.5,57.0718896359,59.7651483783,62.5855037123,65.5389534069,68.6317782696,71.8705555032,75.2621726926,78.8138424528,82.533117766,86.427908045,90.5064959524,94.7775550153,99.2501680696,103.933846576,108.838550849,113.974711232,119.353250289,124.985606022,130.883756201,137.060243835,143.528203847,150.301391002,157.394209163,164.821741921,172.59978467,180.744878199,189.274343869,198.206320444,207.55980267,217.354681667,227.611787228,238.35293212,249.600958465,261.379786322,273.714464551,286.631224084,300.157533709,314.322158482,329.155220899,344.688264958,360.954323239,377.987987147,395.82548048,414.504736459,434.065478399,454.549304183,475.999774723,498.462506601,521.985269071,546.618085651,572.413340501,599.42588982,627.713178509,657.335362328,688.355435828,720.839366309,754.856234114,790.47837953,827.781556626,866.845094356,907.752065258,950.589462113,995.448382951,1042.42422477,1091.61688642,1143.13098104,1197.07605852,1253.56683849,1312.72345426,1374.67170833,1439.54333988,1507.47630496,1578.61506982,1653.11091821,1731.12227303,1812.81503326,1898.36292676,1987.94787972,2081.76040355,2180.0};
    double u[]={24.6,25.7608896338,26.9765623873,28.2496035105,29.5827202534,30.9787476226,32.4406544106,33.971549509,35.57468852,37.2534806797,39.0114961084,40.8524734024,42.7803275849,44.7991584314,46.9132591886,49.1271257042,51.4454659875,53.8732102222,56.4155212501,59.0778055512,61.8657247404,64.7852076079,67.8424627273,71.0439916588,74.3966027752,77.9074257408,81.5839266734,85.4339240214,89.4656051911,93.6875439574,98.1087186972,102.738531483,107.586828076,112.663918867,117.980600798,123.548180329,129.378497477,135.483950995,141.877524746,148.572815305,155.58406088,162.92617159,170.614761171,178.66618018,187.097550769,195.926803094,205.172713448,214.854944187,224.994085548,235.611699434,246.730365267,258.373728006,270.566548433,283.334755804,296.705502996,310.707224245,325.369695618,340.724098334,356.803085072,373.640849413,391.273198553,409.737629456,429.073408587,449.321655424,470.525429897,492.729823963,515.982057496,540.331578707,565.830169299,592.532054585,620.494018808,649.775525893,680.438845906,712.549187478,746.174836476,781.38730122,818.261464553,856.875743088,897.312253967,939.656989493,984.0};
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
    //double u[] = {54.5,57.7293, 61.15, 64.7734, 68.6114, 72.6769, 76.9833, 81.5448, 86.3767, 91.4948, 96.9162, 102.659, 108.742, 115.185, 122.01, 129.24, 136.898, 145.01, 153.602, 162.703, 172.344, 182.556, 193.373, 204.831, 216.968};
    // SPACING = (b/a)^(1/n) , n= 40, b = 136, a =13.6 
    double u[]={54.5,57.0718896359,59.7651483783,62.5855037123,65.5389534069,68.6317782696,71.8705555032,75.2621726926,78.8138424528,82.533117766,86.427908045,90.5064959524,94.7775550153,99.2501680696,103.933846576,108.838550849,113.974711232,119.353250289,124.985606022,130.883756201,137.060243835,143.528203847,150.301391002,157.394209163,164.821741921,172.59978467,180.744878199,189.274343869,198.206320444,207.55980267,217.354681667,227.611787228,238.35293212,249.600958465,261.379786322,273.714464551,286.631224084,300.157533709,314.322158482,329.155220899,344.688264958,360.954323239,377.987987147,395.82548048,414.504736459,434.065478399,454.549304183,475.999774723,498.462506601,521.985269071,546.618085651,572.413340501,599.42588982,627.713178509,657.335362328,688.355435828,720.839366309,754.856234114,790.47837953,827.781556626,866.845094356,907.752065258,950.589462113,995.448382951,1042.42422477,1091.61688642,1143.13098104,1197.07605852,1253.56683849,1312.72345426,1374.67170833,1439.54333988,1507.47630496,1578.61506982,1653.11091821,1731.12227303,1812.81503326,1898.36292676,1987.94787972,2081.76040355,2180.0};
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
    return sqrt(omega_m)*H0*cosh(1.5*H0*sqrt(1-omega_m)*t)*pow(1+z,1.5);///sinh(1.5*H0*sqrt(1-omega_m)*t);
}

double xh2solve(double Ti,double xh2,double xhe1,double xhe2,double xhe3,double dti,double NH1,double Nhe1, double Nhe2,double R,double nhi, double *xrate){
   // double nhe1l = (Y/(4*(1-Y)))*nhi*(xhe1);
    double nhe1l = max((Y/(4*(1-Y)))*nhi*(x_min_fix),(Y/(4*(1-Y)))*nhi*(xhe1));
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
    if (C !=C)
    {
        printf("H1 gone wrong\n");
        printf("gamma %0.12e %0.12e %0.12e\n", fac*vol, Gamma(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l),1/(nh1));
        printf("elec coll %0.12e %0.12e\n", GE(Ti),ne);
        exit(0);
    }
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
    if ((1-J*dti)==0)
    {
        printf("x denom is zero\n");
        exit(0);
    }
    return xh2 + dti*(C-xh2*(C+D))/(1-J*dti);
}

double xhe1solve(double Ti,double xh2i,double xhe1,double xhe2,double xhe3,double dti,double NH1,double Nhe1,double Nhe2,double R,double nhi){
        //double nhe1l = (Y/(4*(1-Y)))*nhi*(xhe1);
        double nhe1l = max((Y/(4*(1-Y)))*nhi*(x_min_fix),(Y/(4*(1-Y)))*nhi*(xhe1));
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
        D = betahe1(Ti)*ne*nhe1l + fac*vol*Gammahe1(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l);
  // printf("C %e D %e\n",C,D );
    if ((1+D*dti)==0)
    {
        printf("xhe1 denom is zero\n");
        exit(0);
    }    
   return nhe1l*(xhe1 + C*dti)/(nhe1l+D*dti);
   //     return xhe1 + (C-D*xhe1)*dti;
}


double xhe2solve(double Ti,double xh2i,double xhe1,double xhe2,double xhe3,double dti,double NH1,double Nhe1,double Nhe2,double R,double nhi, double *hrate){
    //double nhe1l = (Y/(4*(1-Y)))*nhi*(xhe1);
    double nhe1l = max((Y/(4*(1-Y)))*nhi*(x_min_fix),(Y/(4*(1-Y)))*nhi*(xhe1));
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
   if (C !=C)
    {
        printf("He2 gone wrong\n");
        printf("gamma %0.12e %0.12e %0.12e\n", fac*vol, Gammahe1(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l),1/(nhe1l));
        printf("elec coll %0.12e %0.12e, %0.12e\n", xhe1,ne,betahe1(Ti));
        exit(0);
    }
   // D = alphahe2(Ti)*ne;
  // printf("C %e D %e\n",C,D );
    *(hrate) = C-D*xhe2;
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
    if (C !=C)
    {
        printf("He3 gone wrong\n");
        printf("gamma %0.12e %0.12e %0.12e\n", fac*vol, Gammahe2(NH1,Nhe1,Nhe2,nh1,nhe1l,nhe2l),xhe2/(nhe2l));
        printf("elec coll %0.12e %0.12e, %0.12e\n", xhe2,ne,betahe2(Ti));
        exit(0);
    }
  //  printf("C %e D %e\n",C,D );
   *(hrate) = C -D*xhe3; 
   if ((1+D*dti)==0)
    {
        printf("xhe3 denom is zero\n");
        exit(0);
    }
   return (xhe3 + C*dti)/(1+D*dti);
   // return xhe3 + (C-D*xhe3)*dti;
}

double Tsolve(double Ti,double NH1,double Nhe1,double Nhe2,double xh2i,double xhe1,double xhe2,double xhe3,double nhi,double dti,double R,double tsim, double *trate, double *xrate){
    
    
    double fac = -alpha*Ng/pow(3e21,3);
    double fh = fac*1.6e-12; // 1.6e-12 is factor ev to erg
    double vol = 1.0/(4*PI*R*R*cellsize);
    
  //  double K = (2.0/3)*(1.0/(nhi*kb)); // only H # cm^3erg^-1 K    
    
    double nh1 = nhi*(1-xh2i);
    //double nhe1l = (Y/(4*(1-Y)))*nhi*(xhe1);
    double nhe1l = max((Y/(4*(1-Y)))*nhi*(x_min_fix),(Y/(4*(1-Y)))*nhi*(xhe1));
    double nhe2l = (Y/(4*(1-Y)))*nhi*xhe2;
    double nhe3l = (Y/(4*(1-Y)))*nhi*xhe3;
    double ne = xh2i*nhi+(Y/(4*(1-Y)))*nhi*(xhe2+2*xhe3);
    double n = nhi + (Y/(4*(1-Y)))*nhi + ne;

    double K = (2.0/3)*(1.0/(n*kb)); // From Chen
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
    if (K*nu==0 || nu!=nu)
    {
        printf("Knu is zero\n");
        printf("nu %0.12e,%0.12e,%0.12e,%0.12e \n",nu,xh2i,xhe2,xhe3 );
        exit(0);
    }


    
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

    C = c1*ne*xh2i*nhi + c11*ne*nhe2l + c12*ne*nhe3l +c2*ne*nh1 + c21*ne*nhe2l + c3*(nhi*xh2i+nhe2l+4*nhe3l)*ne + c4*ne + 2*Hubble(tsim)*Ti/(K*nu) + *(xrate)*(Ti/n)/(K*nu);
    //double C =0;
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
      //  printf("c5 %0.12e\n", c5);
       // printf("c51 %0.12e\n", c51);
       // printf("c52 %0.12e\n", c52);
        printf("H %0.12e\n", 2*Hubble(tsim)*Ti/(K*nu));
        exit(0);
    }   
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
    
    double dLdt = -(A1+ A2 +A3 +A4 +A11+ A12+ A21 +2*Hubble(tsim)/(K*nu) + *(xrate)*(1/n)/(K*nu));//+(1/(nhi*xh2i))*(*xrate));
    //double dLdt =0;
   // printf("%e\n", K)
	//printf("%e %e\n",fh*vol*Heat(NH1,nh1,nhe1,nhe2)-C, dti);
    *(trate) = L*K*nu;
    return Ti+ dti*L*K*nu;
   //return Ti + nu*dti*L*K/(1-nu*dLdt*K*dti) ;
  // return (L*K*nu*dti + Ti)/(1+2*Hubble(tsim)*dti) ; // anninos implicit
}


long xTiter(double dti,double NH1,double Nhe1, double Nhe2,double R,double nhi,double tsim,double ti,double *Tt, double *xh2t,double *xhe1t,double *xhe2t,double *xhe3t,double *dtrec){
  //  double *Tt =(double*)malloc(countl*sizeof(double));
  //  Tt[0] = T[k][j];
  //  double *xh2t =(double*)malloc(countl*sizeof(double)); 	
   // xh2t[0] = nh2[k][j]/nhi;
   // double *NHt =(double*)malloc(countl*sizeof(double)); 	
   // NHt[0] = NH[k][j];
   // printf("%d %lf\n",k,Tt[0]);
    long storei;
    double frac_change = 0.0;
    //double ne_init;
    double dUU;
    double Trate=0;
    double Xrate=0;
    double he2rate = 0;
    double he3rate = 0; 
    double nerate = 0;
    double ne_init,ne_up;
    for (long i=0;i<countl;i++){
        

        //Tt[i+1] = (double) max(T2_MIN_FIX,Tsolve(Tt[i],NH1,Nhe1,Nhe2,xh2t[i],xhe1t[i],xhe2t[i],xhe3t[i],nhi,dti,R,tsim,&Trate));
        Tt[i+1] = Tsolve(Tt[i],NH1,Nhe1,Nhe2,xh2t[i],xhe1t[i],xhe2t[i],xhe3t[i],nhi,dti,R,tsim,&Trate,&nerate);
      //  printf("T %e\n",Tt[i+1] );
        dUU = (double) fabs(max(T2_MIN_FIX,Tt[i]+Trate*dti)-Tt[i]);
     //   printf("heat %e\n",Heat(NH1,nhi,0,0)*(-alpha*Ng/pow(3e21,3))* 1.6e-12*1.0/(4*PI*R*R*cellsize*kb));
 //   printf("%e\n",fabs(Tt[i+1]-Tt[i])/(Tt[i]+T_MIN));
 //   printf("%e\n",Tt[i+1]);

//if(dti>1e-7){
      //  printf("i %d\n", i);
   //     printf("change %e\n",(double) max(dUU,fabs(Tt[i+1]-Tt[i]))/(Tt[i]+T_MIN) );

       // printf("duu %e %e\n",dUU,fabs(Tt[i+1]-Tt[i]));
	if((double) max(dUU,fabs(Tt[i+1]-Tt[i]))/(Tt[i]+T_MIN)>lim & dti >dtmin){
	//printf("hii");
	return 100000000000;
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

     if((double) max(fabs(xh2t[i+1]-xh2t[i]),Xrate*dti)/(xh2t[i]+X_FM) >lim & dti>dtmin){
	return 100000000001;
}
else{
	frac_change = (double) max(frac_change,fabs(xh2t[i+1]-xh2t[i])/xh2t[i]);
}

//xhe1t[i+1] = xhe1solve(Tt[i+1],xh2t[i],xhe2t[i],xhe3t[i],dti,NH1,Nhe1,Nhe2,R,nhi);
xhe1t[i+1] = xhe1solve(Tt[i+1],xh2t[i+1],xhe1t[i],xhe2t[i],xhe3t[i],dti,NH1,Nhe1,Nhe2,R,nhi);



xhe2t[i+1] = xhe2solve(Tt[i+1],xh2t[i+1],xhe1t[i+1],xhe2t[i],xhe3t[i],dti,NH1,Nhe1,Nhe2,R,nhi,&he2rate);



xhe3t[i+1] = xhe3solve(Tt[i+1],xh2t[i+1],xhe1t[i+1],xhe2t[i+1],xhe3t[i],dti,NH1,Nhe1,Nhe2,R,nhi,&he3rate);

//}

/*check if xhe1+2+3 =1*/

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


	//ne_init = xh2t[i]*nhi;
	ne_init = xh2t[i]*nhi+(Y/(4*(1-Y)))*nhi*(xhe2t[i]+2*xhe3t[i]);
    ne_up = xh2t[i+1]*nhi+(Y/(4*(1-Y)))*nhi*(xhe2t[i+1]+2*xhe3t[i+1]);
    nerate =  Xrate*nhi+(Y/(4*(1-Y)))*nhi*(he2rate+2*he3rate);
    if(fabs(ne_up - ne_init)/ne_init>lim & dti>dtmin){
	return 100000000002;
}
else {
	frac_change = (double) max(frac_change,fabs(ne_up - ne_init)/(ne_init+X_FM));
}

	double tcal = ti + dti;
    ti = tcal ;

//if (dti>1e-7){
  

    if(frac_change <0.5){  // has value 0.05 in paper but 0.5 here
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

	printf("%lf %ld\n",Tt[storei],storei);
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
	long out;
	double *Tt =(double*)malloc(countl*sizeof(double));
	double *xh2t =(double*)malloc(countl*sizeof(double));
    double *xhe1t =(double*)malloc(countl*sizeof(double));
    double *xhe2t =(double*)malloc(countl*sizeof(double));
    double *xhe3t =(double*)malloc(countl*sizeof(double));
	double dt_rec;
	for(int k=0; k<ncells;k++){
    printf("dt new %e %d %d\n",dt[k][j],k,j);
    //printf("T[%d] = %lf\n",k,T[k][j] );
    //dt[k][j+1] = dt[k][j];
    printf("cell %d\n",k);
	//mu = pow((X*(1+nh2[k][j]/nh[k])+Y*(1+nh2[k][j]/nh[k]+xhe2+2*xhe3)/4),(-1));
    Tt[0] = T[k][j];
    xh2t[0] = nh2[k][j]/nh[k];
    xhe2t[0] = nhe2[k][j]/((Y/(4*(1-Y)))*nh[k]);
    xhe3t[0] = nhe3[k][j]/((Y/(4*(1-Y)))*nh[k]);
    xhe1t[0] = 1-xhe2t[0] -xhe3t[0];
	out = xTiter(dt[k][j],NH[k][j],NHe1[k][j],NHe2[k][j],(1+k)*cellsize,nh[k],tsim,ti,Tt,xh2t,xhe1t,xhe2t,xhe3t,&dt_rec);  // equivalent of cool_step routine
	//printf("%d\n",out);
	

    while((out==100000000000 || out==100000000001 || out ==100000000002)){
	dt[k][j] = dt[k][j]/2;
    dt[k][j] = max(dt[k][j],dtmin);
	printf("%ld %e\n",out-100000000000, dt[k][j]);
    /*if (dt[k][j]<1e-8)
    {
        out = xTiter(dt[k][j],NH[k][j],(1+k)*cellsize,nh[k],tsim,ti,Tt,xh2t,&dt_rec);
        break;
    }*/
	out = xTiter(dt[k][j],NH[k][j],NHe1[k][j],NHe2[k][j],(1+k)*cellsize,nh[k],tsim,ti,Tt,xh2t,xhe1t,xhe2t,xhe3t,&dt_rec);

	//printf("%d\n",out);
} 
	printf("%ld\n",out);
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
	double  omega_m = 0.32;
    double  H0 = 67.32*1.05e-6;
    double z =  7;    
    double a = 1/(1+z);
    double ti = (1.0/H0)*(2.0/3)*pow(1-omega_m,-0.5)*asinh(pow((1-omega_m)/omega_m,0.5)*pow(a,1.5));
    printf("ti %lf\n",ti );
	//double ti = 935.7;
	double tf = ti + atof(argv[1]);
//	printf("%lf\n",tf);
	//double epsilon = 5e-2;	
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
	nh[i] = 27.24e-6; // density at z = 7 as calculated by me
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
	//printf("%lf\n",NH[0][0]);
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


FILE *fptr = fopen("dxh1_1myr.txt","w");
FILE *fdptr = fopen("dxhe2_1myr.txt","w");
FILE *fd1ptr = fopen("dxhe3_1myr.txt","w");
FILE *fd2ptr = fopen("dT_1myr.txt","w");
FILE *fd3ptr = fopen("dxhe1_1myr.txt","w");
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

