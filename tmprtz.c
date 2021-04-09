#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h> 
#include<time.h>

#ifndef PI
#define PI 3.1415
#endif
#define cellsize 10
#define boxsize 5000
#define count 100000
#define countl 10000000
#define Ng 3.15e13*1e57
#define alpha -1.5
#define kb 1.38e-16

double cumsum(double *nH, double **nH2, int i, int j){
	double sum = 0;
	for(int l =0;l<i;l++){
	sum = sum + nH[i] -nH2[i][j];
}
return sum;
}

double GE(double Ti){
    double Th1 = 157807; 
    double lam = 2*Th1/Ti;
    
    return 3.15e13*21.11*pow(Ti,-3/2)*exp(-lam/2)*pow(lam,-1.089)/pow(1+pow((lam/0.354),0.874),1.101) ;
}

double ralpha(double Ti){
    double Th1 = 157807; 
    double lam = 2*Th1/Ti;

    return 4*pow(lam,1.503)/pow((1+pow(lam/0.522,0.47)),1.923);
}

double sigma(double f){
    double sigma_0=5.475 ;
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

double PH1(double x,double nh1,double nhe1,double nhe2){
	double qH = exp(-sigma(x)*nh1*cellsize*3e7);
	double qHe1 = 1;	//np.exp(-sigmahe1(x)*nhe1*dr*3e7);
double 	qHe2 = 1;	//np.exp(-sigmahe2(x)*nhe2*dr*3e7)	;
double	tau_tot = sigma(x)*nh1*cellsize*3e7 ;//+sigmahe1(x)*nhe1*dr*3e7+ sigmahe2(x)*nhe2*dr*3e7;
double	D = (1-qH)*qHe2*qHe1 +(1-qHe1)*qH*qHe2 + (1-qHe2)*qH*qHe1;

	return (1-qH)*qHe2*qHe1*(1-exp(-tau_tot))/D;
}

double Gamma(double NH,double nh1,double nhe1,double nhe2){
    
    double u[] ={13.6, 14.4059, 15.2595, 16.1636, 17.1214, 18.1359, 19.2105, 20.3488,21.5545, 22.8317, 24.1846, 25.6176, 27.1356, 28.7435, 30.4466, 32.2507, 34.1617, 36.1859, 38.33, 40.6012, 43.007, 45.5553, 48.2546, 51.1139, 54.1426,57.3507, 60.749, 64.3486, 68.1615, 72.2003, 76.4784, 81.0101, 85.8102, 90.8948, 96.2806, 101.986, 108.029, 114.43, 121.21, 128.392, 136};
    // SPACING = (b/a)^(1/n) , n= 40, b = 136, a =13.6 
    double fsum = 0;
    double x0 = u[0];
    double df = 1;	    
    for(int i=0;i<sizeof(u)/sizeof(u[0]);i++){
     fsum = fsum + df*x0*pow(u[i]/x0,alpha)*exp(-sigma(u[i])*NH)*PH1(u[i],nh1,nhe1,nhe2)/u[i] ;	
}
   

    return fsum;
}

double Heat(double NH,double nh1,double nhe1,double nhe2){
    
    double u[] ={13.6, 14.4059, 15.2595, 16.1636, 17.1214, 18.1359, 19.2105, 20.3488,21.5545, 22.8317, 24.1846, 25.6176, 27.1356, 28.7435, 30.4466, 32.2507, 34.1617, 36.1859, 38.33, 40.6012, 43.007, 45.5553, 48.2546, 51.1139, 54.1426,57.3507, 60.749, 64.3486, 68.1615, 72.2003, 76.4784, 81.0101, 85.8102, 90.8948, 96.2806, 101.986, 108.029, 114.43, 121.21, 128.392, 136};
    // SPACING = (b/a)^(1/n) , n= 40, b = 136, a =13.6 
    double fsum = 0;
    double x0 = u[0];
    double df = 1;	    
    for(int i=0;i<sizeof(u)/sizeof(u[0]);i++){
     fsum = fsum + df*(u[i]-x0)*x0*pow(u[i]/x0,alpha)*exp(-sigma(u[i])*NH)*PH1(u[i],nh1,nhe1,nhe2)/u[i] ;	
}
   

    return fsum;
}

double Hubble(double t){
    double  omega_m = 0.32;
    double  H0 = 67.32*1.05e-6;
    
    return sqrt(1-omega_m)*H0*cosh(1.5*H0*sqrt(1-omega_m)*t)/sinh(1.5*H0*sqrt(1-omega_m)*t);
}

double xh2solve(double Ti,double ne,double xh2,double dt,double NH1,double R,double nhi){
    double nhe2 =0;
    double nhe3 =0;
    double nhe1 =0;
    double xhe2,xhe3;
    xhe2 =0;
    xhe3= 0;
    double fac = Ng/pow(3e21,3);
    double vol = 1.0/(4*PI*R*R*cellsize);
    double nh1 = nhi*(1-xh2);
    double C,D,dbdt,dalphadt,nu,Tnu,dCdx,dDdx,J;
    
    C = GE(Ti)*ne + fac*vol*Gamma(NH1,nh1,nhe1,nhe2)/(nhi*xh2);
    D = ralpha(Ti)*ne;
    //printf("%lf,%lf\n",C,D);
    double X = 1;
    double Y = 0;
   
    dbdt = (1.07732*1e14*exp(-157807/Ti))/(pow((1 + 158648*pow((1/Ti),0.874)),1.101) *pow((1/Ti),1.089) *pow(Ti,(7/2)))
- (1.02402*1e9*exp(-157807/Ti))/(pow((1 + 158648*pow((1/Ti),0.874)),1.101)*pow((1/Ti),1.089)*pow(Ti,(5/2))) 
+ (7.43439*1e8*exp(-157807/Ti))/(pow((1 + 158648*pow((1/Ti),0.874)),1.101)*pow((1/Ti),0.089)*pow(Ti,(3/2))) 
+ (1.0422*1e14*exp(-157807/Ti)*pow((1/Ti),0.785))/(pow((1 + 158648*pow((1/Ti),0.874)),2.101)*pow(Ti,(3/2)));
    
    dalphadt = -((1.10726*1e9*pow((1/Ti),2.503))/(1 + 521.553*pow((1/Ti),0.47))*1.923) + (
 3.4727*1e11*pow((1/Ti),2.973))/pow((1 + 521.553*pow((1/Ti),0.47)),2.923);

    nu = pow((X*(1+xh2)+Y*(1+xh2+xhe2+2*xhe3)/4),(-1));
    Tnu = Ti/nu;
    
    dCdx = nhi*GE(Ti) - ne*Tnu*pow(nu,2)*X*dbdt ;
    dDdx = nhi*ralpha(Ti) -ne*Tnu*pow(nu,2)*X*dalphadt;
    J = dCdx -(C+D) - xh2*(dCdx+dDdx);
    
    return xh2 + dt*(C-xh2*(C+D))/(1-J*dt);
}


double xh2iter(double Ti,double ne,double xh2i,double dt,double NH1,double R,double nhi,double tsim,double ti){
    double *xh2t =(double*)malloc(countl*sizeof(double)); 	
    xh2t[0] = xh2i;
    int storei ;	
    for (int i=0;i<countl;i++){
        xh2t[i+1] = xh2solve(Ti,ne,xh2t[i],dt,NH1,R,nhi);
	//printf("%lf\n",xh2t[i+1]);
        /*
        while (fabs(xh2t[i+1]-xh2t[i])/xh2t[i]>0.1 && dt>1e-5){
            dt = dt/2;
           
            xh2t[i+1] = xh2solve(Ti,ne,xh2t[i],dt,NH1,R,nhi);
            
            printf("%e\n",dt);
            }
	*/
        double tcal = ti+ dt;
        ti = tcal;
        if (tcal >= tsim){
	    storei = i+1;
            break;
	}	
}
	//printf("%lf\n",xh2t[storei]);	
    return xh2t[storei];
}


double Tsolve(double Ti,double NH1,double xh2i,double nhi,double dt,double R,double tsim){
    
    
    double fac = Ng/pow(3e21,3);
    double fh = fac*1.6e-12; // 1.6e-12 is factor ev to erg
    double vol = 1.0/(4*PI*R*R*cellsize);
    double K = (2.0/3)*(1/(nhi+nhi*xh2i)); // only H # cm^3erg^-1 K    
    double nh1 = nhi*(1-xh2i);
    double nhe1 = 0;
    double nhe2 = 0;
    double T5 = Ti/100000.0;
    
    double f1 = 1.778e-29*3.15e13;	//#*6.24e11 # in erg cm^3 Myr^-1 K^-1 #4.05   #
    double f2 = 7.5e-19*3.15e13;	//#*6.24e11 # in erg cm^3 Myr^-1 #1.711e+11
    double f3 = 1.42e-27*3.15e13;	//#*6.24e11# in erg cm^3 K^-1 #489.62
    double gff = 1.5 ;	//# from cen #1.1+ 0.34*(np.exp(-5.5+np.log(T)))**2
    double f4 = 5.65e-36*3.15e13;	//#*6.24e11 #in erg cm^3 Myr^-1 K^-1 #1.289e-6
/*
    f11 = 1
    f12 =8*1.778*1e-29*3.15e13/(1.38e-16)
    f21 =5.54e-17*3.15e13/(1.38e-16)
    f22 =9.1e-27*3.15e13/(1.38e-16)
*/
    double Th1 = 157807;
    //The = 631515
    double lamb = 2*Th1/Ti;
    //lamb2 = 2*The/T
    double omega_m = 0.32;
    double H0 = 67.32*1.05e-6;
    double rho = 1.0;
    //#fh = 1.04e+8
    double Hf = 1.05e-6;
    
    double z = 5.945518;
    
    
    double c1 = (f1)*Ti*pow(lamb,1.965)/pow((1.0+pow((lamb/0.541),0.502)),2.697) ;// recombination cooling coeffecient
    double c2 = (f2)*pow((1+pow(T5,1.0/2)),-1)*exp(-118355/Ti);  //collisional excitation cooling coeffecient
    double c3 = gff*(f3)*pow(Ti,0.5); // Brehmstrallung/ free free emission
    double c4 = (f4)*(Ti-2.73*(1+z))*pow(1+z,4);
    
   
    //double C =  c1*pow((xh2i*nhi),2) +c2*xh2i*nhi*(nhi-xh2i*nhi)+c3*pow((nhi*xh2i),2) +c4*(nhi*xh2i) + 2*Hubble(tsim)*Ti*kb/K;
    double C =0;
       
    double L = (fh*vol*Heat(NH1,nh1,nhe1,nhe2)/kb - C/kb);
    
     //printf("%e, %e\n",c1, c2);
    //derivatives of cooling processes 
    
    double A1 =  -(6.1712*1e10*f1*pow((1/Ti),1.965))/pow((1 + 784.353*pow((1/Ti),0.502)),2.697) + (
 6.79107*1e13*f1*pow((1/Ti),2.467))/pow((1 + 784.353*pow((1/Ti),0.502)),3.697); 
    double A2 = (118355*exp(-118355/Ti)*f2)/((1 + 0.00316228*pow(Ti,0.5))*pow(Ti,2)) - (
 0.00158114*exp(-118355/Ti)*f2)/(pow((1 + 0.00316228*pow(Ti,0.5)),2)*pow(Ti,0.5)); 
    double A3 = (0.5*f3*gff)/pow(Ti,0.5); 
    double A4 =   f4*pow((1+z),4);
    
    
    
    double X = 1;
    double Y = 0;
    double xhe2 = 0; 
    double xhe3 = 0;
    double nu = pow((X*(1+xh2i)+Y*(1+xh2i+xhe2+2*xhe3)/4),(-1));
    //double dLdt = (A1/kb + A2/kb +A3/kb +A4/kb);
    double dLdt =0;
   // printf("%e %e\n", K);
    return Ti + dt*L*K/(1-dLdt*K*dt) ;
}

double Titer(double Ti,double xh2i,double dt,double NH1,double R,double nhi,double tsim,double ti){
    double *Tt =(double*)malloc(countl*sizeof(double));
    Tt[0] = Ti;
    
    int storei;
    for (int i=0;i<countl;i++){
  
        Tt[i+1] = Tsolve(Tt[i],NH1,xh2i,nhi,dt,R,tsim);
      // printf("%lf\n",Tt[i+1]);
     /*   
        while (fabs(Tt[i+1]-Tt[i])/Tt[i]>0.1){ //and np.abs(K*L*dt)/Tt[i]>0.1:
            dt = dt/2;
            Tt[i+1] = Tsolve(Tt[i],NH1,xh2i,nhi,dt,R,tsim);
            //K = Tsolve(Tt[i],NH1,xh2t,nh,dt,R,tsim)[1]
            //L = Tsolve(Tt[i],NH1,xh2t,nh,dt,R,tsim)[2]
            printf("%lf\n",dt);
    }
*/
        //#print(Tt[i+1])
	       
	double tcal = ti + dt;
        ti = tcal ;
        if (tcal >= tsim){
            storei = i+1;
            break ;
} 
   }
	printf("%lf\n",Tt[storei]);
    return Tt[storei];
}

void evolve_eqns(double *nh,double **nh2,double **NH,double **T,int ncells,double ti,double tsim,double Dt, int j){
	double dt = Dt*1e-4;

	for(int k=0; k<ncells;k++){
	nh2[k][j+1] = nh[k]*xh2iter(T[k][j],nh2[k][j],nh2[k][j]/nh[k],dt,NH[k][j],(1+k)*cellsize,nh[k],tsim,ti);	
	//printf("%e\n",nh2[k][j+1]);	
	//*(*(nh2+k)+j+1) = 0.1*(*(*(nh2+k)+j));
	//T[k][j+1] = T[k][j];
	T[k][j+1] = Titer(T[k][j],nh2[k][j]/nh[k],dt,NH[k][j],(1+k)*cellsize,nh[k],tsim,ti);
}
	for(int k=0;k<ncells;k++){
	NH[k][j+1] = cumsum(nh,nh2,k,j+1)*cellsize*3e7;
}
		
}


int main(int argc, char **argv)
{

	/*For each position, solve euler to get nh(xpos,t) and T(xpos,t) where t is tsim= time after each global dt */	
	clock_t t;
    	t = clock();
    	
    
	
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
	double **NH = (double **)malloc(r * sizeof(int *));
    	for (int i=0; i<r; i++)
		NH[i] = (double *)malloc(c * sizeof(int));
	double **T = (double **)malloc(r * sizeof(int *));
    	for (int i=0; i<r; i++)
         	T[i] = (double *)malloc(c * sizeof(int));	

/* Initial values*/
	for (int i = 0; i < ncells; ++i){
	//printf("%d\n",i);
	nh[i] = 7e-5;
	nh2[i][0] = 0.05*nh[i];
	T[i][0] = 20;
	NH[i][0] = cumsum(nh,nh2,i,0)*cellsize*3e7;	
	//printf("%lf\n",NH[i][0]);
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
for (int j=0;j<count;j++){
	tsim = ti +Dt;
	//printf("%d\n",j);
	evolve_eqns(nh,nh2,NH,T,ncells,ti,tsim,Dt,j);
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

printf("%e,",nh2[i][storej]/nh[i]);
}
printf("]");

t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
  
    printf("fun() took %f seconds to execute \n", time_taken);	    


 FILE *fptr = fopen("cwTxh1_01myr.txt","w");
 FILE *fdptr = fopen("cwTT_01myr.txt","w");
     	for (int l = 0; l < ncells; ++l)
    {
   
    	fprintf(fptr,"%lf %0.12e \n", (double)l*cellsize+cellsize,(nh[l]-nh2[l][storej])/nh[l]);
	fprintf(fdptr,"%0.12e \n",T[l][storej]);
    } 
	

}

