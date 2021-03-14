
import numpy as np
from scipy import integrate
#from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math as m
import sys
from scipy.integrate import solve_ivp

p_inf = float("inf")
alpha = -1.5   #source spectral index

# Electron collisional ionisation coeff. HI (cm^3 Myr^-1) (Hui & Gnedin 1997)
def GE(T): 
    Th1 = 157807 
    lam = 2*Th1/T
    
    return np.nan_to_num(3.15e13*21.11*np.power(T,-3/2)*np.exp(-lam/2)*np.power(lam,-1.089)/np.power(1+np.power((lam/0.354),0.874),1.101) )

# Electron collisional ionisation coeff. HeI (cm^3 Myr^-1) (Hui & Gnedin 1997)
def GE2(T):
    Th1 = 285335
    lam = 2*Th1/T

    return  np.nan_to_num(3.15e13*32.38*np.power(T,-3/2)*np.exp(-lam/2)*np.power(lam,-1.146)/np.power(1+np.power((lam/0.416),0.987),1.056) )

# Electron collisional ionisation coeff. HeII (cm^3 Myr^-1) (Hui & Gnedin 1997)
def GE3(T):
    Th1 = 631515 
    lam = 2*Th1/T

    return np.nan_to_num(3.15e13*19.95*np.power(T,-3/2)*np.exp(-lam/2)*np.power(lam,-1.089)/np.power(1+np.power((lam/0.553),0.735),1.275))


# HII recombination coefficient (cm^3 Myr^-1) (Hui & Gnedin 1997)
def ralpha(T):
    Th1 = 157807 
    lam = 2*Th1/T

    return np.nan_to_num(4*np.power(lam,1.503)/np.power((1+np.power(lam/0.522,0.47)),1.923) )

# HeII recombination coefficient (cm^3 Myr^-1) (Hui & Gnedin 1997)
def alpha2(T):
    Th1 = 285335 
    lam = 2*Th1/T
    f = 3e-14*3.15e13
    return np.nan_to_num(f*np.power(lam,0.654))

# HeIII recombination coefficient (cm^3 Myr^-1) (Hui & Gnedin 1997)
def alpha3(T):
    Th1 = 631515 
    lam = 2*Th1/T
    f =2*1.269e-13*3.15e13
    return np.nan_to_num(f*np.power(lam,1.503)/np.power((1.0+np.power(lam/0.522,0.47)),1.923) )


# HI photoionisation cross-section (nm^2) (Hui & Gnedin 1997)
def sigma(f):
    sarr = np.array([6.34629636e-04, 5.43455029e-04, 4.65000442e-04, 3.97554467e-04,
 3.39612511e-04, 2.89889481e-04, 2.47254937e-04, 2.10728637e-04,
 1.79463914e-04, 1.52722910e-04, 1.29870710e-04, 1.10358716e-04,
 9.37102876e-05, 7.95176357e-05, 6.74279672e-05, 5.71366171e-05,
 4.83831488e-05, 4.09433976e-05, 3.46248308e-05, 2.92621389e-05,
 2.47140484e-05, 2.08596489e-05, 1.75953276e-05, 1.48325999e-05,
 1.24960412e-05])
    return sarr
    #print(sigma_0*(np.power((x-1),2)+np.power(yw,2))*np.power(y,0.5*P-5.5)/np.power(1+np.sqrt(y/ya),P))
    #return sigma_0*(np.power((x-1),2)+np.power(yw,2))*np.power(y,0.5*P-5.5)/np.power(1+np.sqrt(y/ya),P)

# HeI photoionisation cross-section (nm^2) (Hui & Gnedin 1997)
def sigmahe1(f):
    """
    sigma_0=9.492e-2 
    nu_0 = 13.61   
    y0 = 0.4434
    x = f/nu_0 - y0
    yw = 2.039
    P = 3.188
    ya = 1.469 
    y1  = 2.136 
    y = np.sqrt(x*x+y1*y1)

    return sigma_0*(np.power((x-1),2)+np.power(yw,2))*np.power(y,0.5*P-5.5)/np.power(1+np.sqrt(y/ya),P)
    """
    sharr = np.array([0.00016854685216437885 ,0.000148789993158018 ,0.00013099114911427456 ,0.00011499931869256949 ,0.00010067429077457867 ,8.788329535309154e-05 ,7.650145649841567e-05 ,6.640961292114752e-05 ,5.7493443719878704e-05 ,4.9644676222502684e-05 ,4.275981898138132e-05 ,3.674123858572162e-05 ,3.1497703315154895e-05 ,2.6944209210561267e-05 ,2.3001555455443907e-05 ,1.9597505097898313e-05 ,1.66668696411452e-05 ,1.4150044063358387e-05 ,1.1993999013529437e-05 ,1.0151026825999545e-05 ,8.578849636350473e-06 ,7.2404561650563485e-06 ,6.103149975955958e-06 ,5.138362304561513e-06 ,4.321235699091948e-06])
    return sharr

# HeII photoionisation cross-section (nm^2) (Hui & Gnedin 1997)
def sigmahe2(f):
    """
    sigma_0=1.369 
    nu_0 = 1.72   
    y0 = 0
    x = f/nu_0 - y0
    yw = 0
    P = 2.963
    ya = 32.88 
    y1  = 0 
    y = np.sqrt(x*x+y1*y1)
    return sigma_0*(np.power((x-1),2)+np.power(yw,2))*np.power(y,0.5*P-5.5)/np.power(1+np.sqrt(y/ya),P)
    """
    return np.array([0.00015810273336613407 ,0.0001353875222712522 ,0.00011584018753273866 ,9.90347370024959e-05 ,8.460007463839189e-05 ,7.221221856716023e-05 ,6.15904883675312e-05 ,5.249107238281784e-05 ,4.470210275702134e-05 ,3.804072002671639e-05 ,3.234820283585642e-05 ,2.7487482060700154e-05 ,2.3340503712192593e-05 ,1.9805406394181356e-05 ,1.6793925421831326e-05 ,1.4230370938590393e-05 ,1.2050063960375295e-05 ,1.0196930671508321e-05 ,8.623190146281607e-06 ,7.287585093569676e-06 ,6.154799152296611e-06 ,5.1948130461697165e-06 ,4.381815551758566e-06 ,3.6937627146712386e-06 ,3.111846405602817e-06 ])

# Probability of absorption by HI in the current cell
def PH1(x,nh1,nhe1,nhe2):
	qH = np.exp(-sigma(x)*nh1*dr*3e7)
	qHe1 = np.exp(-sigmahe1(x)*nhe1*dr*3e7)
	qHe2 = np.exp(-sigmahe2(x)*nhe2*dr*3e7)	
	tau_tot = sigma(x)*nh1*dr*3e7 +sigmahe1(x)*nhe1*dr*3e7+ sigmahe2(x)*nhe2*dr*3e7
	D = (1-qH)*qHe2*qHe1 +(1-qHe1)*qH*qHe2 + (1-qHe2)*qH*qHe1
	#print(D)
	return (1-qH)*qHe2*qHe1*(1-np.exp(-tau_tot))/D

# Probability of absorption by HeI in the current cell
def PHe1(x,nh1,nhe1,nhe2):
	qH = np.exp(-sigma(x)*nh1*dr*3e7)
	qHe1 = np.exp(-sigmahe1(x)*nhe1*dr*3e7)
	qHe2 = np.exp(-sigmahe2(x)*nhe2*dr*3e7)
	tau_tot = sigma(x)*nh1*dr*3e7 +sigmahe1(x)*nhe1*dr*3e7+ sigmahe2(x)*nhe2*dr*3e7
	D = (1-qH)*qHe2*qHe1 +(1-qHe1)*qH*qHe2 + (1-qHe2)*qH*qHe1
	return (qH)*qHe2*(1-qHe1)*(1-np.exp(-tau_tot))/D

# Probability of absorption by HeII in the current cell
def PHe2(x,nh1,nhe1,nhe2):
	qH = np.exp(-sigma(x)*nh1*dr*3e7)
	qHe1 = np.exp(-sigmahe1(x)*nhe1*dr*3e7)
	qHe2 = np.exp(-sigmahe2(x)*nhe2*dr*3e7)
	tau_tot = sigma(x)*nh1*dr*3e7 +sigmahe1(x)*nhe1*dr*3e7+ sigmahe2(x)*nhe2*dr*3e7
	D = (1-qH)*qHe2*qHe1 +(1-qHe1)*qH*qHe2 + (1-qHe2)*qH*qHe1
	return (qH)*(1-qHe2)*qHe1*(1-np.exp(-tau_tot))/D

# Functions returning effective Luminosity for HI/HeI/HeII in the current cell (dimensionless)
def L(x,NH,nh1,nhe1,nhe2):
    #Bolton isothermal
    #alpha = 0 
    x0 = 13.6
    np.seterr(all='ignore')
    #return np.power(x/x0,alpha)*np.exp(-sigma(x)*NH)*(1-np.exp(-nh1*dr*sigma(x)*3e7))/x 
    h = x0*np.power(x/x0,alpha)*np.exp(-sigma(x)*NH)*PH1(x,nh1,nhe1,nhe2)/x 
    return h, (x-x0)*h
def Lhe1(x,NHe1,nh1,nhe1,nhe2):
    #alpha =-1.5
    #Bolton isothermal
    #alpha = 0 
    x0 = 54.5#24.6
    np.seterr(all='ignore')
    #return np.power(x/x0,alpha)*np.exp(-sigma(x)*NH)*(1-np.exp(-nh1*dr*sigma(x)*3e7))/x
    h = x0*np.power(x/x0,alpha)*np.exp(-sigmahe1(x)*NHe1)*PHe1(x,nh1,nhe1,nhe2)/x  
    return h, (x-x0)*h
def Lhe2(x,NHe2,nh1,nhe1,nhe2):
    #alpha =-1.5
    #Bolton isothermal
    #alpha = 0 
    x0 = 54.5
    np.seterr(all='ignore')
    #return np.power(x/x0,alpha)*np.exp(-sigma(x)*NH)*(1-np.exp(-nh1*dr*sigma(x)*3e7))/x 
    p = x0*np.power(x/x0,alpha)*np.exp(-sigmahe2(x)*NHe2)*PHe2(x,nh1,nhe1,nhe2)/x
    return p, (x-x0)*p

# Functions calculating photoionisation rates and photoionisation heating rates (dimensionless)
def Gamma(NH,nh1,nhe1,nhe2):
    #alpha = -1.5
    #x0 = 13.6
    #L = lambda x:np.power(x/x0,alpha)*np.exp(-sigma(x)*NH)/x
    #fsum = integrate.quad(L,13.6,500) 
    #u = np.array([13.6,27.2,54.4,108.8,217.6])
    u =np.array([13.6, 14.4059, 15.2595, 16.1636, 17.1214, 18.1359, 19.2105, 20.3488,21.5545, 22.8317, 24.1846, 25.6176, 27.1356, 28.7435, 30.4466, 32.2507, 34.1617, 36.1859, 38.33, 40.6012, 43.007, 45.5553, 48.2546, 51.1139, 54.1426 ])#57.3507, 60.749, 64.3486, 68.1615, 72.2003, 76.4784, 81.0101, 85.8102, 90.8948, 96.2806, 101.986, 108.029, 114.43, 121.21, 128.392, 136])
    # SPACING = (b/a)^(1/n) , n= 40, b = 136, a =13.6 
    y1,y2 = L(u,NH,nh1,nhe1,nhe2)  
    fsum = integrate.simps(y1,u)    #	Integral/Summation
    fsumh = integrate.simps(y2,u)
    np.seterr(all='warn')
    return fsum, fsumh
def Gamma2(NHe1,nh1,nhe1,nhe2):
    #alpha = -1.5
    #x0 = 13.6
    #L = lambda x:np.power(x/x0,alpha)*np.exp(-sigma(x)*NH)/x
    #fsum = integrate.quad(L,13.6,500) 
    #u = np.array([13.6,27.2,54.4,108.8,217.6])
    
    u = np.array([54.5,57.7293, 61.15, 64.7734, 68.6114, 72.6769, 76.9833, 81.5448, 86.3767, 91.4948, 96.9162, 102.659, 108.742, 115.185, 122.01, 129.24, 136.898, 145.01, 153.602, 162.703, 172.344, 182.556, 193.373, 204.831, 216.968])
    y1,y2 = Lhe1(u,NHe1,nh1,nhe1,nhe2) 
    fsum = integrate.simps(y1,u)    #	Integral
    fsumh = integrate.simps(y2,u)
    np.seterr(all='warn')
    return fsum, fsumh
def Gamma3(NHe2,nh1,nhe1,nhe2):
    #alpha = -1.5
    #x0 = 13.6
    #L = lambda x:np.power(x/x0,alpha)*np.exp(-sigma(x)*NH)/x
    #fsum = integrate.quad(L,13.6,500) 
    #u = np.array([13.6,27.2,54.4,108.8,217.6])
    #u =np.array([13.6, 14.4059, 15.2595, 16.1636, 17.1214, 18.1359, 19.2105, 20.3488,21.5545, 22.8317, 24.1846, 25.6176, 27.1356, 28.7435, 30.4466, 32.2507, 34.1617, 36.1859, 38.33, 40.6012, 43.007, 45.5553, 48.2546, 51.1139, 54.1426 ])#57.3507, 60.749, 64.3486, 68.1615, 72.2003, 76.4784, 81.0101, 85.8102, 90.8948, 96.2806, 101.986, 108.029, 114.43, 121.21, 128.392, 136])
    
    u = np.array([54.5,57.7293, 61.15, 64.7734, 68.6114, 72.6769, 76.9833, 81.5448, 86.3767, 91.4948, 96.9162, 102.659, 108.742, 115.185, 122.01, 129.24, 136.898, 145.01, 153.602, 162.703, 172.344, 182.556, 193.373, 204.831, 216.968])
    y1,y2 = Lhe2(u,NHe2,nh1,nhe1,nhe2) 
    fsum = integrate.simps(y1,u)    #	Integral
    fsumh = integrate.simps(y2,u)
    np.seterr(all='warn')
    return fsum, fsumh

# Cooling processes - function returning Lambda_cool
def Cool(T,nH,nH2,nHe2,nHe3,nE):
    T5 = T/100000
    f1=4.05
    f2 =1.711e+11
    f3 =489.62
    f4 = 1.289e-6
    f11 = 1
    f12 =8*1.778*1e-29*3.15e13/(1.38e-16)
    f21 =5.54e-17*3.15e13/(1.38e-16)
    f22 =9.1e-27*3.15e13/(1.38e-16)
    Th1 = 157807
    The = 631515
    lamb = 2*Th1/T
    lamb2 = 2*The/T
    omega_m =0.32
    H0 =67.32*1.05e-6
    rho =1
    fh = 1.04e+8
    Hf = 1.05e-6
    
    z = 5.945518
    nH1= nH -nH2
    nHe1 = Y/(4*(1-Y))*nH - nHe2 - nHe3
    # Source of coefficients: c4-Bolton 2007,{c1,c2,c3,c11,c12,c21}- Hui and Gnedin 1997, c22 - Cen 1992
    c1 = f1*T*np.power(lamb,1.965)/np.power((1.0+np.power((lamb/0.541),0.502)),2.697)		# recombination cooling coefficient HI 
    c2 = f2*np.power((1+np.power(T5,1.0/2)),-1)*np.exp(-118348/T)		# collisional excitation cooling coefficient HI
    c3 = f3*np.power(T,0.5)		# Brehmstrallung/ free free emission
    c4 = f4*(T-2.73*(1+z))*np.power(1+z,4)		# inverse compton cooling 
    c11 = f11*T*alpha2(T)		#recomb HeII 
    c12 = f12*T*np.power(lamb2,1.965)/np.power((1.0+np.power((lamb2/0.541),0.502)),2.697)		#recomb HeIII
    c21 = f21*np.power(T,-0.397)*np.power((1+np.power(T5,1.0/2)),-1)*np.exp(-473638/T)			# collisional excitation cooling coefficient HeI--change to HeII
    c22 = f22*np.power(T,-0.1687)*np.power((1+np.power(T5,1.0/2)),-1)*np.exp(-13179/T)		# collisional excitation cooling coefficient HeII--change to HeI
    #Ltot = c1*nH2*nH2 +c2*nH2*nH+c3*nH2*nH2+c4*nH2		# Hydrogen only
    Ltot = c1*nH2*nE +c2*nE*nH1+c3*(nH2+nHe2+4*nHe3)*nE+c4*nE +c11*nHe2*nE +c12*nHe3*nE+ c21*nE*nHe1+c22*nE*nHe2  # last term assuming typo in Cen1992 paper as units don't match otherwise                                                 
    np.seterr(all='warn')
    return np.nan_to_num(Ltot)

# Hubble parameter (Bolton 2007) verify
def Hubble(t):
    omega_m = 0.32
    H0 = 67.32*1.05e-6
    
    return np.sqrt(1-omega_m)*H0*np.cosh(1.5*H0*np.sqrt(1-omega_m)*t)/np.sinh(1.5*H0*np.sqrt(1-omega_m)*t)


# Rate equations
def model(t,y,NH1,NHe1,NHe2,R,nh):
    Ng =3.15e13*1e57#*7.7e54		# emitted no of photons (Myr^-1)
    #Bolton 
    #Ng = 3.15e13*5e53
    fac = Ng/pow(3e21,3)
    fh = fac*1.6e-12/(1.38e-16) #convert hv from ev to erg, kB in ergK^-1
    dydx =np.zeros(4)   
    Y = 0.24		# cosmological He fraction
    vol = 1/(4*np.pi*R*R*dr)
    #y[0] = z0[0]
    #y[1] = z0[1]
    #y[2] = z0[2]
    #y[3] = z0[3]
    # 0--H2, 1-- He2, 2--He3, 3--T
    
    # Closing conditions
    nh1 =   nh -y[0]
    nhe1 = Y/(4*(1-Y))*nh - y[1] - y[2]
    ne = y[0] + y[1] + 2*y[2]# helium + y[2]
    T = y[3]
    #print(nhe1)
    #print(T)
    
    #H only
    #dydx[0] = np.nan_to_num(nh1*ne*GE(T) -ne*y[0]*ralpha(T) + fac*vol*Gamma(NH1,nh1,nhe1,y[1]))
    #dydx[1] = 0
    #dydx[1] = (2.0/3)*(fh*vol*Heat(NH1,nh1)-Cool(T,nh,y[0],ne))/nh -2*Hubble(t)*T -T/(nh+ne)*dydx[0]
    #print(nh1,nhe1,y[1])
    
	
    # original and complete eqns with H and He
    
    #dydx[0] = np.nan_to_num((nh -y[0])*(y[0] + y[1] + 2*y[2])*GE(y[3]) - (y[0] + y[1] + 2*y[2])*y[0]*ralpha(y[3]) + fac*vol*Gamma(NH1,(nh -y[0]),Y/(4*(1-Y))*nh - y[1] - y[2],y[1]))
    #dydx[1] = np.nan_to_num(fac*vol*Gamma2(NHe1,(nh -y[0]),Y/(4*(1-Y))*nh - y[1] - y[2],y[1])+(Y/(4*(1-Y))*nh - y[1] - y[2])*(y[0] + y[1] + 2*y[2])*GE2(y[3])+y[2]*(y[0] + y[1] + 2*y[2])*alpha3(y[3]) - fac*vol*Gamma3(NHe2,nh -y[0],Y/(4*(1-Y))*nh - y[1] - y[2],y[1])-y[1]*(y[0] + y[1] + 2*y[2])*GE3(y[3])-y[1]*(y[0] + y[1] + 2*y[2])*alpha2(y[3]))   # He II equation
    #dydx[2] = np.nan_to_num(fac*vol*Gamma3(NHe2,(nh -y[0]),Y/(4*(1-Y))*nh - y[1] - y[2],y[1])+y[1]*(y[0] + y[1] + 2*y[2])*GE3(y[3])-y[2]*(y[0] + y[1] + 2*y[2])*alpha3(y[3])) # HeIII abundance
    
    dydx[0] = np.nan_to_num(nh1*ne*GE(T)-ne*y[0]*ralpha(T) + fac*vol*Gamma(NH1,nh1,nhe1,y[1])[0])
    dydx[1] = np.nan_to_num(fac*vol*Gamma2(NHe1,nh1,nhe1,y[1])[0]+nhe1*ne*GE2(T)+y[2]*ne*alpha3(T) - fac*vol*Gamma3(NHe2,nh1,nhe1,y[1])[0]-y[1]*ne*GE3(T)-y[1]*ne*alpha2(T))   # He II equation
    dydx[2] = np.nan_to_num(fac*vol*Gamma3(NHe2,nh1,nhe1,y[1])[0]+y[1]*ne*GE3(T)-y[2]*ne*alpha3(T)) 
    dydx[3] = np.nan_to_num((2.0/3)*(fh*vol*(Gamma(NH1,nh1,nhe1,y[1])[1]+Gamma2(NHe1,nh1,nhe1,y[1])[1]+Gamma3(NHe2,nh1,nhe1,y[1])[1])-Cool(T,nh,y[0],y[1],y[2],ne))/(nh+4*y[1]+4*y[2]+ne) -2*Hubble(t)*T -T/(nh+ne+0.315*nh)*(dydx[0]+dydx[1]+2*dydx[2]) )
    
    """
    # Davies 2014 arxiv eqns -- no collisional ionisation
    dydx[0] = np.nan_to_num(-ne*y[0]*ralpha(T) + fac*vol*Gamma(NH1,nh1,nhe1,y[1])[0])
    dydx[1] = np.nan_to_num( fac*vol*Gamma2(NHe1,nh1,nhe1,y[1])[0]+ y[2]*ne*alpha3(T) - fac*vol*Gamma3(NHe2,nh1,nhe1,y[1])[0] - y[1]*ne*alpha2(T))   # He II equation
    dydx[2] = np.nan_to_num(fac*vol*(Gamma3(NHe2,nh1,nhe1,y[1]))[0]-y[2]*ne*alpha3(T)) # HeIII abundance
    dydx[3] = np.nan_to_num((2.0/3)*(fh*vol*(Gamma(NH1,nh1,nhe1,y[1])[0]+Gamma2(NHe1,nh1,nhe1,y[1])[1]+Gamma3(NHe2,nh1,nhe1,y[1])[1])-Cool(T,nh,y[0],y[1],y[2],ne))/(nh+4*y[1]+4*y[2]+ne) -2*Hubble(t)*T -T/(nh+ne+0.315*nh)*(dydx[0]+dydx[1]+2*dydx[2]) )
    #dydx[1] = 
    """
    #print((fac*vol*Gamma(NH1,nh1,nhe1,y[1])))    
    return dydx 


def evol_cells(lenx,n2,nhe2,nhe3,T,NHcell,NHe1,NHe2,R,nh,ti,tsim):
	Y =0.24
	nl2  = np.zeros(lenx)
	nhel2 = np.zeros(lenx)
	nhel3 = np.zeros(lenx)
	Tl = np.zeros(lenx)
	lNHcell = np.zeros(lenx)
	lNHe1 = np.zeros(lenx)
	lNHe2 = np.zeros(lenx)
	
	# Solving rate equations in each cell at global time = tsim
	for k in range(lenx):
		z0 = np.array([n2[k],nhe2[k],nhe3[k],T[k]])
		nl2[k] = solve_ivp(model,(ti,tsim),z0,method='Radau',args=(NHcell[k],NHe1[k],NHe2[k],R[k],nh[k]),t_eval=[tsim],max_step=tsim-ti).y[0,:]
		nhel2[k] = solve_ivp(model,(ti,tsim),z0,method='Radau',args=(NHcell[k],NHe1[k],NHe2[k],R[k],nh[k]),t_eval=[tsim],max_step=tsim-ti).y[1,:]
		nhel3[k] = solve_ivp(model,(ti,tsim),z0,method='Radau',args=(NHcell[k],NHe1[k],NHe2[k],R[k],nh[k]),t_eval=[tsim],max_step=tsim-ti).y[2,:]
		Tl[k] = solve_ivp(model,(ti,tsim),z0,method='Radau',args=(NHcell[k],NHe1[k],NHe2[k],R[k],nh[k]),t_eval=[tsim],max_step=tsim-ti).y[3,:]   
	
	# Computing optical depth upto each cell for next global timestep
	for l in range(lenx):
		lNHcell[l] = np.cumsum(nh-n2)[l]*dr*3e7
		lNHe1[l] = np.cumsum((Y/(4*(1-Y)))*nh - nhe2 - nhe3)[l]*dr*3e7
		lNHe2[l] = np.cumsum(nhe2)[l]*dr*3e7
        
      
	return nl2,nhel2,nhel3,Tl,lNHcell,lNHe1,lNHe2    




# Seting up grid 
#x = np.arange(0,10000,10)
#Bolton 
#Davies 
x = np.arange(0,5000.0,10)
#x = np.arange(0,160000,78.125)
z0 = np.zeros(4)
dr = x[1] - x[0]     # in units of kpc
print(dr)
print(len(x))

#Source 
Nion    = 3.15e13*1e57#3.15e13*7.7e57#4 # in units of Myr^-1
# Bolton
#Nion    = 3.15e13*5e53

# initial conditions and defining arrays
count   = 100000		# random upper limit on no of iterations/global time steps
Y = 0.24
IF= 0
#Y = 0
nh = np.zeros(len(x))		# in units cm^-3
n2  = np.zeros((count,len(x))) 
nhe2 = np.zeros((count,len(x)))
nhe3 = np.zeros((count,len(x)))
T = np.zeros((count,len(x)))		# in units K
NHcell = np.zeros((count,len(x)))		# in units nm^-2
NHe1 = np.zeros((count,len(x)))		
NHe2 = np.zeros((count,len(x)))
R = np.zeros(len(x))	# distance from source in kpc
tc = np.zeros(count)    # to store tsim
davt = np.zeros(count)  # dt from Davies 2016
#for i in range(len(x)):
    
    #Bolton
    #T[0][i] = 1.5e4
    #nh[i] = 7e-5
    #n2[0][i] = 0.00*nh[i]
    #nhe2[0][i] = 0
    #nhe3[0][i] = 0

    # Cosmological simulation 
    #T[0][i] = np.loadtxt("./10TI.txt",float)[i]#1e4
    #nh[i] = np.loadtxt("./10nh.txt",float)[i]#115e-6
    #n2[0][i] = nh[i]-np.loadtxt("./10nhI.txt",float)[i]#0.0*nh[i]

#initial density,temperature and optical depth on grid 	
for i in range(len(x)):
    T[0][i] = 10 #1.5e4
    nh[i] = 7e-5#1e-3#7e-5
    n2[0][i] = 0.05*nh[i]
    nhe2[0][i] = 0
    nhe3[0][i] = 0 #(Y/(4*(1-Y)))*nh[i] 
    NHcell[0][i] = np.cumsum(nh-n2[0])[i]*dr*3e7
    NHe1[0][i] = np.cumsum((Y/(4*(1-Y)))*nh - nhe2[0] - nhe3[0])[i]*dr*3e7
    NHe2[0][i]=np.cumsum(nhe2[0])[i]*dr*3e7
    R[i] = (1+i)*dr + 10

# determining IF position
#for i in range(len(x)):
#	if n2[0][i]/nh[i] >0.5:
#		IF = i 
#RIF =   (1+IF)*dr +10		# IF radius in units of kpc  
#print(RIF)
#dt = (np.power(3e21,3)/3e7)*np.pi*RIF*RIF*NHcell[0][IF]/Nion  # Davies timestep in units of Myr
#print(dt)
#dt = 1e-5




#global time step (Myr) - Bolton 2007 
epsilon = 0.1#0.05
dt  = 3261.6*(epsilon/0.1)*(dr/10.0)*1e-6 
#dt = 1e-3
print(dt)

# initial and final time in Myr
ti =  935.7#752
tf = ti + 10
#dt = 1e-5
z0 = np.zeros(2)
# j is counter for global time step
for j in range(count):
	tsim = ti + dt 
	t = np.linspace(ti,tsim)
	print(tsim,dt)
	tc[j] = tsim
	
        # evolve all cells for time tsim and save the column density NHcell
	n2[j+1],nhe2[j+1],nhe3[j+1],T[j+1],NHcell[j+1],NHe1[j+1],NHe2[j+1] = evol_cells(len(x),n2[j],nhe2[j],nhe3[j],T[j],NHcell[j],NHe1[j],NHe2[j],R,nh,ti,tsim)
	#print(NHcell[j+1][l])
	
	# update time for next iteration
	ti = tsim
	
	# determining IF position
	#for i in range(len(x)):
	#	if n2[j+1][i]/nh[i] >0.5:
	#		IF = i 
	#print(NHcell[j+1][1])
	#if NHcell[j+1][1]<0:
	#	print('negative')#sys.exit()
	#print(i)
	#RIF =  (1+IF)*dr+10   # in units of kpc 
	
	# If change in total optical depth > 10%, timestep is reduced by half, else unchanged
	if np.abs((np.cumsum(-NHcell[j+1]+NHcell[j])[len(NHcell)-1])/(np.cumsum(NHcell[j])[len(NHcell)-1]))>0.1 :
		dt = dt/2
	else :
		dt = 3261.6*(epsilon/0.1)*(dr/10.0)*1e-6
		#if tsim> 935.7+0.1:		#assuming equilibrium has been reached based on NH values
		#	dt = 1e-1
	
	# time step from Davies 2016		
	#davt[j] = (np.power(3e21,3)/3e7)*np.pi*RIF*RIF*NHcell[j+1][IF]/Nion   
	#dt = davt[j]
	#print("timestep", dt)
    
	# stopping criterion
	if np.abs(tf-tsim)<=dt:
		print(n2[j+1])
		break
	#sys.exit()        

""" SAVING SOLUTION """
#np.savetxt("n2Res_test_dr_10kpc.txt",n2[j+1],delimiter=',')
#np.savetxt("TRes_test_dr_10kpc.txt",T[j+1],delimiter=',')
"""
np.savetxt("bhe_01myr.txt",y4_plot,delimiter=',')
np.savetxt("bhe2_01myr.txt",y2_plot,delimiter=',')
np.savetxt("bhe3_01myr.txt",y3_plot,delimiter=',')
np.savetxt("bT_01myr.txt",y1_plot,delimiter=',')
"""

""" PLOTTING SOLUTION"""


y_plot = np.zeros(len(x))
y1_plot = np.zeros(len(x))
y = np.zeros(len(x))
y2_plot = np.zeros(len(x))
y3_plot = np.zeros(len(x))
y4_plot = np.zeros(len(x))
Y = 0.24
f = Y/(4*(1-Y))
# calculating abundances using Defn of Davies 2016
for i in range(len(x)):
    #y2_plot[i] = (f*nh[i]-nhe2[j+1][i])/(f*nh[i])
    y2_plot[i]= nhe2[j+1][i]/(f*nh[i])
    y_plot[i] = (nh[i]-n2[j+1][i])/(nh[i])
    y3_plot[i] = nhe3[j+1][i]/(f*nh[i])
    y1_plot[i] = T[j+1][i]
    y4_plot[i] = (f*nh[i]-nhe2[j+1][i]-nhe3[j+1][i])/(f*nh[i])


"""
plt.semilogy(x,y_plot,label='xHI',color='k')
plt.semilogy(x,y4_plot,label='xHeI',color='yellow')
plt.semilogy(x,y2_plot,label='xHeII',color='b')
plt.semilogy(x,y3_plot,label='xHeIII',color='g')
#plt.semilogy(x,y1_plot,label='T')
plt.xlabel('R(kpc)')
plt.legend()
plt.show()
"""
fig, ax1 = plt.subplots()
#color = 'tab:black'
ax1.set_xlabel('R(kpc)')
ax1.set_ylabel('log(xi)' ,color='k')
ax1.semilogy(x, y_plot,label='HI',color='k')
ax1.semilogy(x,y4_plot,label='HeI',linestyle='--',color='y')
ax1.semilogy(x,y2_plot,label='HeII',linestyle='--',color='b')
ax1.semilogy(x,y3_plot,label='HeIII',linestyle='--',color='g')
ax1.legend(loc="lower right")
#ax1.set_xlim(0,650)
#ax1.set_ylim(1.1)
ax1.tick_params(axis='y', labelcolor='k')
"""
"""
ax2 = ax1.twinx()   # instantiate a second axes that shares the same x-axis

#plt.title('t= 14 Myr, cell size = 90kpc, Ti =10K')
color = 'tab:red'
ax2.set_ylabel('log(T(K))', color=color)    # we already handled the x-label with ax1
ax2.semilogy(x, y1_plot, color=color)
#ax2.set_xlim(0,900)
#ax1.ylim(min(y_plot),max(y_plot))
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout() # otherwise the right y-label is slightly clipped


"""

plt.savefig("10myr_tau_total.png",bbox_inches ='tight')
#np.savetxt("10times.txt",tc,delimiter=',')
#np.savetxt("10T@times_d10.txt",T[:,10],delimiter=',')
#np.savetxt("10T@times_d100.txt",T[:,100],delimiter=',')
#np.savetxt("timesteps.txt",davt,delimiter=',')

plt.show()






#collisional ionisation term 


#temp at fixed radius
#ionisation timescale
#change timestep



































