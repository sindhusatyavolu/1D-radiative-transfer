import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.ticker

#x = np.loadtxt("tminxh1_01myr.txt",float)[:,0]
#xh1 = np.loadtxt("tminxh1_01myr.txt",float)[:,1]
#T = np.loadtxt("tminT_01myr.txt",float)

x = np.loadtxt("t2dxh1_1myr.txt",float)[:,0]
xh1 = np.loadtxt("t2dxh1_1myr.txt",float)[:,1]
xhe1 = np.loadtxt("t2dxhe1_1myr.txt",float)
xhe2 = np.loadtxt("t2dxhe2_1myr.txt",float)
xhe3 = np.loadtxt("t2dxhe3_1myr.txt",float)
T = np.loadtxt("t2dT_1myr.txt",float)
dx = np.loadtxt("davies_xh1.csv",float,delimiter=',')[:,0]
dxh1 = np.loadtxt("davies_xh1.csv",float,delimiter=',')[:,1]
dx2 = np.loadtxt("davies_he1.csv",float,delimiter=',')[:,0]
dxhe1 = np.loadtxt("davies_he1.csv",float,delimiter=',')[:,1]
dx3 = np.loadtxt("davies_he2.csv",float,delimiter=',')[:,0]
dxhe2 = np.loadtxt("davies_he2.csv",float,delimiter=',')[:,1]
dx4 = np.loadtxt("davies_he3.csv",float,delimiter=',')[:,0]
dxhe3 = np.loadtxt("davies_he3.csv",float,delimiter=',')[:,1]
dx5 = np.loadtxt("davies_T.csv",float,delimiter=',')[:,0]
dT = np.loadtxt("davies_T.csv",float,delimiter=',')[:,1]

#t = np.loadtxt("Tat5thcell",float)[:,0]
#Tt = np.loadtxt("Tat5thcell",float)[:,1]
#T = np.array([3.037909e+04,3.360990e+04,5.126823e+05,6.939858e+04,1.527144e+06,5.964988e+05,7.674479e+05,2.183442e+04,5.290894e+04,1.264389e+05,3.557291e+05,2.494593e+04,1.206141e+06,1.386834e+05,1.964146e+04,1.152309e+04,1.000000e-02,3.675239e+04,1.000000e-02,1.549435e+04,1.687758e+04,2.741256e+04,1.000000e-02,1.572319e+04,1.458598e+04,1.299459e+04,1.098889e+04,1.006253e+04,8.511984e+03,8.812355e+03,7.742132e+03,7.342575e+03,6.777501e+03,6.790780e+03,6.372566e+03,5.679961e+03,2.124439e+03,2.221086e+03,2.291388e+03,2.342590e+03,2.379088e+03,2.404975e+03,2.423892e+03,2.434151e+03,3.213162e+03,6.844212e+02,7.181170e+02,7.468682e+02,7.711705e+02,7.915196e+02,9.816877e+01,9.830859e+01,9.843579e+01,9.855167e+01,9.865742e+01,9.875406e+01,9.884249e+01,9.892353e+01,9.899789e+01,9.906622e+01,9.912906e+01,9.918695e+01,9.924032e+01,9.928959e+01,9.933511e+01,9.937723e+01,9.941622e+01,9.945237e+01,9.948591e+01,9.951705e+01,9.954600e+01,9.957292e+01,9.959800e+01,9.962136e+01,9.964315e+01,9.966349e+01,9.968248e+01,9.970023e+01,9.971684e+01,9.973239e+01,9.974695e+01,9.976059e+01,9.977339e+01,9.978540e+01,9.979668e+01,9.980728e+01,9.981725e+01,9.982662e+01,9.983545e+01,9.984376e+01,9.985159e+01,9.985897e+01,9.986593e+01,9.987250e+01,9.987870e+01,9.988456e+01,9.989010e+01,9.989533e+01,9.990028e+01,9.990497e+01])
fig,ax = plt.subplots()

#print(T)

#ax.loglog(x, T/(xh1*7e-5), color="black")
#ax.set_xlabel("t(Myr)")
#ax.set_xlabel("R(kpc)",fontsize=14)
#ax.set_ylabel("$\Gamma_{HI}$")

#ax.set_ylabel("log(T(K))")
#plt.title("Temperature profile at a distance 50kpc")
#logfmt = matplotlib.ticker.LogFormatterExponent(base=10.0, labelOnlyBase=True)
ax.semilogy(x, xh1,label='xHI',color="black")
#ax.semilogy(dx, dxh1,label='dxHI',color="black",linestyle='--')
#ax[0].yaxis.set_major_formatter(logfmt)

ax.semilogy(x,xhe1,label='HeI',color = "orange")
#ax.semilogy(dx2,dxhe1,label='HeI',color = "orange",linestyle='--')
ax.semilogy(x,xhe2,label='HeII',color="blue")
#ax.semilogy(dx3,dxhe2,label='HeII',color="blue",linestyle='--')
#ax[0].set_ylabel("log($x_{HeII}$)",color="black",fontsize=14)
#ax[1].set_xlim(0,1000)
ax.semilogy(x,xhe3,label='He III',color="green")
#ax.semilogy(dx4,dxhe3,label='He III',color="green",linestyle='--')
# set x-axis label
ax.set_xlabel("R(kpc)",fontsize=14)
# set y-axis label
#ax[0].set_ylabel("log($x_{HI}$)",color="black",fontsize=14)
#ax.set_xlim(0,500)
# twin object for two different y-axis on the sample plot
plt.legend(loc="lower right")
ax2=ax.twinx()
# make a plot with different y-axis using second axis object
ax2.semilogy(x,T,color="red")
#ax2.semilogy(dx5,dT,color="red",linestyle='--')
ax2.set_ylabel("T(K)",color="red",fontsize=14)
#ax2.set_ylim(np.min(T),np.max(T))
#ax2.set_xlim(0,150)
#plt.title("alpha corrected, 0.01Myr, run time= 1189s $ Photon rate = 5*10^{53} s^{-1}$")
#plt.title("Mnochromatic source, 1Myr,  $ Photon rate = 5*10^{53} s^{-1}$, resolution = 10kpc") #took 4800 secs to run
#plt.title("1Myr,  $ Photon rate = 5*10^{53} s^{-1}$, resolution = 10kpc, only recombination cooling")
#plt.title(" Explicit scheme for T,Photon rate $= 5*10^{53} s^{-1}$,$T_{i}= 10 K$, subcycling x,T")
plt.title("Davies 2016 parameters including timestep, $t_{Q}=10Myr$ dR 10kpc")
plt.legend(loc="lower right")

plt.show()
# save the plot as a file
#fig.savefig('two_different_y_axis_for_single_python_plot_with_twinx.jpg',
#            format='jpeg',
#            dpi=100,
#            bbox_inches='tight')
