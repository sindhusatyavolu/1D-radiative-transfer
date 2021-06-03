from __future__ import print_function
import numpy as np  


def freqbin(b,a,nbins,nbinint):
	for i in range(nbinint):
		print(a*np.power(b/a,i/(nbins-1)),end=",")
	return 0

freqbin(13.6*10,13.6,20.0,20)	