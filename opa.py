#opacity interpolation
#using opacity table from OPAL - table 66 w/ solar characteristics and X=0.7 Y=0.3
#import modules
def opa(logT,logrho):
	import numpy as np
	import math
	from scipy import interpolate

#first ask desired log T and log rho
#remember log rho must be converted to log R (R is rho/(T/1e6)**3)

	rho=10.**logrho
	T=10.**logT
	R=rho/((T*1.0e-6)**3.)
	logR=math.log10(R)
	#print R,T,logT,logR
	
#import the table
	tab=np.genfromtxt('opac_table.txt',skiprows=2,filling_values=1)
	Rtab = tab[0,1:]
	Ttab = tab[1:,0]
	ktab = tab[1:,1:]
	f = interpolate.interp2d(Rtab, Ttab, ktab)
	znew = f(logR, logT)
	#print znew
	k=10.**znew
	return k
	#print 'Kappa is: ', k
	#print 'Log_kappa is:', logk