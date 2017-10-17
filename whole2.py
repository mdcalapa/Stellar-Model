#derivs: takes in initial guess from core and end surface
#townsend for 1 solar mass logP,logT, R, L 1.73285653E+01,7.17183363E+00,60154047125.9,6.7222794e33
#Solar mass		Mo	1.99	33	g
#Solar radius		Ro	6.96	10	cm
#Solar luminosity	Lo	3.9	33	erg s-1
#townsend columns: log10_L=4, log10_R=5, log10_Tc=7, log10_Pc=9
#import sys
#sys.path.remove('/usr/local/scisoft/packages/python/lib/python2.7/site-packages')

def whole2():
	import numpy as np
	import math
	import opa as opa
	import load1
	import load2
	import f_derivs2
	from scipy.integrate import odeint
	import en
	import rho
	import call
	import copy
	

	x0 = input("Please enter the initial radius, luminosity, logP, logT: ")
	x0=np.array(x0)
	delRo,delLo,delPo,delTo=call.call(x0)
	diff=np.array([delRo,delLo,delPo,delTo])
	w=0
	while np.max(abs(diff)) >= 0.1:
		w=1+w
		print w
		J = np.ndarray(shape=(4,4),dtype=float)
		for i in range(0,4,1):
			h=x0[i]*0.001
			tmp=x0[i]+h
			xnew=copy.copy(x0)
			xnew[i]=tmp
			#print "xnew: ", xnew
			delRn,delTn,delLn,delPn=call.call(xnew)
			pdr = (delRo-delRn)/h
			#print pdr
			pdl = (delLo-delLn)/h
			pdp = (delPo-delPn)/h
			pdt = (delTo-delTn)/h
			J[0,i]=pdr
			J[1,i]=pdl
			J[2,i]=pdp
			J[3,i]=pdt
		#print J
		Jinv = np.linalg.inv(J)
		dx = np.dot(-1.0*Jinv,diff)
		print "dx:",dx
		dx=0.001*dx
		x0=x0+dx
		print "guess:", x0
		x0=np.array(x0)
		delRo,delLo,delPo,delTo=call.call(x0)
		diff=np.array([delRo,delLo,delPo,delTo])
	f_derivs2.f_derivs2(x0)