#Boundary conditions for surface of star for input to derivs
#Guess M, L, and R

	 #def load2():
def load2(rg,lg):
	import numpy as np
	import math
	import opa as opa
	from scipy.optimize import fmin #minimizer for initial guess
	#minimizing function (P_1 - P_2)
	def func(x):
		G=6.674e-8 #cm3g-1s-2
		R=8.315e7
		sig=5.67e-5
		mass=4e33
		X=0.7
		Y=0.28
		#print (abs(x[2])
		mu = 2./(1+(3.*X)+(Y/2.))
		if (4.*math.pi*(abs(x[0])**2.)*sig) <= 0 or abs(x[2]) < 1e30:
			q = 99
		else:
			T_f=(abs(x[2])/(4.*math.pi*(abs(x[0])**2.)*sig))**(1./4.)
			#print T_f,abs(x[2])
			logrho=math.log10(abs(x[1]))
			logT=math.log10(T_f)
			k=opa.opa(logT,logrho)
			q= abs(((G*mass)/abs(x[0])**2.)*(2./3.)*(1./k) - (R/mu)*abs(x[1])*T_f)/(((G*mass)/abs(x[0])**2.)*(2./3.)*(1./k))
		return q

	#initial guesses for R, rho, L
	#x0=[9e10,1e-7,6.4e33]
	x0=[rg,1e-7,lg]
	res=fmin(func,x0,maxiter=100,disp=0)
	G=6.674e-8 #cm3g-1s-2
	R=8.315e7
	sig=5.67e-5
	mass=4e33
	rad=abs(res[0])
	X=0.7
	Y=0.28
	mu = 2./(1+(3.*X)+(Y/2.))
	rhos=abs(res[1])
	logrho=math.log10(rhos)
	L=abs(res[2])
	T_f=(L/(4.*math.pi*(rad**2.)*sig))**(1./4.)
	logT=math.log10(T_f)
	k=opa.opa(logT,logrho)
	P_1=((G*mass)/rad**2.)*(2./3.)*(1./k)
	P_2=(R/mu)*rhos*T_f

	P_f=(P_1+P_2)/2.
	#print "rf: ",rad
	#print "lf: ",L
	#print "Pf: ",P_f
	#print "Tf: ",T_f
	return rad,L,P_f,T_f,rhos

#best so far:rad=8.67e10, L=6.33e37, rhos=3.43e-8, logP=18.1 logT =8
#higher surface r, lower surface T, higher surface L, lower surface P
#rad=8.67e11, rhos=1.31e-9,L=1.33e38,logP=18.1 logT =8
#still higher r, now higher T, lower L, higher P

#higher r, higher t, lower l, higher p
#w/ 18, 8: lower r, highr t, higher l, higher p