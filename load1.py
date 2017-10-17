#Boundary conditions for center of star for input to derivs
#Guess T_core and P_core
def load1(Pg,Tg):
	 #def load1():
	import numpy as np
	import math
	import opa as opa #opacity interpolation
	import rho
	#rho_core calculated from logP, logT, X, Y
	import en #nuclear energy generation
	G=6.674e-8#cm3g-1s-2
	mtot=4e33
	dm=(1e-6*mtot)
	a = 7.56e-15#ergs cm-3 K-4
	c= 3e10#cm/s
	P = Pg
	T = Tg
	#P = input("Please enter the log P_c: ")
	#T = input("Please enter the log T_c: ")
	X=0.7
	#Chemical composition from OPAL opacity table used
	Y=0.28
	rho_c=rho.rho(P,T,X,Y)
	if rho_c <=0:
		rho_c = 74
	inT=10**T
	#Little step based on dm
	r=(3./(4.*math.pi*rho_c))**(1./3.)*dm**(1./3.)
	l=en.en(X,Y,rho_c,inT)*dm
	P_1=-((3*G)/(8*math.pi))*((4*math.pi*rho_c)/3.)**(4./3.)*dm**(2./3.)+(10.**P)
	log10rho=math.log10(rho_c)
	delrad=(3*opa.opa(T,log10rho)*l*(10.**P))/(16.*a*math.pi*c*G*dm*(10.**T)**4.)
	#Determine if T_conv or T_rad is more appropriate
	if delrad >= 0.4:
		lnT=-(math.pi/6.)**(1./3.)*G*((0.4*rho_c**(4./3.))/(10.**P))*dm**(2./3.)+math.log((10.**T))
		T_1=math.exp(lnT)
	else:
		T4=(-1./(2.*a*c))*(3./(4.*math.pi))**(2./3.)*opa.opa(T,log10rho)*en.en(X,Y,rho_c,inT)*rho_c**(4./3.)*dm**(2./3.) + (10.**T)**4.
		T_1=T4**(0.25)
	#print "r1: ",r
	#print "l1: ",l
	#print "P1: ",P_1
	#print "T1: ",T_1
	return r,l,P_1,T_1,rho_c