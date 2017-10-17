#calc linear density given logP log T X Y Z
def rho(P,T,X,Y):
	import numpy as np
	import math
#define constants
	a = 7.56e-15 #ergs cm-3 K-4
	R = 8.315e7 #erg K-1 g-1
	P = 10.**P
	T = 10.**T
	mu = 2./(1.+(3.*X)+(Y/2.))
	rho = (mu/(R*T))*(P - ((a*(T**4.))/3.))
	Pg = (R/mu)*rho*T
	Pt = Pg + ((a*T**4.)/3.)
	#beta = Pg/Pt
	return rho