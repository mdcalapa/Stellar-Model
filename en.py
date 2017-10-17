def en(inx,iny,inrho,inT):
	import math
	#pp-chain
	T6=inT/(1.0e6)
	T9=T6/1000.
	rho=inrho
	X=inx
	Y=iny
	#estimate psi based on temperature graph in Kippenhahn
	if T6/10. < 2 and T6/10. >=1:
		psi=1+(T6/10.)*0.25
	elif T6/10. < 1:
		psi = 1
	elif T6/10. >= 2 and T6/10. < 3:
		psi=1.6 - 0.72*(T6/10.)
	else:
		psi=1.6

	g=1+(3.82*(T9))+(1.51*(T9**2.))+(0.144*(T9**3.))-(0.0114*(T9**4.))
	epsp=(2.57e4)*psi*g*rho*(X**2.)*(T9**(-2./3.))*math.exp(-3.381/((T9**(1./3.))))

	#CNO
	Xcno=0.7*(1.0-(X+Y))
	g2=(1-(2*T9)+(3.41*(T9**2.))-(2.43*(T9**3.)))
	epsc=(8.24e25)*g2*Xcno*X*rho*(T9**(-2./3.))*math.exp(-15.231*(T9**(-1./3.))-(T9/0.8)**2.)
	return epsp+epsc