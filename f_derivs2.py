#derivs: takes in initial guess from core and end surface
#sys.path.remove('/usr/local/scisoft/packages/python/lib/python2.7/site-packages')
def f_derivs2(x0):
	import numpy as np
	import math
	import opa as opa
	import load1
	import load2
	from scipy.integrate import odeint
	import en
	import rho

	sig=5.67e-5
	G=6.674e-8 #cm3g-1s-2
	c= 3e10 #cm/s
	a = 7.56e-15 #ergs cm-3 K-4
	X=0.7
	Y=0.28
	mtot=4e33 #total mass, small step in for surface
	r,l,P,T,rho_c=load1.load1(x0[2],x0[3])
	y0=[r,l,P,T]
	rf,lf,Pf,Tf,rhof=load2.load2(x0[0],x0[1])
	y02=[rf,lf,Pf,Tf]
	m=1e-6
	mf=mtot-(m*mtot)

	def derivatives(y,m):
		X=0.7
		Y=0.28
		sig=5.67e-5
		G=6.674e-8 #cm3g-1s-2
		c= 3e10 #cm/s
		a = 7.56e-15 #ergs cm-3 K-4
		r=y[0]
		l=y[1]
		P=y[2]
		T=y[3]
		delrad = (3.*opa.opa(np.log10(T),np.log10(rho.rho(math.log10(P),math.log10(T),X,Y)))*l*(P))/((16.*a*math.pi*c*G*(m))*(T)**4.)
		#print "del:",delrad
		if delrad >= 0.4:
			dela=0.4
		else:
			dela=delrad
		#print "del_final:", dela
		#print "Opacity:", opa.opa(np.log10(T),np.log10(rho.rho(math.log10(P),math.log10(T),X,Y)))

		drdm=1./(4.*math.pi*(r**2.)*rho.rho(math.log10(P),math.log10(T),X,Y))
		dpdm= -(G*(m))/(4.*math.pi*(r**4.))
		dLdm= en.en(X,Y,rho.rho(math.log10(P),math.log10(T),X,Y),T)
		dTdm=-(G*(m)*T*dela)/(4.*math.pi*(r**4.)*P)
		#print "drdm:", drdm
		#print "dpdm:", dpdm
		#print "dLdm:",  dLdm
		#print "dTdm:",dTdm
		return [drdm,dLdm,dpdm,dTdm]

#iterate from center to 0.7M
	dm=1e-6
	msteps=10.**(np.linspace(np.log10(dm*mtot),np.log10(0.7*mtot),500))
	result=odeint(derivatives,y0,msteps,mxstep=20000)
	result_r,result_L,result_p,result_T=result.T
	#import matplotlib.pyplot as plt
	#plt.plot(msteps,result_r,label='Outward Integration')
	#plt.show()
	

#iterate from M to 0.3M
	mstepsin1=np.arange(1,0.99,-1e-7)*mtot
	mstepsin2=10.**(np.arange(np.log10(0.99*mtot),np.log10((0.7)*mtot),-1e-3))
	mstepsin=np.append(mstepsin1,mstepsin2)
	result2=odeint(derivatives,y02,mstepsin,mxstep=20000)
	result_r2,result_L2,result_p2,result_T2=result2.T

	print "initial:",y0
	print "final:",y02
	print "r:",result_r[499],result_r2[100151]
	print "L:",result_L[499],result_L2[100151]
	print "P:",result_p[499],result_p2[100151]
	print "T:",result_T[499],result_T2[100151]

	import matplotlib.pyplot as plt
	plt.plot(msteps,result_r,label='Outward Integration')
	plt.plot(mstepsin,result_r2,label='Inward Integration')
	plt.xlabel('Mass (g)')
	plt.ylabel('Radius (cm)')
	plt.legend(loc=2)
	plt.title('Radius')
	plt.show()
#plt.subplot(2,2,2)
	plt.plot(msteps,result_T,label='Outward Integration')
	plt.plot(mstepsin,result_T2,label='Inward Integration')
	plt.xlabel('Mass (g)')
	plt.ylabel('Temperature (K)')
	plt.legend(loc=2)
	plt.title('Temperature')
	plt.show()
#plt.subplot(2,2,3)
	plt.plot(msteps,result_L,label='Outward Integration')
	plt.plot(mstepsin,result_L2,label='Inward Integration')
	plt.xlabel('Mass (g)')
	plt.ylabel('Luminosity (erg/s)')
	plt.legend(loc=2)
	plt.title('Luminosity')	
	plt.show()
#plt.subplot(2,2,4)
	plt.plot(msteps,result_p,label='Outward Integration')
	plt.plot(mstepsin,result_p2,label='Inward Integration')
	plt.xlabel('Mass (g)')
	plt.ylabel('Pressure (bar)')
	plt.legend(loc=2)
	plt.title('Pressure')
	plt.show()