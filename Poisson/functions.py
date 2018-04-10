from math import *
import matplotlib.pyplot as plt
import numpy as np

def Poisson_Ratio(x,m1,m2):
	n=max(20,m1+m2)
	m1=float(m1)
	m2=float(m2)
	out=exp(-m1-m2)
	for i in range(1,n+1):
		r0=float(i)
		out=out+exp(r0*x*(1+log(m1/(r0*x)))-m1+r0*(1+log(m2/r0))-m2)/r0 
	return out/(2.*pi*sqrt(x))
	
if __name__ == "__main__":
	x=np.arange(0.1,6.1,0.01)
	y1=np.zeros(len(x))
	y2=np.zeros(len(x))
	y3=np.zeros(len(x))
	for i in np.arange(len(x)):
		y1[i]=Poisson_Ratio(x[i],10,10)
		y2[i]=Poisson_Ratio(x[i],20,10)
		y3[i]=Poisson_Ratio(x[i],10,20)
	plt.figure()
	plt.plot(x,y1)
	plt.plot(x,y2)
	plt.plot(x,y3)
	plt.show()
	quit()
	
