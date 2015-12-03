
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate




def simulate(N=10000, mu=0, si=1, a=-5,b=5):
	f 	= lambda x: (1.0 / math.sqrt(2*math.pi*pow(si,2))    )*math.exp(-pow(x-mu,2)/(2*pow(si,2)))
	u 	= lambda x: 1.0 / (b-a)
	xs 	= list()
	for i in range(N):
		x 		= np.random.uniform(a,b)
		F,U  	= f(x), u(x)
		r 	 	= f(x) / (u(x) + f(x) )
		if np.random.uniform(0,1) > r:
			xs.append(x)
	c 	= 0.06
	g 	= lambda x: c*((u(x) / (f(x)+u(x))))
	F 	= plt.figure()
	ax 	= F.add_subplot(1,2,1)
	ax.hist(xs,bins=50, normed=1)
	ax2 	= F.add_subplot(1,2,2)
	xs.sort()
	N 		= float(len(xs))
	ax2.plot(xs,[i/N  for i in range(len(xs))])
	
	plt.show()
	return xs


if __name__ == "__main__":
	X 	= simulate()
