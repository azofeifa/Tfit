
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

def get_w(X, si=1, mu=0) :
	w 	= 0.01
	converged 	= False
	T 	= 1000
	t 	= 0
	a,b = min(X[:,0]),max(X[:,0])
	f 	= lambda x: (1.0 / math.sqrt(2*math.pi*pow(si,2))    )*math.exp(-pow(x-mu,2)/(2*pow(si,2)))
	u 	= lambda x: 1.0 / (b-a)
	g 	= lambda x: ((u(x) / (f(x)+u(x))))
	c 	= integrate.quad(g, a, b)[0] #find normalizing factor via numerial integrator
	p 	= lambda x: (1/c)*g(x)
	N 	= X.shape[0]
	prevll 	= np.inf
	while t < T and not converged:
		NULL 	= 0
		DEPLETE = 0
		EX 	= 0.
		ll 	= 0.
		for i in range(N):
			n,d 	= u(X[i,0])*(1-w),p(X[i,0])*w
			ll+=math.log(n/(n+d))*X[i,1] + math.log(d/(n+d))*X[i,1]
			NULL+=n/(n+d)*X[i,1]

			DEPLETE+=d/(n+d)*X[i,1]
		
		w 	= DEPLETE / (NULL + DEPLETE)
		if abs(prevll - ll) < 0.01:
			converged=True 
		prevll 	= ll
		t+=1
	return w



def simulate(N=10000, mu=0, si=1, a=-5,b=5, w=0.2, bins=200):
	f 	= lambda x: (1.0 / math.sqrt(2*math.pi*pow(si,2))    )*math.exp(-pow(x-mu,2)/(2*pow(si,2)))
	u 	= lambda x: 1.0 / (b-a)
	xs 	= list()
	for i in range(N):
		x 		= np.random.uniform(a,b)
		F,U  	= f(x), u(x)
		r 	 	= f(x) / (u(x) + f(x) )
		if np.random.uniform(0,1) < w:
			if np.random.uniform(0,1) > r:
				xs.append(x)
		else:
			xs.append(np.random.uniform(a,b))
	
	F 	= plt.figure()
	ax 	= F.add_subplot(1,1,1)
	counts,edges 	= np.histogram(xs, bins=bins)
	X 	= np.zeros((bins,2))
	X[:,0] 	= edges[1:]
	X[:,1] 	= counts
	W 		= get_w(X)
	ax.hist(xs,bins=bins, normed=1, alpha=0.5)
	xs 	= np.linspace(a,b,100)
	g 	= lambda x: ((u(x) / (f(x)+u(x))))
	c 	= integrate.quad(g, a, b)[0]

	p 	= lambda x: (1/c)*g(x)
	cov = lambda x: W*p(x) + (1-W)*u(x)
	ax.set_title("Depletion to Noise (EM): " + str(W) + ", actual: " + str(w)  )
	ax.plot(xs, map(cov, xs), linewidth=4.0, linestyle="--")
	plt.show()
	return xs




if __name__ == "__main__":
	mu,si 	= 0,100
	a,b 	= -1500,1500
	f 	= lambda x: (1.0 / math.sqrt(2*math.pi*pow(si,2))    )*math.exp(-pow(x-mu,2)/(2*pow(si,2)))
	u 	= lambda x: 1.0 / (b-a)
	g 	= lambda x: ((u(x) / (f(x)+u(x))))
	c 	= integrate.quad(g, a, b)[0] #find normalizing factor via numerial integrator
	print c
