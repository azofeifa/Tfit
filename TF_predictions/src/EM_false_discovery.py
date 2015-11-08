import numpy as np
import matplotlib.pyplot as plt
import math
def simulate(w=0.5, s=10, N=1000, SHOW=True):
	XS 	= list()
	for n in range(N):
		U 	= np.random.uniform(0,1)
		if U < w:
			x 	= np.random.normal(0,s,1)[0]
		else:
			x 	= np.random.uniform(-100,100)
		XS.append(x)
	if SHOW:
		plt.hist(XS, bins=200)
		plt.show()
	counts,edges 	= np.histogram(XS, bins=200)
	edges 			= edges[1:]
	X 				= np.zeros((len(counts), 2))
	X[:,0] 			= edges
	X[:,1] 			= counts
	return X
def normal(x, mu, si):
	return (1.0 / (math.sqrt(2*math.pi) * si ))* math.exp(-pow(x-mu,2)/(2*pow(si,2)))
def uniform(x,a,b):
	return 1.0 / (b-a)

import math
def get_w(d, mu=0, s=100, bins=100):
	counts,edges = np.histogram(d, bins=bins)
	X = np.zeros((bins, 2))
	edges=(edges[1:]+ edges[:-1])/2.

	X[:,0]=edges
	X[:,0]*=1000
	X[:,1]=counts
	def normal(x, mu, si):
		return (1.0 / (math.sqrt(2*math.pi) * si ))* math.exp(-pow(x-mu,2)/(2*pow(si,2)))
	def uniform(x,a,b):
		return 1.0 / (b-a)
	w =	0.5
	t = 0
	a,b = min(X[:,0]), max(X[:,0])
	max_iterations=100
	while t < max_iterations:
		EN, EU = 0,0
		ll 	= 0
		for i in range(X.shape[0]):

			x,y = X[i,:]
			rn= w*normal(x, mu, s) / (w*normal(x, mu, s) + (1-w)*uniform(x, a,b))
			ru = (1-w)*uniform(x, a,b) / (w*normal(x, mu, s) + (1-w)*uniform(x, a,b))
			EN+=(rn*y)
			EU+=(ru*y)
			ll+=math.log(rn + ru)*y
		w = EN / (EU+ EN)
		t+=1
	return w
def load_FILE(FILE):
	G = {}
	with open(FILE) as FH:
		header = True
		for line in FH:
			if not header:
				TF, uni_pval, cent_pval, bi, d = line.split("\t")[:5]

				d = [float(x) for x in d.strip("\n").split(",") if x ]
				w=get_w(d,bins=200)
				if  w > 0.26:
					plt.hist(d,bins=100)
					plt.show()
				if w < 0.1:
					print TF, "BAD"
					plt.hist(d,bins=100)
					plt.show()
					
				G[TF] = (uni_pval, cent_pval, bi,w)
			else:
				header=False
	return G



if __name__ == "__main__":
	X 	= simulate(w=0.3, N=10000, SHOW=False)
	F 	= load_FILE("/Users/joazofeifa/Lab/ENCODE/HCT116/Allen2014_DMSO2_3-1_bidirectional_hits_intervals.txt")
	print fit(X)

