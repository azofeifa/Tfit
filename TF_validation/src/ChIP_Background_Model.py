import numpy as np
import random
import matplotlib.pyplot as plt
import math
import time
def choose(n,k):
	if k > n:
		return 0
	r=1
	N=0
	for d in range(1,int(k+1)):
		r*=(n-N)
		r/=d
		N+=1
	return r


def prob(k,w,l,M):
	p=float(2*w) / float(l- 2*w)
	vl 	= choose(M,k)*pow(p,k)*pow(1.0-p,M-k  )
	return vl
def normal(x,mu, si):
	return (1.0 / (math.sqrt(2*math.pi)*si  ))*math.exp(-pow(x-mu,2)/(2*pow(si,2)))


def sim_background(N=1000, l=1000,SN=1, w=5, SIM=5000 ):
	count=list()	
	for t in range (SIM):
		D 	= np.random.randint(0,l, N)
		counts,edges 	= np.histogram(D, bins=l,range=(0,l))
		edges 			= edges[1:]
		for i in range(SN):
			x 	= np.random.randint(w ,l-w)
			start, stop 	= x-w, x+w
			DD 	= 0
			for ct, edge in zip(counts, edges):
				if edge > stop:
					break
				if start<=edge<=stop:
					DD+=ct
			count.append(DD)
	counts,edges 	= np.histogram(count, bins=max(count)-min(count) )
	counts 			= [float(ct) / sum(counts) for ct in counts]
	edges 			= edges[:-1]
	F 				= plt.figure(figsize=(15,10))
	ax 				= F.add_subplot(1,1,1)
	ax.set_xlim(0,max(count))
	ks 	= range(min(count), max(count))
	ax.plot(edges, [prob(int(k), w,l, N) for k in edges], color="green")
	ax.bar(edges, counts, width=1,alpha=0.5)
	xs 	= np.linspace(min(count), max(count), 1000)
	mu 	= np.mean(count)
	si 	= np.std(count)
	ax.plot(xs, [normal(x, mu, si ) for x in xs])
	plt.show()
if __name__== "__main__":
	sim_background()

