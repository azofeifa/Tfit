import numpy as np
import matplotlib.pyplot as plt
def runOne(mu=0, s=1, l=5, lr=100, ll=-100, we=0.5,wl=0.25, wr=0.25, pie=0.5, pil=0.1, pir=0.9, N=1000, SHOW=False , bins=200, noise=False ):
	print s
	forward 	 = list(np.random.normal(mu, s, int(N*we*pie)) + np.random.exponential(l, int(N*we*pie) ))
	forward 	+= list(np.random.uniform(mu, lr, int(N*wr*pir)))
	
	reverse 	 = list(np.random.normal(mu, s, int(N*we*(1-pie))) - np.random.exponential(l, int(N*we*(1-pie) )))
	reverse 	+= list(np.random.uniform(ll, mu, int(N*wl*(1-pil))))

	#simulate noise?
	if noise:
		forward += list(np.random.uniform(ll-50, lr+50, int(N*0.05) ))
		reverse += list(np.random.uniform(ll-50, lr+50, int(N*0.05) ))



	X 			= np.zeros((bins, 3))
	X[:,0] 		= np.linspace(min(reverse), max(forward) ,bins)
	for j,f in enumerate((forward, reverse)):
		for x in f:
			i=0
			while i < X.shape[0] and X[i,0] <= x:
				i+=1
			X[i-1, j+1]+=1
	if SHOW:
		plt.figure(figsize=(15,10))
		plt.bar(X[:,0], X[:,1], width = (X[-1,0]- X[0,0]) / bins)
		plt.bar(X[:,0], -X[:,2], color="red", width = (X[-1,0]- X[0,0]) / bins)
		plt.show()
	return X
def BIN(forward, reverse, bins, SHOW=False):
	start, stop 	= min((min(forward)[0],min(reverse)[0])),max((max(forward)[0],max(reverse)[0]))
	X 				= np.zeros((bins, 3))
	X[:,0] 			= np.linspace(start, stop, bins)
	for j,f in enumerate((forward, reverse)):
		for x,c in f:
			i  	= 0
			while i < X.shape[0] and X[i,0] <= x:
				i+=1
			X[i-1, j+1]+=c
	#want to filter out the zeros
	N 				= 0
	for i in range(X.shape[0]):
		if X[i,1] or X[i,2]:
			N+=1
	nX 				= np.zeros((N,3))
	j 				= 0
	for i in range(X.shape[0]):
		if X[i,1] or X[i,2]:
			nX[j,:] 	= [X[i,0], X[i,1], X[i,2]]
			j+=1


	if SHOW:
		plt.figure(figsize=(15,10))
		plt.bar(X[:,0], X[:,1], width = (X[-1,0]- X[0,0]) / bins)
		plt.bar(X[:,0], -X[:,2], color="red", width = (X[-1,0]- X[0,0]) / bins)
		plt.show()
	return nX



	






