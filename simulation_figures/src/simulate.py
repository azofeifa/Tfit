import numpy as np
import matplotlib.pyplot as plt
import math as m
def runOne(mu=0, s=1, l=5, lr=100, ll=-100, we=0.5,wl=0.25, wr=0.25, pie=0.5, pil=0.1, pir=0.9, N=1000, SHOW=False , bins=200, noise=False, foot_print = 0 ):

	forward 	 = list(np.random.normal(mu+foot_print, s, int(N*we*pie)) + np.random.exponential(l, int(N*we*pie) ))
	forward 	+= list(np.random.uniform(mu+foot_print, lr, int(N*wr)))
	
	reverse 	 = list(np.random.normal(mu-foot_print, s, int(N*we*(1-pie))) - np.random.exponential(l, int(N*we*(1-pie) )))
	reverse 	+= list(np.random.uniform(ll, mu-foot_print, int(N*wl)))

	#simulate noise?
	if noise:
		forward += list(np.random.uniform(ll-130, lr+130, int(N*0.05) ))
		reverse += list(np.random.uniform(ll-130, lr+130, int(N*0.05) ))



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

def runMany(K=1, RGE=(0,500), N=1000, SHOW=False, bins=500):
	delta 	= float(RGE[1]-RGE[0])/float(K)
	forward, reverse 	= list(), list()
	for k in range(K):
		mu 			 = delta*k
		forward 	+= list(np.random.normal(mu ,3, int(N*0.5)) + np.random.exponential(10, int(N*0.5)))
		reverse 	+= list(np.random.normal(mu , 3, int(N*0.5)) - np.random.exponential(10, int(N*0.5)))

	forward += list(np.random.uniform(RGE[0], RGE[1], int(N*0.15) ))
	reverse += list(np.random.uniform(RGE[0], RGE[1], int(N*0.15) ))
	X 			= np.zeros((bins, 3))
	X[:,0] 		= np.linspace(min(reverse), max(forward) ,bins)
	for j,f in enumerate((forward, reverse)):
		for x in f:
			i=0
			while i < X.shape[0] and X[i,0] <= x:
				i+=1
			X[i-1, j+1]+=1
	X[:,0]-=X[0,0]	
	if SHOW:
		plt.figure(figsize=(15,10))
		X[:,0]*=100

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
def weighted_mean(x,y):
	N 	= sum(y)
	return sum([x[i]*y[i] for i in range(len(x))])/float(N)
def weird_variance(x,y,mu,direct):
	S 	= 0
	N 	= 0
	for i in range(len(x)):
		if direct==1 and x[i] < mu:
			S+=pow(x[i]-mu,2)*y[i]
			N+=y[i]
			pass
		elif direct==-1 and x[i] > mu:
			S+=pow(x[i]-mu,2)*y[i]
			N+=y[i]
	return S/N

if __name__=="__main__":
	X 	= runOne(mu=0, s=20, l=5, lr=100, ll=-100, we=0.5,wl=0.25, wr=0.25, pie=0.5, pil=0.1, pir=0.9, N=10000, SHOW=False , bins=200, noise=False )
	mu 	= 0.5*(weighted_mean(X[:,0], X[:,1])+weighted_mean(X[:,0], X[:,2]))
	l 	= 0.5*(weighted_mean(X[:,0], X[:,1])-weighted_mean(X[:,0], X[:,2]))
	print mu,l
	print m.sqrt(weird_variance(X[:,0], X[:,1], mu, 1))
	print m.sqrt(weird_variance(X[:,0], X[:,2], mu, -1))
	
	bins 	= len(X)
	plt.bar(X[:,0], X[:,1],width = (X[-1,0]- X[0,0]) / bins,alpha=0.2)
	plt.bar(X[:,0], -X[:,2], color="red",width = (X[-1,0]- X[0,0]) / bins,alpha=0.2)

	plt.scatter([mu], [0], s=20)
	plt.show()



	






