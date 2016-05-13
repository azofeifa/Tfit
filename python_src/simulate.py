import numpy as np
import matplotlib.pyplot as plt
import math as m
from matplotlib import rc
import model
font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 22}

rc('font', **font)
def runOne(mu=0, s=1, l=5, lr=100, ll=-100, we=0.5,wl=0.25, wr=0.25, pie=0.5, pil=0.1, pir=0.9, N=1000, SHOW=False , bins=200, noise=False, foot_print = 0 ):

	forward 	 = list(np.random.normal(mu+foot_print, s, int(N*we*pie)) + np.random.exponential(l, int(N*we*pie) ))
	forward 	+= list(np.random.uniform(mu+foot_print, lr, int(N*wr*pir)))
	
	reverse 	 = list(np.random.normal(mu-foot_print, s, int(N*we*(1-pie))) - np.random.exponential(l, int(N*we*(1-pie) )))
	reverse 	+= list(np.random.uniform(ll, mu-foot_print, int(N*wl*(1-pil))))

	#simulate noise?
	if noise:
		forward += list(np.random.uniform(ll-150, lr+150, int(N*0.05) ))
		reverse += list(np.random.uniform(ll-150, lr+150, int(N*0.05) ))



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
def density(xs, wb=0.5, we=0.5):
	B 	=	model.component_bidir(0.0, 2.0, 1.0/5.0, wb ,0.5, None,foot_print=0)
	U 	= 	model.component_elongation(0,100, we, 1.0, None, None, None, 0)
	ysf = [B.pdf(x,1) + U.pdf(x,1)for x in xs]
	ysr = [-(B.pdf(x,-1) + U.pdf(x,-1))for x in xs]
	return ysf, ysr


if __name__=="__main__":
	wb=0.5
	we=0.5
	X 	= runOne(mu=0, s=2, l=5, lr=100, ll=-100, we=wb,wl=0.0, wr=we, pie=0.5, pil=0.1, pir=0.9, N=1000, SHOW=False , bins=90, noise=False )
	bins 	= 100.0

	xs 	= np.linspace(X[0,0], X[-1,0]+1.0,1000)
	ysf, ysr 	= density(xs, wb=wb, we=we)
	F 		= plt.figure()
	ax 		= F.add_subplot(111)
	ax.bar(X[:,0], 0.75*X[:,1]/np.sum(X[:,1:]),width = (X[-1,0]- X[0,0]) / len(X), edgecolor="white",alpha=0.5,label="forward strand")
	ax.bar(X[:,0], -0.75*X[:,2]/np.sum(X[:,1:]), color="red",width = (X[-1,0]- X[0,0]) / len(X), edgecolor="white",alpha=0.5,label="reverse strand")
	ax.plot(xs, ysf, lw=2.0, color="black", label=r'$p(g|\hat{\Theta}$)')
	ax.plot(xs, ysr, lw=2.0, color="black")
	ax.set_xlabel("Relative Genomic Coordinate")
	ax.set_ylabel("Frequency/Density")
	ax.legend(loc="best",fontsize=16)
	ax.set_yticklabels([0.4,0.3,0.2,0.1,0.0,0.1,0.2,0.3,0.4])
	plt.tight_layout()
	plt.show()



	






