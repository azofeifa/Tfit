import numpy as np
from model import component_bidir
import matplotlib.pyplot as plt
import math
def simulate(N, pi, W, B,X):
	forward 		= np.random.uniform(-W, W,int(N*pi+1))
	reverse 		= np.random.uniform(-W, W, int(N*(1-pi)+1) )
	f_cts,f_edges 	= np.histogram(forward, bins=B, range=(-W,W) ) 	
	r_cts,r_edges 	= np.histogram(reverse, bins=B, range=(-W,W) ) 	
	X[:,0] 			= (f_edges[1:] + f_edges[:-1])/2.
	X[:,1] 			= f_cts
	X[:,2] 			= r_cts
def LOG(x):
	if x > 0:
		return math.log(x)
	else:
		return -100000
def likelihood(X):
	mu 	= X[X.shape[0]/2,0]
	N_reverse 	= np.sum(X[:,2])
	N_forward 	= np.sum(X[:,1])

	N 	= N_reverse + N_forward
	pi 	= N_forward / N
	mean_forward 	= np.average(X[:,0], weights=X[:,1])
	mean_reverse 	= np.average(X[:,0], weights=X[:,2])

	lam =  1./ (0.5*( mean_forward - mean_reverse) )
	if lam < 0:
		lam 	= 1
	var = math.sqrt(np.sum([ pow(X[i,0]-mu,2)*(X[i,1]+X[i,2]) for i in range(X.shape[0]) ])/N)
	vl 		= 1.0 / (X[0,-1] - X[0,0])
	U_ll 	= LOG(vl*pi)*N_forward + LOG(vl*(1-pi))*N_reverse
	U_BIC 	= -2*U_ll  + 1*math.log(N_forward + N_reverse)


	EMG 	= component_bidir(mu, var, lam, 1.0,pi , None,foot_print=0)
	xs 		= np.linspace(X[0,0], X[-1,0], 1000)


	LL 		= sum([ LOG(EMG.pdf(X[k,0],1))*X[k,1]  for k in range(0,X.shape[0]) ])
	LL 		+=sum([  LOG(EMG.pdf(X[k,0],-1))*X[k,2]  for k in range(0,X.shape[0]) ])

	EMG_BIC = -2*LL + 3*math.log(N)

					
	return  U_BIC/EMG_BIC
def bayes_normal(N, W,alpha):
	SI 		= 10
	G 		= {}
	PS 		= {}
	G[0] 	= 1.0
	PS[0] 	= np.random.uniform(-W,W)
	X 		= np.zeros((N,))

	X[0] 	= np.random.normal(PS[0],SI)
	for i in range(1,N):
		KS 	= np.array([G[k]/(i+alpha) for k in G] + [alpha / (i+alpha)])
		KS 	/=sum(KS)
		k 	= np.random.multinomial(1,KS).argmax()
		if k not in G:
			G[k] 	= 0
			PS[k] 	= np.random.uniform(-W,W)
		G[k]+=1
		X[i] 		= np.random.normal(PS[k],SI)

	return X
def bayes_np(N,W,pi,alpha, X, B):
	forward 		= bayes_normal(int(N*pi),W,alpha)
	reverse 		= bayes_normal(int(N*(1.0-pi)),W,alpha)
	f_cts,f_edges 	= np.histogram(forward, bins=B, range=(-W,W) ) 	
	r_cts,r_edges 	= np.histogram(reverse, bins=B, range=(-W,W) ) 	
	X[:,0] 			= (f_edges[1:] + f_edges[:-1])/2.
	X[:,1] 			= f_cts
	X[:,2] 			= r_cts
def simulate_many_np(N,n,W,pi, alpha, X, B ):
	BIC 	= list()
	X 		= np.zeros((B,3))
	for i in range(100):
		bayes_np(n,W,pi, alpha, X, B)
		current 	= likelihood(X)
		if current > 0.5:
			BIC.append(current)

	return BIC
def simulate_many_uniform(N,n,W,pi, alpha, X, B ):
	BIC 	= list()
	X 		= np.zeros((B,3))
	for i in range(100):
		simulate(n,pi,W,B, X)

		current 	= likelihood(X)
		if current > 0.5:
			BIC.append(current)

	return BIC


if __name__ == "__main__":
	W 		= 100
	B 		= 100
	N 		= 400
	alpha 	= 1.0
	pi 		= 0.5
	X 		= np.zeros((B,3))
	X2 		= np.zeros((B,3))
	
	bayes_np(N,W,pi, alpha, X, B)
	simulate(N,0.5,W,B, X2)
	
	F 		= plt.figure()
	ax1 	= F.add_subplot(2,2,1)
	ax1.set_title("Bayes NP Model")
	ax1.bar(X[:,0], X[:,1], edgecolor="w", width=(X[-1,0]- X[0,0])/X.shape[0] )
	ax1.bar(X[:,0], -X[:,2], color="r", edgecolor="w", width=(X[-1,0]- X[0,0])/X.shape[0])
	ax2 	= F.add_subplot(2,2,2)
	ax2.set_title("uniform model")
	ax2.bar(X2[:,0], X2[:,1], edgecolor="w", width=(X2[-1,0]- X2[0,0])/X2.shape[0] )
	ax2.bar(X2[:,0], -X2[:,2], color="r", edgecolor="w", width=(X2[-1,0]- X2[0,0])/X2.shape[0])
	
	BIC_np 	= simulate_many_np(500,N,W,pi, alpha, X, B)
	BIC_uni 	= simulate_many_uniform(500,N,W,pi, alpha, X, B)
	ax3 	= F.add_subplot(2,2,3)
	ax3.set_title("Null of NP")
	ax3.hist(BIC_np,bins=30)

	ax4 	= F.add_subplot(2,2,4)
	ax4.set_title("Null of Uniform")
	ax4.hist(BIC_uni,bins=30)

	plt.show()

	

