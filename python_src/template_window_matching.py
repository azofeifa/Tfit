import load,matplotlib.pyplot as plt
import numpy as np
import model

import math
def LOG(x):
	if x >0:
		return math.log(x)
	return -np.inf
def window(X, std=10, lam=0.1, step_size=1, norm_to_max=True):
	i 				= 0
	win 			= (std + 1.0 / lam)
	mns 			= list()
	forward,reverse = list(), list()
	while i < X.shape[0]:
		j,k 	= i,i
		while j < X.shape[0] and (X[j,0] - X[i,0]) < win:
			j+=1
		while k >=0 and (X[i,0] - X[k,0]) < win:
			k-=1
		forward.append( np.sum(X[i:j, 1]) )
		reverse.append( np.sum(X[k:i, 2]) )
		i+=1
	scores 	= [LOG(x) + LOG(y) for x,y in zip(forward, reverse)]
	
	if norm_to_max:
		scores 	= [x / max(scores) for x in scores]
	return scores
	
def bayes_factor(X, std=10, lam=0.1,step_size=1, norm_to_max=True):
	KS 	= list()
	i 	= 0
	win 			= (3*std + 1.0 / lam)
	M1S 			= list()
	M2S 			= list()
	while i < X.shape[0]:
		j,k 	= i,i
		EMG 	= model.component_bidir(X[i,0], std, lam, 1.0, 0.5, None)

		while j < X.shape[0] and (X[j,0] - X[i,0]) < win:
			j+=1
		while k >=0 and (X[i,0] - X[k,0]) < win:
			k-=1
		if j < X.shape[0] and k>=0:
			M1 		= sum([LOG(x)*y for x,y in zip(map(lambda x: EMG.pdf(x,1), X[k:j,0]), X[k:j,1])])
			M1 		+=sum([LOG(x)*y for x,y in zip(map(lambda x: EMG.pdf(x,-1), X[k:j,0]), X[k:j,2])])
			l 		= 1 / (X[-1,0] - X[0,0])
			M2 		= sum([LOG(l)*y for y in X[k:j,1]]) + sum([LOG(l)*y for y in X[k:j,2]])
			if M2==0:
				KS.append(1)
			else:
				KS.append(M2/M1)
			M1S.append(M1/np.sum(X[:,1:]))
			M2S.append(M2/np.sum(X[:,1:]))
		i+=1
	if norm_to_max:
		KS 	= [k/ max(KS) for k in KS]
	return KS
def find_peaks(data):
	peaks =list()
	for i in range(1, len(data)-1):
		if data[i-1][1] < data[i][1] > data[i+1][1]:
			peaks.append((data[i][1], data[i][0]))
	peaks.sort()
	peaks 	= peaks[::-1]
	peaks 	= [(y,x) for x,y in peaks]
	return peaks

def center(coverage_scores, bayes_ks):
	return [(x+y)/2. for x,y in zip(coverage_scores, bayes_ks)]
def draw(X, coverage_scores, bayes_ks, hybrid,starts):
	F 	= plt.figure(figsize=(15,10))
	ax1 = F.add_subplot(2,1,1)
	ax1.bar(X[:,0], X[:,1]/np.sum(X[:,1:]))
	ax1.bar(X[:,0], -X[:,2]/np.sum(X[:,1:]))
	ax1.grid()
	
	ax2 = F.add_subplot(2,1,2)
	ax2.bar([x for x,y in starts[:3]], [y for x,y in starts[:3]] )
	ax2.plot(np.linspace(X[0,0], X[-1,0], len(coverage_scores)), coverage_scores, label="Log Coverage Score" )
	ax2.plot(np.linspace(X[0,0], X[-1,0], len(bayes_ks)), bayes_ks, label="Bayes Factor" )
	ax2.plot(np.linspace(X[0,0], X[-1,0], len(hybrid)), hybrid, label="center", linestyle="--")
	ax2.grid()
	
	plt.show()

def compute_possible_EM_starts(X, std=1, lam=0.1):
	coverage_scores 	= window(X,std=std,lam=lam,step_size=1)
	bayes_ks 			= bayes_factor(X,std=std,lam=lam,step_size=1)
	hybrid 				= center(coverage_scores, bayes_ks)
	starts 				= find_peaks([(x,y) for x,y in zip(np.linspace(X[0,0], X[-1,0], len(hybrid)), hybrid)])
	return coverage_scores, bayes_ks, hybrid, starts
def sample(X, k, std=1, lam=0.1):
	coverage_scores, bayes_ks, hybrid, starts 	= compute_possible_EM_starts(X,std=std, lam=lam)
	keeps 			= list()
	for i in range(k):
		j 			= np.random.geometric(0.8)-1
		keeps.append(starts[j][0])
		starts 		= starts[:j] + starts[j+1:]
	return keeps




	
	

if __name__=="__main__":
	X 	= load.grab_specific_region("chr1",6229860,6303055, SHOW=False, bins=500 )
	X[:,0]/=100.
	X[:,0]-=X[0,0]

	coverage_scores, bayes_ks, hybrid, starts 	= compute_possible_EM_starts(X,std=1,lam=0.1)
	#draw(X, coverage_scores, bayes_ks, hybrid,starts)
	clf = model.EMGU(noise=True, K=3,noise_max=0.01,
		moveUniformSupport=5,
		max_it=50,seed=True)
	clf.fit(X)
	clf.draw(X)
	
	