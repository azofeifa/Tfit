import load
from model import component_bidir
import math
import numpy as np
import matplotlib.pyplot as plt
def LOG(x):
	if x!= 0:
		return math.log(x)
	return -np.inf
def get_sigma(X, mu, i,j, lam,foot_print=0):
	forward=0
	reverse=0
	N_f 	= 0
	N_r 	= 0
	for k in range(i,j):
		N_f+=X[k,1]
		N_r+=X[k,2]
		forward+=pow(X[k,0]-mu-foot_print,2 )*X[k,1]
		reverse+=pow(X[k,0]-mu+foot_print,2 )*X[k,2]
	if N_f and N_r:
		si 	=  0.5*((math.sqrt(forward/N_f)) +  math.sqrt(reverse/N_r))
		return si - (1. / lam)
	return None
def get_lambda(X, i,j,foot_print=0):
	forward=0
	reverse=0
	N_f 	= 0
	N_r 	= 0
	for k in range(i,j):
		N_f+=(X[k,1])
		N_r+=(X[k,2])
		forward+=(X[k,0]-foot_print)*X[k,1]
		reverse+=(X[k,0]+foot_print)*X[k,2]
	if N_f and N_r:
		return 1.0 / (0.5*((forward/N_f) - (reverse/N_r) ));
	return None

def compute_ll(X, i,j, mu,si, l, pi, SHOW=False, foot_print=0):
	EMG 	= component_bidir(mu, si, l, 1.0,pi , None,foot_print=foot_print)
	LL 		= sum([ LOG(EMG.pdf(X[k,0],1))*X[k,1]  for k in range(i,j) ])
	LL 		+=sum([ LOG(EMG.pdf(X[k,0],-1))*X[k,2]  for k in range(i,j) ])
	if SHOW:
		plt.bar(X[i:j,0], X[i:j,1]/np.sum(X[i:j,1]))
		plt.bar(X[i:j,0], -X[i:j,2]/np.sum(X[i:j,2]) )
		xs 	= np.linspace(X[i,0],X[j,0],100)
		plt.plot(xs, [EMG.pdf(x,1) for x in xs])
		plt.plot(xs, [-EMG.pdf(x,-1) for x in xs])
		
		plt.show()
		
	return LL
	pass


def run_MM(X, window=500, scale=100):
	ratios_x, ratios_y, densities_s 	= list(),list(),list()
	foot_print=0
	for w in (1500,):
		window 	= w / scale
		ratio_x 	= list()
		ratio_y 	= list()
		densities 	= list()
		for i in range(X.shape[0]):
			#get j,k
			j 	= i
			k 	= i
			while j > 0 and abs(X[j,0]-X[i,0]) < window:
				j-=1
			while k <(X.shape[0]-1) and abs(X[k,0]-X[i,0]) < window:
				k+=1

			mu 	= X[i,0]
			lam = get_lambda(X, j,k,foot_print=foot_print)
			si 	= get_sigma(X, mu, j,k, lam,foot_print=foot_print)
			print lam, si, np.sum(X[j:k,1:]) ,j,k
			if (lam!= None and si!= None and lam >0 and si > 0 and np.sum(X[j:k,1:]) ):
				EMG_ll 	=  compute_ll(X, j,k, mu,si, lam, 0.5,foot_print=foot_print)
				EMG_BIC = -2*EMG_ll + 3*math.log(np.sum(X[j:k,1:]))

				vl 		= 1.0 / (X[k,0]-X[j,0])
				pi 		= np.sum(X[j:k, 1])/ np.sum(X[j:k, 1:])
				U_ll 	= LOG(vl*pi)*np.sum(X[j:k, 1]) + LOG(vl*(1-pi))*np.sum(X[j:k, 2])
				U_BIC 	= -2*U_ll + 1*math.log(np.sum(X[j:k, 1:]))
				ratio_x.append(X[i,0])
				ratio_y.append(U_BIC/EMG_BIC)
				densities.append(vl*np.sum(X[j:k, 1]) + vl*np.sum(X[j:k,2]) )
			else:
				ratio_x.append(X[i,0])
				ratio_y.append(0)
				densities.append(0 )
				
		densities_s.append(densities)
		ratios_x.append(ratio_x)
		ratios_y.append(ratio_y)
	F 	= plt.figure(figsize=(15,10))
	ax1 = F.add_subplot(3,1,1)
	ax2 = F.add_subplot(3,1,2)
	ax3 = F.add_subplot(3,1,3)
	ax1.bar(X[:,0],X[:,1],width=(X[-1,0]-X[0,0])/X.shape[0])
	ax1.bar(X[:,0],-X[:,2],width=(X[-1,0]-X[0,0])/X.shape[0])
	ax1.set_xlim(X[0,0],X[-1,0])
	ax1.grid()
	colors 	= ("r", "g", "b", "black")
	for i,(ratio_x, ratio_y, densities)  in enumerate(zip(ratios_x,ratios_y,densities_s)):

		ax2.scatter(ratio_x, ratio_y,c=colors[i])
		ax2.set_xlim(X[0,0],X[-1,0])
		ax2.grid()
		ax3.scatter(ratio_x, densities,c=colors[i])
		ax3.set_xlim(X[0,0],X[-1,0])
		ax3.grid()
		
	plt.show()












if __name__ == "__main__":

	IN  = "/Users/joazofeifa/Lab/gro_seq_files/Danko_2014/"
	#chr1:87,691,254-87,695,004
	#88,319,575-88320266
	#92,308,146-92,315,100
	X 	= load.grab_specific_region("chr1",88319575,88320266, SHOW=False, bins=500, 
		pos_file=IN+"GSE66031_ac16.unt.all_plus.bw.bedGraph.mod.bedGraph", neg_file=IN+"GSE66031_ac16.unt.all_minus.bw.bedGraph.mod.bedGraph" )
	X[:,0]-=min(X[:,0])
	scale = 100
	window = 1500
	X[:,0]/=scale
	run_MM(X, window=window, scale=scale )
	