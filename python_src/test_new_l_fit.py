import simulate
from model import component_bidir, component_elongation
import numpy as np
import matplotlib.pyplot as plt
import load
import math as m
def grab_data(forward, reverse):
	D 	= (list(),list())
	for i,FILE in enumerate((forward, reverse)):
		with open(FILE) as FH:
			for line in FH:
				chrom,start, stop, cov 	= line.strip("\n").split("\t")
				x 						= (int(stop) + int(start ))/2.
				D[i].append((x, float(cov)))
	X 	= simulate.BIN(D[0], D[1], 200)
	X[:,0] 	-= X[0,0]
	X[:,0] /=100.
	return X
def LOG(x):
	if x ==0:
		return -np.inf
	return m.log(x)
def prior_l_norm(l, mu, si):

	return 1.0 / m.sqrt(2*m.pi)


def move_l(X, parameters,penality=2):
	forward_as 	= [mu+si + (1.0 / l) for mu, si, l, w, pi in parameters]
	reverse_bs 	= [mu-si - (1.0 / l) for mu, si, l, w, pi in parameters]
	forward_is 	= list()
	reverse_is 	= list()
	for f in range(len(parameters)):
		for i in range(X.shape[0]):
			if X[i,0] > forward_as[f]:
				break
		for j in range(X.shape[0]):
			if X[j,0] > reverse_bs[f]:
				break
		forward_is.append(i)	
		reverse_is.append(j)
	for f,i in enumerate(forward_is):
		lls 	= list()
		if f < len(forward_is)-1:
			j 	= reverse_is[f+1]
		else:
			j 	= X.shape[0]-1
		N 		= sum(X[i:j,1])
		oll 		= sum(X[i:j,1])*m.log(1. / (X[j,0]-X[i,0] ))
		NULL 	= -2*oll + LOG(N)
		for l in range(i,j):
			N_left , N_right 	= sum(X[i:l, 1]),sum(X[l:j, 1])
			w_left 				= N_left / (N_right + N_left)
			vl_left 			= LOG(w_left/(X[l,0] - X[i,0]))
			vl_right 			= LOG((1.0-w_left)/(X[j,0] - X[l,0]))
		
			ll 					= vl_left*N_left + vl_right*N_right
			print oll, ll
			BIC 				= -2*ll + penality*LOG(N)
			lls.append(NULL/BIC)
		F 	= plt.figure(figsize=(15,10))
		ax1 	= F.add_subplot(311)
		ax2 	= F.add_subplot(312)
		ax3 	= F.add_subplot(313)
		ax1.bar(X[:,0], X[:,1])
		ax1.bar(X[:,0], -X[:,2])
		ax1.scatter([X[i,0],X[j,0]], [0,0], s=30 )

		ax2.bar(X[i:j, 0], X[i:j, 1])
		ax2.grid()
		ax3.plot(X[i:j, 0], lls)
		ax3.grid()
		plt.show()
	for f,j in enumerate(reverse_is):
		lls 	= list()
		if f == 0:
			i 	= 0
		elif f>0:
			i 	= forward_is[f-1]
		N 		= sum(X[i:j,2])
		oll 		= sum(X[i:j,2])*m.log(1. / (X[j,0]-X[i,0] ))
		NULL 	= -2*oll + LOG(N)
		for l in range(i,j):
			N_left , N_right 	= sum(X[i:l, 2]),sum(X[l:j, 2])
			w_left 				= N_left / (N_right + N_left)
			vl_left 			= LOG(w_left/(X[-1,0] - X[0,0]))
			vl_right 			= LOG((1.0-w_left)/(X[j,0] - X[l,0]))
			ll 					= vl_left*N_left + vl_right*N_right
			BIC 				= -2*ll + penality*LOG(N)
			lls.append(NULL/BIC) 
		F 	= plt.figure(figsize=(15,10))
		ax1 	= F.add_subplot(311)
		ax2 	= F.add_subplot(312)
		ax3 	= F.add_subplot(313)
		ax1.bar(X[:,0], X[:,1])
		ax1.bar(X[:,0], -X[:,2])
		ax1.scatter([X[i,0],X[j,0]], [0,0], s=30 )
		ax3.grid()

		ax2.bar(X[i:j, 0], -X[i:j, 2])
		ax2.grid()
		ax3.plot(X[i:j, 0], lls)
		plt.show()
	



	



if __name__=="__main__":
	# D 	= simulate.runOne(mu=0, s=1, l=5, lr=100, ll=-100, we=0.5,wl=0.25, 
	# 	wr=0.25, pie=0.5, pil=0.1, pir=0.9, N=1000, SHOW=False , bins=200, noise=True )
	forward 	= "/Users/joazofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/chr1_6229860_6303055_DMSO2_3.pos.BedGraph"
	reverse 	= "/Users/joazofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/chr1_6229860_6303055_DMSO2_3.neg.BedGraph"
	X 			= grab_data(forward, reverse)
	parameters 	= [	(298.119448,1.036515,0.228770,0.111111,0.423449), 
					(358.097039,3.538068,0.129734,0.111111,0.391122), 
					(660.517058,0.063181,0.187830,0.111111,0.379944)]
	move_l(X, parameters, penality=1000)












