import numpy as np
import matplotlib.pyplot as plt
import time, sys
sys.path.append("/Users/azofeifa/Lab/")
sys.path.append("/Users/joazofeifa/Lab/")

from interval_searcher import intervals, node

def simulate(L1=10, L2=100,  N1=100, N2=100, N3=100, D=1000000 , sigma=1 ):
	xe, xm 		= 0,0
	motif, eRNA = list(),list()
	motif, eRNA = np.random.uniform(0, D,  N1 ),np.random.uniform(0, D, N2 )
	motif, eRNA = [(x-L1, x+L1) for x in motif],[(x-L2, x+L2) for x in eRNA]
	coocur_pts 	= np.random.uniform(0,D, N3 )
	coocur_pts 	= [(x, x+np.random.normal(0, sigma)) for x in coocur_pts]
	motif, eRNA = motif + [(y - L1, y+ L1) for x,y in coocur_pts], eRNA + [(x - L2, x+ L2) for x,y in coocur_pts]
	
	return motif, eRNA

def compute_overlaps(x,y,OUT=""):
	y.sort()
	T 	= node.tree(y)
	j,N,O 	= 0,len(y),list()
	#FHW.open(OUT,"w")
	for start, stop in x:
		FINDS 	= T.searchInterval((start, stop))
		cx 		= (stop+start) / 2.
		for st,sp in FINDS:
			cy  = (sp+st) / 2.
			O.append((cx-cy))
	return O
def show():
	F 	= plt.figure(figsize=(15,10))
	ax1 = F.add_subplot(2,2,1)
	ax2 = F.add_subplot(2,2,2)
	ax3 = F.add_subplot(2,2,3)
	ax4 = F.add_subplot(2,2,4)
	
	N1,N2,N3,sigma 	= 5000,5000,0,1	
	ax1.set_title("N1: "+ str(N1) + ", N2: "+ str(N2)+ ", N3: "+ str(N3) + ", sigma: " + str(sigma))
	arrivals_motif,arrivals_eRNA 	= simulate(L1=10, L2=10, N1=N1,N2=N2,N3=N3, sigma=sigma)
	D 	= compute_overlaps(arrivals_motif, arrivals_eRNA)
	ax1.hist(D, bins=100)
	ax1.grid()

	N1,N2,N3,sigma 	= 100,1000,100,1	
	ax2.set_title("N1: "+ str(N1) + ", N2: "+ str(N2)+ ", N3: "+ str(N3) + ", sigma: " + str(sigma))
	arrivals_motif,arrivals_eRNA 	= simulate(L1=10, L2=10, N1=N1,N2=N2,N3=N3, sigma=sigma)
	D 	= compute_overlaps(arrivals_motif, arrivals_eRNA)
	ax2.hist(D, bins=100)
	ax2.grid()
	
	N1,N2,N3,sigma 	= 10000,1000,10000,5
	ax3.set_title("N1: "+ str(N1) + ", N2: "+ str(N2)+ ", N3: "+ str(N3) + ", sigma: " + str(sigma))
	arrivals_motif,arrivals_eRNA 	= simulate(L1=10, L2=10, N1=N1,N2=N2,N3=N3, sigma=sigma)
	D 	= compute_overlaps(arrivals_motif, arrivals_eRNA)
	ax3.hist(D, bins=100)
	ax3.grid()

	N1,N2,N3,sigma 	= 10,1000,100,1	
	ax4.set_title("N1: "+ str(N1) + ", N2: "+ str(N2)+ ", N3: "+ str(N3) + ", sigma: " + str(sigma))
	arrivals_motif,arrivals_eRNA 	= simulate(L1=10, L2=10, N1=N1,N2=N2,N3=N3, sigma=sigma)
	D 	= compute_overlaps(arrivals_motif, arrivals_eRNA)
	ax4.hist(D, bins=100)
	ax4.grid()
	
	plt.show()
	

show()