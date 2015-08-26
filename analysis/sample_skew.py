import matplotlib.pyplot as plt
import math as m 
import numpy as np
import time
def load_intervals(FILE,pad):
	G 	= {}
	with open(FILE) as FH:
		for line in FH:
			chrom,start, stop 	= line.split("\t")[:3]
			start, stop 		= int(start), int(stop)
			if chrom not in G:
				G[chrom] 		= list()
			G[chrom].append((start-pad, stop+pad, list(), list() ))
	for chrom in G:
		G[chrom].sort()
	return G

def insert(FILE, i, G):
	prevchrom=None
	t 		= 0
	with open(FILE) as FH:
		for line in FH:
			chrom,start, stop, cov 	= line.strip("\n").split("\t")
			x, cov 					= (float(start) + float(stop)) / 2., float(cov)
			if (chrom != prevchrom):
				j 	= 0
				if chrom in G:
					N 	= len(G[chrom])
				else:
					N 	= 0
				t+=1
				if (t > 1):
					break
			while j < N and G[chrom][j][1] < x:
				j+=1
			if (j < N and G[chrom][j][0]<= x <=G[chrom][j][1]  ):
				G[chrom][j][i].append((x,cov))
			prevchrom=chrom
	return G

def sample_mean(D):
	S,N 	= 0,0
	for x,y in D:
		S+=(x*y)
		N+=y
	if N:
		return S/N
	return None
def sample_mean2(f,r):
	F,R 	= sample_mean(f), sample_mean(r)

	if (F is None or R is None):
		return None
	PI 		= sum([y for x,y in f]) / (sum([y for x,y in f]) + sum([y for x,y in r]))  
	PI 		= 0.5
	return sample_mean(f)*PI + sample_mean(r)*(1.-PI) 
def sample_variance(D, mu):
	
	S,N 	= 0,0
	for x,y in D:
		S+=pow(x-mu,2)*y
		N+=y
	return S/N

def sample_skew(D,mu):
	if mu is None:
		return None
	sigma 		= sample_variance(D,mu)
	if sigma is None or sigma==0:
		return None
	S,N 	= 0,0
	for x,y in D:
		S+=pow(x-mu,3)*y
		N+=y

	return (S/N)/pow(sigma, 3./2.)
def density(f,r, st, sp, hamming=False,mu=None):
	S1,S2 	= 0.,0.
	N1,N2 	= 0.,0.
	if mu is None:
		mu 	= sample_mean2(f,r)
	for i,D in enumerate((f,r)):
		for x,y in D:
			if x > mu :
				if hamming and y>0:
					S1+=1
				else:
					S1+=y
				N1+=1
			else:
				if hamming and y>0:
					S2+=1
				else:
					S2+=y
				N2+=1
	try:
		vl1, vl2 	= (S1/(sp-st)),(S2/(sp-st))
		if (vl1==0 or vl2==2):
			return None,None
		return vl1, vl2
	except:
		return None, None
def density2(f, r, st, sp):
	S=0
	for D in (f,r):
		for x,y in D:
			S+=y

	return m.log(S/(sp-st))
	pass
def N_DIST(x, mu, sigma):
	return (1.0 / m.sqrt(2*m.pi)*sigma)*m.exp(-pow(x-mu,2)/(2*pow(sigma,2)))

def GMM(D, Y,K):
	minX, maxX 	= min(D), max(D)

	mus 		= [np.random.uniform(minX, maxX) for k in range(K)]
	sigmas 		= [10. for k in range(K)]
	pis 		= [1. / K for k in range(K) ]
	R 			= np.zeros((len(D), K))
	t 			= 0
	T 			= 1000
	while t < T:
		for k in range(K):
			mu,sig,pi 		= mus[k], sigmas[k], pis[k]
			for i,x in enumerate(D):
				R[i,k] 	= N_DIST(x, mu, sig)*pi
		#normalize
		for i in range(len(D)):
			N  	= sum(R[i,:])
			if N:
				R[i,:]/=N
		#new parameters
		PI 				= 0
		for k in range(K):
			WN 			= sum(R[:,k]*Y)
			mus[k] 		= sum([x*R[i,k]*Y[i] for i, x in enumerate(D)])/WN
			sigmas[k] 	= m.sqrt(sum([pow(x-mus[k],2)*R[i,k]*Y[i] for i, x in enumerate(D)])/WN)
			PI 			+=WN
		for k in range(K):
			pis[k] 		= (sum(R[:,k]*Y) )/(PI )


		print t,pis
		t+=1
	plt.hist(D, weights=Y, bins=70, normed=1, alpha=0.5)
	plt.plot(np.linspace(minX, maxX, 1000), [ sum([ N_DIST(x, mus[k], sigmas[k])* pis[k] for k in range(K)])  for x in np.linspace(minX, maxX, 1000)],linewidth=2. )
	plt.show()
	
def histogram(G, pad):
	ss_f 		= [sample_skew(f,sample_mean2(f,r)) for chrom in G for st,sp, f, r in G[chrom] if sample_skew(f,sample_mean2(f,r)) is not None and sample_skew(r,sample_mean2(f,r)) is not None   ]
	ss_r 		= [sample_skew(r,sample_mean2(f,r)) for chrom in G for st,sp, f, r in G[chrom] if sample_skew(r,sample_mean2(f,r)) is not None and sample_skew(f,sample_mean2(f,r)) is not None  ]
	counts,edges = np.histogram(ss_f+ss_r, 550)
	edges 		= (edges[:-1] + edges[1:])/2.
	#GMM(edges, counts, 3)
	lens 		= [abs(y-x-pad*2) for chrom in G for x,y ,f,r in G[chrom] if (abs(y-x) < 7000)]
	densities 	= [ (density(f, r, st, sp),sample_skew(f,sample_mean2(f,r)), sample_skew(r,sample_mean2(f,r)))  for chrom in G for st, sp, f, r in G[chrom] if len(f) and len(r) ]
	densities_a = [ (density2(f,r, st, sp), sample_skew(f,sample_mean2(f,r) )) for chrom in G for st, sp, f, r in G[chrom] if len(f) and len(r) and sample_skew(f,sample_mean2(f,r)) is not None ]
	densities_b = [ (density2(f,r, st, sp), sample_skew(r,sample_mean2(f,r) )) for chrom in G for st, sp, f, r in G[chrom] if len(f) and len(r) and sample_skew(r,sample_mean2(f,r)) is not None ]

	F 			= plt.figure(figsize=(15,10))
	ax1 		= F.add_subplot(2,2,1)
	ax2 		= F.add_subplot(2,2,2)
	ax3 		= F.add_subplot(2,2,3)
	ax4 		= F.add_subplot(2,2,4)
	
	ax1.hist(ss_f, alpha=0.6, label="forward", bins=75, normed=1)
	ax1.hist(ss_r, alpha=0.6, label="reverse", bins=75, normed=1)
	ax1.set_xlabel("Sample Skew")
	ax1.set_ylabel("Frequency")
	ax1.grid()
	ax1.legend()
	DDD=[ sample_mean2(f,r) - ((sp+st)/2.) for chrom in G for st, sp, f,r in G[chrom] if len(f) and len(r)]
	ax2.hist(DDD ,bins=50, label="Mean: " + str(np.mean(DDD) ) + "\nStd: " + str(np.std(DDD)) )
	ax2.legend()
	ax2.set_xlabel("Displacement of Model Mean from Center of Peak Call")
	ax2.set_ylabel("Frequency")
	ax2.grid()

	ax3.scatter([m.log(x) for (x ,y), sf,sr in densities  if x is not None and y is not None ],[m.log(y) for (x,y), sf,sr in densities if x is not None and y is not None  ],alpha=0.2, color="blue")
	ax3.set_xlabel("Forward Strand Density")
	ax3.set_ylabel("Reverse Strand Density")


	ax3.grid()

	ax4.grid()
	ax4.scatter([x for x,y in densities_a],[y for x,y in densities_a],alpha=0.2, color="blue", label="Forward")
	ax4.scatter([x for x,y in densities_b],[y for x,y in densities_b],alpha=0.2, color="red", label="Reverse")
	# # ax4.scatter([ pow(x-y,2) for (x ,y), sf,sr in densities ],[ pow(sf-sr,2) for (x ,y), sf,sr in densities ])
	

	

	ax4.set_xlabel("Read Density")
	ax4.set_ylabel("Skew")
	ax4.legend()


	# plt.tight_layout()
	plt.show()
	pass

def load(FILE):
	G 	= {}
	with open(FILE) as FH:
		for line in FH:
			chrom,start, stop, cov 	= line.strip("\n").split("\t")
			x,cov 					= (float(stop) + float(start)) / 2. , float(cov)
			if chrom not in G:
				G[chrom]	= list()
			G[chrom].append((x,cov))
	return G

def BIN2(f,r,res, bins=None, normed=False):
	if not len(f) or not len(r):
		return None
	ab 			= f, r
	minX,maxX 	= min(min(ab[0])[0], min(ab[0])[0]),max(max(ab[1])[0], max(ab[1])[0])
	if bins is None:
		bins 		= (maxX-minX)/float(res)
	X 			= np.zeros((bins, 3))
	X[:,0] 		= np.linspace(minX, maxX, int(bins))
	for i,d in enumerate(ab):
		j 		= 0
		for x,y in d:
			while j < bins and X[j,0] < x:
				j+=1
			if (j < bins and x < X[j,0]):
				X[j,i+1]+=y
	
	return X

def BIN(f,r,res):
	G 	= {}
	for chrom in f:
		if chrom in r:
			
			G[chrom] 	= BIN2(f[chrom], r[chrom], res)
	return G

			



def search(G, delta=10, window=1000, skew=0.5):
	A 	= {}
	FHW 	= open("/Users/joeyazo/Desktop/Lab/gro_seq_files/HCT116/EMG_out_files/bidirectional_hits_intervals.bed","w")

	for chrom in G:
		X 		= G[chrom]
		A[chrom] 	= list()
		j, N 	= 0, X.shape[0]
		#first thing look at density corratlations
		start, stop 	= 0,0
		for i in range(N):
			lower, upper 	= max(X[i,0]-window, X[0,0]), min(X[i,0]+window, X[-1,0])
			while j < N and X[j,0] < lower:
				j+=1
			k 	= j
			f,r 	= list(), list()
			while k < N and X[k,0] < upper:
				f.append((X[k,0], X[k,1]))
				r.append((X[k,0], X[k,2]))
				k+=1
			if X[i,0] < upper and X[i,0] > lower:
				df, dr 	= sum([y for x,y in f if x > X[i,0] and y > 0 ])/float(upper-X[i,0]),sum([y for x,y in r if x < X[i,0] and y > 0 ])/float(X[i,0]-lower)
				if df > 5. and dr > 5.:
					mean 	= sample_mean2(f,r)
					if mean is not None and abs(mean-X[i,0]) < 50:
						sf,sr 	= sample_skew(f, mean), sample_skew(r, mean)
						if (sf>0 and sr < 0):
							plt.title("df: " + str(df) + ", dr: " + str(dr))
							plt.scatter([mean], [0])
							plt.bar([x for x,y in f], [y for x,y in f])
							plt.bar([x for x,y in r], [-y for x,y in r])
							plt.show()
								
	FHW.close()



if __name__ == "__main__":
	forward 	= "/Users/joeyazo/Desktop/Lab/gro_seq_files/HCT116/bed_graph_files/test_DMSO2_3.pos.BedGraph"
	reverse 	= "/Users/joeyazo/Desktop/Lab/gro_seq_files/HCT116/bed_graph_files/test_DMSO2_3.neg.BedGraph"
	ESTIMATE 	= True
	SEARCH 		= False
	if ESTIMATE:
		intervals 	= "/Users/joeyazo/Desktop/Lab/ENCODE/HCT116/H3K27ac/peak_files/HCT-116_H3K27Ac_narrowPeak.bed"
		intervals 	= "/Users/joeyazo/Desktop/Lab/gro_seq_files/HCT116/EMG_out_files/bidirectional_hits_intervals.bed"
		pad 		= 0
		G 			= load_intervals(intervals, pad)
		G 			= insert(forward, 2, G)
		G 			= insert(reverse, 3, G)
		histogram(G,pad)
	if SEARCH:
		f 			= load(forward)
		r 			= load(reverse)
		G 			= BIN(f, r, 10)
		search(G, delta=50, window=500, skew=0.5)

	