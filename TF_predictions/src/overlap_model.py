import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("/Users/azofeifa/Lab/")
sys.path.append("/Users/joazofeifa/Lab/")
from interval_searcher import intervals, node
import numpy as np
import matplotlib.pyplot as plt
import time
import math
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


def prob_single(k,w,l,M):
	p=float(w) / float(l-w)
	vl 	= choose(M,k)*pow(p,k)*pow(1.0-p,M-k  )
		
	return vl

def prob_all(K,w,l,N1, N2):
	return pow(sum([prob_single(k,w,l,N2) for k in range(0, K)]), N1)

def normal(x,mu, si):
	return (1.0 / (math.sqrt(2*math.pi)*si  ))*math.exp(-pow(x-mu,2)/(2*pow(si,2)))


def run(N1,N2,T=100,l=100,alpha=1,beta=1):
	A_stats 	= {}
	B_stats 	= {}
	AOs, BOs 	= list(),list()
	for t in range(T):
		A 	= [(x,x+alpha) for x in np.random.uniform(0, l-alpha, N1)]
		B 	= [(x,x+beta) for x in np.random.uniform(0, l-beta, N2)]
		A.sort()
		B.sort()
		TA 	= node.tree(A)
		TB 	= node.tree(B)
		AO,BO 	= 0,0
		for a_st, a_sp in A:
			FINDS 	= TB.searchInterval((a_st, a_sp))
			if len(FINDS) not in A_stats:
				A_stats[len(FINDS)] 	= 0
			A_stats[len(FINDS)]+=1
			AO+=len(FINDS)

		for a_st, a_sp in B:
			FINDS 	= TA.searchInterval((a_st, a_sp))
			if len(FINDS) not in B_stats:
				B_stats[len(FINDS)] 	= 0
			B_stats[len(FINDS)]+=1
			BO+=len(FINDS)
		AOs.append(AO)
		BOs.append(BO)

	F 	= plt.figure(figsize=(15,10))
	ax1 = F.add_subplot(2,2,1)
	ax1.set_title("List A; N: " + str(N1) )
	ax1.hist([b for b in A_stats],weights=np.array([A_stats[b] for b in A_stats]) / float(sum([A_stats[b] for b in A_stats])), alpha=0.3)
	ax1.scatter([b for b in A_stats],[prob_single(b,alpha+beta, l,N2) for b in A_stats]  )
	ax1.set_xticks(range(0, max(A_stats.keys()) +1) )

	ax2 = F.add_subplot(2,2,2)
	ax2.set_title("List B; N: " + str(N2) )
	ax2.hist([b for b in B_stats],weights=np.array([B_stats[b] for b in B_stats]) / float(sum([B_stats[b] for b in B_stats])), alpha=0.3)
	ax2.set_xticks(range(0, max(B_stats.keys())+1))
	ax2.set_xticklabels([str(i) for i in range(0,max(B_stats.keys())+1)])
	ax2.scatter([b for b in B_stats],[prob_single(b,alpha+beta, l,N1) for b in B_stats]  )
	

	ax3 = F.add_subplot(2,2,3)
	ax3.set_title("Total Number of Overlapping Events on A")
	
	counts,edges 	= np.histogram(AOs,bins=max(AOs)-min(AOs))
	edges 			= edges[:-1]
	counts 			=[float(ct)/float(sum(counts)) for ct in counts]
	
	ax3.bar(edges,counts , alpha=0.3)
	#ax3.set_xticks(range(0, max(AOs) ) )
	xs 	= np.linspace(min(AOs), max(AOs))
	mu 	= np.mean(AOs)
	std = np.std(AOs)
	ax3.scatter(AOs, [ prob_single(a,alpha+beta, l,N2*N1) for a in AOs])
	
	ax3.plot(xs,[ normal(x, mu, std) for x in xs])
	
	ax4 = F.add_subplot(2,2,4)
	ax4.set_title("Total Number of Overlapping Events on B")
	
	counts,edges 	= np.histogram(BOs,bins=max(BOs)-min(BOs))
	edges 			= edges[:-1]
	counts 			=[float(ct)/float(sum(counts)) for ct in counts]
	
	ax4.bar(edges,counts , alpha=0.3)
	#ax3.set_xticks(range(0, max(AOs) ) )
	xs 	= np.linspace(min(BOs), max(BOs))
	mu 	= np.mean(BOs)
	std = np.std(BOs)
	ax4.scatter(BOs, [ prob_single(a,alpha+beta, l,N2*N1) for a in BOs])
	
	ax4.plot(xs,[ normal(x, mu, std) for x in xs])
	



	

	plt.show()

run(20,50,T=30000, l=100, alpha=1, beta=1)


		