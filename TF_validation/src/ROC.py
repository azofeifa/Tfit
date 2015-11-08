import load
import matplotlib.pyplot as plt
import time
import math 
import numpy as np

from scipy.special import erf
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

def normal(x,mu, si):
	return (1.0 / (math.sqrt(2*math.pi)*si  ))*math.exp(-pow(x-mu,2)/(2*pow(si,2)))
def normal_cdf(x,mu,si):
	return 0.5*(1+erf((x-mu)/(si*math.sqrt(2)) ))
def binomial(k,w,l,M):
	p 	= float(2*w) / float(l- 2*w)
	vl 	= choose(M,k)*pow(p,k)*pow(1.0-p,M-k  )
	return vl


def draw_background(ax,l, N):
	ks 	= range(1, 3000)
	w 	= 1500 
	p 	= float(2*w) / float(l- 2*w)
	mu 	= N*p

	si  = math.sqrt(mu*(1-p))

	pdfs 	= np.array([normal(k,mu, si) for k in ks])
	Spdf 	= sum(pdfs)
	ax.plot(ks, [ sum(pdfs[:i])/ Spdf for i in range(len(pdfs)) ] )
	ax.set_xscale("log")


def plot_ROC(A,B, compares=tuple()):
	AG, AB 	= A
	BG , BB = B

	F 	= plt.figure(figsize=(15,10))
	
	#=========================
	#comparing three things
	#1. Bidir and Motif location < D
	#2. Bidir and Motif location > D
	#3. Just Motif 

	BM_T 	 	= list()
	JM  		= list()
	w 		= 1500
	l,n 	= AB[compares[0]]

	p 		= float(2*w) / float(l- 2*w)
	mu 		= n*p
	si  	= 400
	
	for chrom in AG:
		for a in AG[chrom]:
			if a.get_density()>100:
				TF_CT 		= 0
				if compares[0] in a.raw:
					TF_CT 	= a.raw[compares[0]]
				if compares[1] in a.TFS:
					d 	= min([(abs(t),x,y) for t,x,y,z in a.TFS[compares[1]]])

					BM_T.append((math.log(TF_CT+1,10), d[0],  int(1.0-normal_cdf(TF_CT, mu, si) < pow(10,-2) ), d[1], d[2]     )  )
				else:
					d 	= 5000

					BM_T.append((math.log(TF_CT+1,10), d,  int( normal_cdf(TF_CT, mu, si) > 0.999999999 ), 1.0, 1.0     )  )
					
	ax 			= F.add_subplot(2,2,1)
	ax2 		= F.add_subplot(2,2,2)
	ax3 		= F.add_subplot(2,2,3)

	pvs 		= [math.log(pv,10) for x,y,z,pv,q in BM_T if pv]
	min_pv 		= min(pvs)
	max_pv 		= max(pvs)
	ax.set_title("Significant ChIP-seq Signal correlates with smaller i and motif distance")
	
	ax.hist([y for x,y,z,pv,q in BM_T if z==1 and y!=5000 ], bins=100 ,alpha=0.25, label="Significant\nChIP Signal"  )

	ax.hist([y for x,y,z,pv,q in BM_T if z==0 and y!=5000 ], bins=100 ,alpha=0.25, label="Insignificant\nChIP Signal" )
	ax.grid()
	ax.legend()
	
	DS 	= np.linspace(0,5050,1000)
	TPS,FPS = list(),list()

	P 	= float(sum([ 1 for ct, d, ct_p, pv, qv in BM_T if ct_p==1 ]))
	N 	= float(sum([ 1 for ct, d, ct_p, pv, qv in BM_T if ct_p==0 ]))
	
	for D in DS:
		TP 	= 0.0
		FP 	= 0.0
		for ct, d, ct_p, pv, qv in BM_T :
			if d < D and ct_p==1:
				TP+=1
			elif d > D and ct_p==0:
				FP+=1

		TPS.append(TP/P)
		FPS.append(1-(FP/N))

	print 1.0 / len(TPS)
	print 
	print FPS
	

	ax2.plot(FPS, TPS, label="AUC: " + str(sum([TPS[i+1]*(  (FPS[i+1]-FPS[i]))  for i,(x,y) in enumerate(zip(FPS, TPS)) if i +1 < len(TPS) ]) )[:5])
	ax2.plot([0,1],[0,1])
	ax2.legend(loc=(0.8,0.1))
	ax2.set_xlim([0,1])
	ax2.set_ylim([0,1])




	# for D in DS:
	# 	TPS, FPS 	= list(),list()

	# 	pos, pos2, neg 	= 0.0,0.0,0.0
		

	# 	for tf_ct, d, hit in BM_T:
	# 		if d <= D and hit:
	# 			pos+=1.0
	# 		elif d > D and not hit:
	# 			neg+=1
	# 	TPS.append(pos/POSITIVES)
	# 	FPS.append( 1-(neg/NEGATIVES ))
	# 	ax2.plot(FPS, TPS)


	#ax2.plot([0,1],[0,1])

	plt.show()












if __name__ == "__main__":
	OUT1 		= "/Users/joazofeifa/Lab/TF_predictions/distances_crude/Allen2014_rawCHIP_motif_distances.tsv"	
	OUT2 		= "/Users/joazofeifa/Lab/TF_predictions/FIMO_OUT/MAX/fimo_rawCHIP.txt"
	OUT3 		= "/Users/joazofeifa/Lab/TF_predictions/FIMO_OUT/Sp1/fimo_rawCHIP.txt"
	TEST 		= False
	DMSO2_3 	= load.load_bidir(OUT1, RAW=True, test=TEST)
#	MAX 		= load.read_in_motifs(OUT2, RAW=True, test=TEST)
	SP1 		= load.read_in_motifs(OUT3, RAW=True, test=TEST)


	plot_ROC(DMSO2_3,SP1 , compares=("SP1","SP1_f1"  ) )