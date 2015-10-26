import expression_correlation as ec
import time
import numpy as np
from scipy.stats import hypergeom
from operator import mul    # or mul=lambda x,y:x*y
from fractions import Fraction
import math
from math import factorial as fac
from scipy import special
class dnase:
	def __init__(self,chrom, start, stop, info, FILTER):
		self.chrom,self.start, self.stop 	= chrom, start, stop
		self.TFS 	= {}
		self._parse_info(info, FILTER)
		self.eRNA 	= False
	def _parse_info(self, info, FILTER):
		for TF in info.split(","):
			tf, pval 	= TF.split(":")
			tf 			= tf.strip("_fimo_out")
			pval 		= float(pval)
			tf 			= tf.split("_")[0]
				
			if tf in FILTER :
				if tf not in self.TFS:
					self.TFS[tf]=list()
				self.TFS[tf].append(pval)

def load_DNAse(FILE, FILTER={}, test=False):
	G 	= {}
	t 	= 0
	with open(FILE) as FH:
		for line in FH:
			if test and t > 10000:
				break
			chrom,start, stop, info = line.strip("\n").split("\t")
			if chrom not in G:
				G[chrom]=list()
			G[chrom].append((int(start), dnase(chrom,int(start),int(stop), info, FILTER)))
			t+=1
	for chrom in G:
		G[chrom].sort()
		G[chrom] 	= [y for x,y in G[chrom]]

	return G

def label_DNase(B,D, attr="eRNA"):
	O 	= 0
	for chrom in D:
		if chrom in B:
			b,DNASE 	= B[chrom], D[chrom]
			j,N 		= 0,len(b)
			for d in DNASE:
				while j < N and b[j].stop < d.start:
					j+=1
				if j < N and b[j].start < d.stop:
					O+=1
					setattr(d, attr, True)
def write_out(D, B, OUT, TFS):
	FHW 	= open(OUT, "w")
	FHW.write("#Marginals\n")
	#write out marginals...
	M 		= dict([(tf, [0, 0]) for tf in TFS])
	MM 		= dict([(tf, [0, 0]) for tf in TFS])
	for chrom in D:
		for d in D[chrom]:
			for tf in d.TFS:
				if not d.eRNA:
					M[tf][0]+=1
				else:
					M[tf][1]+=1
	for chrom in B:
		for b in B[chrom]:
			for TF in b.TFS:
				for tf in TF.split(","):
					tf 	= tf.split("_")[0]
					if not b.DNAse:
						MM[tf][0]+=1
					else:
						MM[tf][1]+=1
						

	for tf in M:
		FHW.write(tf + "\t" + str(M[tf][0])+ "," + str(M[tf][1]) + "\t" + str(MM[tf][0])+","+str(MM[tf][1])+ "\n"  )
	#okay now we want pairwise
	IDS 	= dict([ (tf, i) for i,tf in enumerate(TFS)])
	FHW.write("#" + ",".join([str(IDS[tf]) + ":" + tf for tf in IDS ])+ "\n")
	A 		= np.zeros((len(IDS), len(IDS), 2))
	for chrom in D:
		for b in D[chrom]:
			current 	= b.TFS.keys()
			for i  in range(len(current)) :
				for j in range(i+1, len(current)) :
					u,v 	= IDS[current[i]], IDS[current[j]]
					if not b.eRNA:
						A[u,v,0]+=1
						A[v,u,0]+=1					
					else:
						A[u,v,1]+=1
						A[v,u,0]+=1
						
	for u in range(A.shape[0]):
		line=""
		for v in range(A.shape[1]):
			line+=str(A[u,v][0])+":"+ str(A[u,v][1])+","
		line=line.strip(",")
		FHW.write(line+"\n")
	
def load_counts(FILE):
	M 	= {}
	MARGINAL 	= True
	L 			= list()
	with open(FILE) as FH:
		header=True
		for line in FH:
			if line[0]=="#" and 'Marginals' not in line:
				toIDS 	= dict([(pair.split(":")[1],int(pair.split(":")[0])) for pair in line[1:].strip("\n").split(",")])
				fromIDS = dict([(int(pair.split(":")[0]),(pair.split(":")[1])) for pair in line[1:].strip("\n").split(",")])

				MARGINAL=False
			elif MARGINAL and not header:
				TF, DNAse, eRNA 	= line.strip("\n").split("\t")
				M[TF]=[float(d) for d in DNAse.split(",")]+[float(d) for d in eRNA.split(",")]
			elif not header:
				line_array 	= line.strip("\n").split(",")
				L.append([(float(pair.split(":")[0]),float(pair.split(":")[1]) ) for pair in line_array])
			header = False
	A 	= np.zeros((len(L), len(L), 2))
	for i,l in enumerate(L):
		for j,(none, there) in enumerate(l):
			A[i,j,:] = (none, there)
	return A, M, toIDS, fromIDS
def choose2(n,k):
	if k > n:
		return 0
	r=1
	N=0
	for d in range(1,int(k+1)):
		r*=(n-N)
		r/=d
		N+=1
	return r
def choose(n,k):
	lgn1 = special.gammaln(n+1)
	lgk1 = special.gammaln(k+1)
	lgnk1 = special.gammaln(n-k+1)
	return lgn1 - (lgnk1 + lgk1)
class HG:
	def __init__(self,KS, n):
		self.KS 	= KS
		self.N 		= sum(self.KS)
		self.n 		= n
		assert sum(self.KS)
	def pmf(self, kis):

		TOP 	= np.sum([ choose(self.KS[i], kis[i]) for i,k in enumerate(kis) ])
		try:
			return math.exp(TOP - (choose(self.N,self.n)))
		except:
			return 1
	def expectation(self):
		return [(self.n*k) / self.N for k in self.KS]

def compute_pvalues(A, M, toIDS, fromIDS):

	NORM 	= 100.0

	W 		= 100000 
	eRNA 	= 30000  

	for i in range(A.shape[0]):
		for j in range(i+1,A.shape[0]):
			#first want to test if eRNA motif only
			tfa, tfb 	= fromIDS[i], fromIDS[j]
			tfa_eRNA 	= M[tfa][1] 
			tfb_eRNA 	= M[tfb][1] 
			
			tfa_tfb_eRNA= A[i,j,1]
			tfa_nottfba_eRNA= tfa_eRNA-tfa_tfb_eRNA
			p=100*float(tfa_nottfba_eRNA) /tfa_eRNA
			if p > 99:
				print tfa, "-", tfb
				print "All", tfa, "eRNA overlaps:", tfa_eRNA
				print "eRNA overlaps between", tfa, "-", tfb, ":", tfa_tfb_eRNA
				print "eRNA overlaps without", tfb, ":",tfa_nottfba_eRNA, str(100*float(tfa_nottfba_eRNA) /tfa_eRNA)[:3] + "%"
				print "----------------------------------"

				

				time.sleep(0.1)

if __name__ == "__main__":
	WRITE 		= False
	LOAD 		= True
	W 			= 14800
	eRNA 		= 2872
	eRNA_motif1_motif_2 = 6000
	motif_1_only 		= 190
	motif1_all 			= 1000
	motif2_all 			= 50000
	# hg 			= HG([motif1_all, motif2_all, W-(motif1_all+motif2_all)], eRNA)
	# print "------------------------------"
	# print "testing eRNA motif1 motif2 overlap"
	# print hg.pmf((eRNA_motif1_motif_2,eRNA_motif1_motif_2,eRNA_motif1_motif_2))
	# print hg.expectation()
	# print "------------------------------"
	print "testing eRNA motif1 only"
	hg 			= HG([motif1_all, W-motif1_all],eRNA)
	print hg.pmf((motif_1_only,eRNA-motif_1_only))
	print hg.expectation()
	print "------------------------------"
	hg 			= hypergeom(W, motif1_all,eRNA)
	print hg.pmf((motif_1_only,motif_1_only))
	print hg.mean()
	
	


	

	if LOAD:
		IN 		= "/Users/joazofeifa/Lab/TF_predictions/DNase_motif_eRNA_counts.tsv"
		A, M, toIDS, fromIDS=load_counts(IN)
		compute_pvalues(A, M, toIDS, fromIDS)
	if WRITE:
		OUT 		= "/Users/joazofeifa/Lab/TF_predictions/DNase_motif_eRNA_counts.tsv"
		eRNA_motif  = "/Users/joazofeifa/Lab/TF_predictions/assignments_refined/Allen2014_DMSO2_3-1_0.05"
		DNase_motif = "/Users/joazofeifa/Lab/TF_predictions/DNAse_motif_overlap.bed"
		B  			= ec.load(eRNA_motif,test=False)
		TFS 		= dict([(tf.split("_")[0],1) for chrom in B for b in B[chrom] for TF in b.TFS for tf in TF.split(",") ])
		D 			= load_DNAse(DNase_motif, FILTER=TFS, test=False)
		label_DNase(B,D, attr="eRNA")
		label_DNase(D,B, attr="DNAse")
		write_out(D, B, OUT, TFS)
	pass