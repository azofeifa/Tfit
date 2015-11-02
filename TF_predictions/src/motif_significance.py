import time
import math,re,sys
from scipy.special import erf
import matplotlib.pyplot as plt
import assign_TF as at
sys.path.append("/Users/azofeifa/Lab/")
sys.path.append("/Users/joazofeifa/Lab/")
from interval_searcher import intervals, node

class m:
	def __init__(self, parent, d, p_val, q_val, strand):
		self.parent = parent
		self.d , self.p_val, self.q_val, self.strand 	= d, p_val, q_val, strand
class b:
	def __init__(self, chrom,start, stop, ps, motifs):
		self.chrom,self.start , self.stop 	= chrom, int(start),int(stop)
		self.TFS 			  				= {}
		self._construct_parameters(ps)
		self._construct_motifs(motifs)
	def _construct_motifs(self, line):
		if line:
			motifs 	= line.split(",")
			for i,M in enumerate(motifs):
				motif, d 	= M.split("=")
				ds 			= d.split(";")

				for instance in ds:
					distance, p_val, q_val, strand 	= instance.split("_")
					m_instant 						= m(self,float(distance), float(p_val), float(q_val), strand)
					if motif not in self.TFS:
						self.TFS[motif]=list()

				self.TFS[motif].append(m_instant)
	def _construct_parameters(self, ps):
		self.mu,self.si,self.l, self.w, self.pi, self.N, self.fp, self.ll 	= [float(p) for p in ps.split("_")]

def load_distance_crude_file(FILE ,test=True, FILTER=None ):
	M,G  		= {}, {}
	interval 	= False
	t 			= 0
	with open(FILE ) as FH:
		for line in FH:

			if line[0]!= "-" and not interval:
				motif, total 	= line.strip("\n").split("\t")
				M[motif] 		= float(total)
			elif line[0]== "-" and not interval :
				interval=True
			elif interval:
				chrom,start, stop, ps, motifs 	= line.strip("\n").split("\t")
				start, stop 					= int(start), int(stop)
				if FILTER is None or not FILTER[chrom].searchInterval((start, stop)):
					if chrom not in G:
						print chrom,
						G[chrom]=list()
					G[chrom].append(b(chrom,start, stop, ps, motifs))
					if test and t > 5000:
						break
					t+=1
	print
	return G,M
def normal_pdf(x):
	pass
def normal_cdf(x, mu=0, si=1):
	return 0.5*(1+ erf((x-mu)/(si*math.sqrt(si))))
def get_significant_overlap(G, M):
	l 	= 0
	B 	= 0
	O 	= {}
	for chrom in G:
		MIN,MAX = 0,0	
		for b in G[chrom]:
			if MIN is None or b.start < MIN:
				MIN=b.start
			if MAX is None or b.stop > MAX:
				MAX=b.stop
			for TF in b.TFS:
				if TF not in O:
					O[TF]=0
				for d in b.TFS[TF]:
					if abs(d.d) < 10:
						O[TF]+=1
			B+=1
		l+=(MAX-MIN)
	w 	= 20
	p 	= w / (l*2.0)
	P 	= {}
	for motif in M:
		m 	= M[motif]
		n 	= m*B
		mu 	= n*p
		std = math.sqrt(mu*(1-p))
		k 	= O[motif]
		pv 	= normal_cdf(k, mu=mu, si=std)
		if pv > 0.999:
			print "ON", motif, k, mu, pv
		else:
			print "OFF", motif, k, mu, pv

		P[motif]=(k,mu,std, n,m,B, pv)
	return P
		
def load_RNA_seq(FILE):
	G 	= {}
	R 	= {}
	with open(FILE) as FH:
		for line in FH:
			gene, chrom,start, stop, cov 	= re.split("\s+", line.strip("\n"))[:5]
			G[gene] 	= (chrom,start, stop, float(cov) )
			if chrom not in R:
				R[chrom]=list()
			R[chrom].append((int(start), int(stop), gene))
	for chrom in R:
		R[chrom] 	= node.tree(R[chrom])
	return G,R
def seperate(R,P):
	N 	= 0.0
	T 	= 0.0
	ON 	= {}
	OFF = {}
	for motif in P:
		M 	= motif.split("_")[0]
		if M in R:
			k,mu,std,n,m,B ,pv 	= P[motif]
			#cov 	= R[M][3]
			cov 	= R[M]
			if pv > 0.9999:
				if M not in ON:
					ON[M]  = cov
				else:
					ON[M]=max(ON[M], cov )

			else:
				if M not in OFF:
					OFF[M] 	= cov
				else:
					OFF[M] 	= min(cov,OFF[M])
					
	plt.boxplot((OFF.values(), ON.values()) )
	plt.show()

def make_GRO_seq(R,FILES,OUT):
	FINDS 	= {}
	for F in FILES:
		t=0
		print F
		with open(F) as FH:
			for line in FH:

				chrom,start, stop, cov = line.strip("\n").split("\t")
				if chrom in R:
					F 	= R[chrom].searchInterval((int(start), int(stop)))
					if F:
						for f in F:
							M 	= f[2]
							if M not in FINDS:
								FINDS[M]=list()
							FINDS[M].append((int(start), float(cov)))
	FHW 	= open(OUT, "w")
	for M in FINDS:
		try:
			S 	= sum([cov for start, cov in FINDS[M]])
			S 	/= (FINDS[M][-1][0]-FINDS[M][0][0])
			FHW.write(M+"\t" + str(S) + "\n")
		except:
			print FINDS[M]
	return FINDS
def read_GRO(FILE):
	G={}
	with open(FILE) as FH:
		for line in FH:
			M,D 	= line.strip("\n").split("\t")
			G[M] 	= float(D)
	return G


if __name__ == "__main__":
	MAKE_GRO 	= False

	FILE 	= "/Users/joazofeifa/Lab/TF_predictions/distances_crude/Allen2014_DMSO2_3_motif_2000"
	RNA 	= "/Users/joazofeifa/Lab/TF_predictions/FIMOGenes.txt"
	OUT 	= "/Users/joazofeifa/Lab/TF_predictions/FIMO_GRO.txt"
	R,RT 	= load_RNA_seq(RNA)
	TSS_T 	= at.get_TSS_tree("/Users/joazofeifa/Lab/genome_files/TSS.bed")
		
	if MAKE_GRO:
		BD 		= "/Users/joazofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/"
		make_GRO_seq(RT,(BD+"DMSO2_3.pos.BedGraph", BD+"DMSO2_3.neg.BedGraph"),OUT)

	GRO 	= read_GRO(OUT)

	G,M 	= load_distance_crude_file(FILE, test=False, FILTER=None)
	P 		= get_significant_overlap(G,M)
	seperate(GRO,P)


