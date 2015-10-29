
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
					if M not in self.TFS:
						self.TFS[M]=list()

				self.TFS[M].append(m_instant)
	def _construct_parameters(self, ps):
		self.mu,self.si,self.l, self.w, self.pi, self.N, self.fp, self.ll 	= [float(p) for p in ps.split("_")]

def load_distance_crude_file(FILE ,test=True):
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
				if chrom not in G:
					G[chrom]=list()
				G[chrom].append(b(chrom,start, stop, ps, motifs))
				t+=1
	return G

if __name__ == "__main__":
	DC 		= "/Users/joazofeifa/Lab/TF_predictions/distances_crude/"
	FILE  	= DC + "15_WT_t0_Saponaro2014-2_2000"
	FILE 	= "/Users/joazofeifa/Desktop/motif_distances.tsv"
	G 		= load_distance_crude_file(FILE)
