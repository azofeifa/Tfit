import ROC_from_k_model_fits as roc
import node,time,math
import matplotlib.pyplot as plt
def load_DNase(FILE):
	G 	= {}
	with open(FILE) as FH:
		for line in FH:
			chrom,start, stop  	= line.strip("\n").split("\t")[:3]
			if chrom not in G:
				G[chrom]=list()
			G[chrom].append((int(start), int(stop)))
	for chrom in G:
		G[chrom].sort()
	return G

def load_k_model_fits(FILE):
	S 	= None
	G 	= {}
	with open(FILE) as FH:
		for line in FH:
			if line[0]==">":
				if S is not None:
					if S.chrom not in G:
						G[S.chrom]=list()
					G[S.chrom].append((S.start, S))
				S 	= roc.segment(line)
			elif S is not None:
				S.add_fit(line)
	if S is not None and S.check():
		G[S.chrom].append((S.start, S))
	for chrom in G:
		G[chrom].sort()
		G[chrom] 	= [y for x,y in G[chrom]]
	return G
def load_gene_counts(FILE):
	G 	= {}
	with open(FILE) as FH:
		for line in FH:
			name, chrom,start, stop, pos_neg 	= line.strip("\n").split('\t')
			if chrom not in G:
				G[chrom]=list()
			G[chrom].append((int(start), int(stop), (name, pos_neg) ))
	for chrom in G:
		G[chrom].sort()
		G[chrom]=node.tree(G[chrom])
	return G
def get_DNase_overlaps(A,B):
	t,T 	= 0.0,0.0
	for chrom in A:
		if chrom in B:
			a,b 	= A[chrom],B[chrom]
			j,N 	= 0,len(b)
			for a_s in a:
				while j < N and b[j][1] < a_s.start:
					j+=1
				if j < N and b[j][0] < a_s.stop:
					a_s.DNase 	= True
					t+=1
				T+=1
			A[chrom]=a
	return A,B
def get_nearest_gene(A,R,W=50000):
	for chrom in A:
		a 	= A[chrom]
		if chrom in R:
			r 	= R[chrom]
			for a_s in a:
				#FINDS 	= r.searchInterval((a_s.start, a_s.stop))
				FINDS 	= r.searchInterval((a_s.start-W, a_s.stop+W))
				if FINDS:
					center 		= (a_s.start + a_s.stop) /2.
					nearest 	= min([(abs( center - ((  start + stop) / 2.)), FINDS[i])  for i,(start, stop, info) in enumerate(FINDS) ])
					a_s.nearest_gene 	= nearest
	return A
def write_out(A, OUT):
	FHW 	= open(OUT,"w")
	for chrom in A:
		a 	= A[chrom]
		for a_s in a:
			nearest_gene 	= "NA"
			if a_s.nearest_gene is not None:
				d,nearest_gene 	=  a_s.nearest_gene
				nearest_gene= str(nearest_gene[0]) +"-"+str(nearest_gene[1]) +","+ str(d) + ","+ nearest_gene[2][0] +  "," + nearest_gene[2][1]
			FHW.write(chrom+"\t" + str(a_s.start)+"\t"+str(a_s.stop) +"\t" +str(a_s.K[0])+"," + str(a_s.K[1]) +"\t"+ str(a_s.DNase)+"\t" + nearest_gene + "\t" + 
				str(a_s.fN)+","+ str(a_s.rN)+ "\n"  )	
class segment:
	def __init__(self,chrom,start, stop, lls, DNase, nearest_gene, gro_N, penality=700):
		self.chrom,self.start, self.stop 	= chrom,float(start), float(stop)
		self.lln,self.llm 					= [float(x) for x in lls.split(",")]
		self.DNase 							= (DNase=="True")
		#2984839-3032270,52616.0,NM_004042,87.0,122.0
		gene_array 	= nearest_gene.split(",")
		self.gene_dist,self.gene_name 		= float(gene_array[1] ),gene_array[2]
		self.gene_pos,self.gene_neg  		= float(gene_array[3]), float(gene_array[4])
		self.gro_pos_N,self.gro_neg_N 		= [float(x) for x in gro_N.split(",")]
		self.gro_N 							= self.gro_pos_N + self.gro_neg_N
		self.gro_N 							/= (float(gene_array[0].split("-")[1]) - float(gene_array[0].split("-")[0]))
		null, model 						= self.lln*(-2) + math.log(self.gro_N), self.llm*(-2) + penality*math.log(self.gro_N)
		self.bidir 							= model < null

def load(FILE):
	G 	= list()
	with open(FILE) as FH:
		for line in FH:
			chrom,start,stop, lls,DNase, nearest_gene, gro_N 	= line.strip("\n").split("\t")
			if nearest_gene!= "NA":
				S 	= segment(chrom,start, stop, lls, DNase, nearest_gene, gro_N)
				G.append(S)
	return G

def boxplot(G):
	DNAse_no_bidir 	= [math.log(s.gene_pos + s.gene_neg,10) for s in G if s.DNase and not s.bidir or not s.DNase    ]
	DNAse_bidir 	= [math.log(s.gene_pos + s.gene_neg,10) for s in G if  s.DNase and  s.bidir]
	plt.hist(DNAse_no_bidir,alpha=0.5,label="No Bidir: " + str(len(DNAse_no_bidir)),bins=50, normed=1)
	plt.hist(DNAse_bidir,alpha=0.5,label="Bidir: " + str(len(DNAse_bidir)),bins=50, normed=1)
	plt.legend()
	plt.show()



if __name__ == "__main__":
	WRITE 			= True
	if WRITE:
		ELF1_file 		= "/Users/joazofeifa/Lab/nearest_gene_expression/ChIP_bidir_files/ELF1_ChIP-2_K_models_MLE.tsv"
		ATF_file 		= "/Users/joazofeifa/Lab/nearest_gene_expression/ChIP_bidir_files/ATF_ChIP-2_K_models_MLE.tsv"
		CBX3_file 		= "/Users/joazofeifa/Lab/nearest_gene_expression/ChIP_bidir_files/CBX3_ChIP-2_K_models_MLE.tsv"
		Rad21_file 		= "/Users/joazofeifa/Lab/nearest_gene_expression/ChIP_bidir_files/Rad21_ChIP-2_K_models_MLE.tsv"
	
		DNase_file 		= "/Users/joazofeifa/Lab/nearest_gene_expression/wgEncodeUwDnaseHct116PkRep1.narrowPeak.bed"
		gene_cts_file 	= "/Users/joazofeifa/Lab/genome_files/DMSO2_3_Gene_Counts.tsv"
	
		ELF1_OUT 		= "/Users/joazofeifa/Lab/nearest_gene_expression/ELF1_DNAse_nearest_gene.tsv"
		ATF_OUT 		="/Users/joazofeifa/Lab/nearest_gene_expression/ATF_DNAse_nearest_gene.tsv"
		CBX3_OUT 		="/Users/joazofeifa/Lab/nearest_gene_expression/CBX3_DNAse_nearest_gene.tsv"
		Rad21_OUT 		="/Users/joazofeifa/Lab/nearest_gene_expression/Rad21_DNAse_nearest_gene.tsv"

		G 				= load_k_model_fits(Rad21_file)
		DNase 			= load_DNase(DNase_file)
		gene_cts 		= load_gene_counts(gene_cts_file)
		G, DNase 		= get_DNase_overlaps(G, DNase)
		G 				= get_nearest_gene(G, gene_cts)
	
		write_out(G, Rad21_OUT)
	ELF1_IN 			= "/Users/joazofeifa/Lab/nearest_gene_expression/ELF1_DNAse_nearest_gene.tsv"
	ATF_IN 				= "/Users/joazofeifa/Lab/nearest_gene_expression/ATF_DNAse_nearest_gene.tsv"
	CBX3_IN 			= "/Users/joazofeifa/Lab/nearest_gene_expression/CBX3_DNAse_nearest_gene.tsv"
	Rad21_IN 			= "/Users/joazofeifa/Lab/nearest_gene_expression/Rad21_DNAse_nearest_gene.tsv"

	G 					= load(ATF_IN)
	OUT 				= open("/Users/joazofeifa/something.bed", "w")
	for s in G:
		OUT.write(s.chrom+"\t" + str(int(s.start)) + "\t" + str(int(s.stop)) + "\t" + str(s.bidir)+"\n" )




	
	boxplot(G)
