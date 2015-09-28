import matplotlib.pyplot as plt
import time
import math 
import numpy as np
def refseq(FILE):
	G 	= {}
	with open(FILE) as FH:
		header=True
		for line in FH:
			if not header:
				start, stop 	= line.split("\t")[4:6]
				chrom 			= line.split("\t")[2]
				if chrom not in G:
					G[chrom] 	= list()
				G[chrom].append((int(start)-5000, int(stop)+10000 ))
			else:
				header=False 
	#merge
	A={}
	for chrom in G:
		G[chrom].sort()
		A[chrom] 	= list()	
		j,N 	=0,len(G[chrom])
		while j < N:
			o_st, o_sp 	=  G[chrom][j]
			while j < N and G[chrom][j][1] > o_st and G[chrom][j][0] < o_sp:
				o_st, o_sp 	= min(G[chrom][j][0],o_st ),max(G[chrom][j][1], o_sp )	
				j+=1
			A[chrom].append((o_st, o_sp))
	GAPS={}
	for chrom in A:
		i 	= 0
		while i < (len(A[chrom])-1):
			start 	= A[chrom][i][1]
			stop 	= A[chrom][i+1][0]
			if chrom not in GAPS:
				GAPS[chrom]=list()
			GAPS[chrom].append((start, stop))
			i+=1

	FHW=open("/Users/joazofeifa/Lab/genome_files/NO_ANNOTATION.bed", "w")
	for chrom  in GAPS:
		for start, stop in GAPS[chrom]:
			FHW.write(chrom + "\t" + str(start) + "\t" + str(stop) +"\n")
	FHW.close()
def getNONE(FILE):
	D 	= {}
	for line in open(FILE):
		chrom,start, stop 	= line.strip("\n").split("\t")
		if chrom not in D:
			D[chrom]=list()
		D[chrom].append((int(start), int(stop)))
	return D
def get_distributions(FILE,NONE):
	D 	= {}
	prevChrom 	= {}
	t 	= 0
	FH 	= open(FILE)
	for line in FH:

		chrom,start, stop, cov 	= line.strip("\n").split("\t")
		if chrom!=prevChrom:
			print chrom
			if chrom in NONE:
				j,N 	= 0,len(NONE[chrom])
			else:
				j,N 	= 0,0
			if t > 1:
				break
			t+=1
		while j < N and  NONE[chrom][j][1] < float(start):
			j+=1
		if j < N and  NONE[chrom][j][0] < float(stop):
			if chrom not in D:
				D[chrom] 	= {}
			if NONE[chrom][j][0] not in D[chrom]:
				D[chrom][ NONE[chrom][j][0]] 	= [NONE[chrom][j][1], 0]
			D[chrom][ NONE[chrom][j][0]][1]+=((float(stop)-float(start) )*float(cov))

		prevChrom=chrom
	FH.close()
	return D
def compare(HCT, DANK, PUC, LI):
	HCT 	= [float(HCT[chrom][start][1])/float( HCT[chrom][start][0]-start)  for chrom in HCT for start in HCT[chrom]  ]
	DANK 	= [float(DANK[chrom][start][1])/float(DANK[chrom][start][0]-start)  for chrom in DANK for start in DANK[chrom] ]
	PUC 	= [float(PUC[chrom][start][1])/float(PUC[chrom][start][0]-start)  for chrom in PUC for start in PUC[chrom] ]
	LI 		= [float(LI[chrom][start][1])/float(LI[chrom][start][0]-start)  for chrom in LI for start in LI[chrom] ]
	
	plt.hist([math.log(h) for h in HCT if h], alpha=0.5,label ="ALLEN, " + str(np.mean(HCT)) + "," + str(np.var(HCT)) , bins=50)
	plt.hist([math.log(h) for h in DANK if h], alpha=0.5, label="DANKO, "+ str(np.mean(DANK))+ "," + str(np.var(DANK)) ,bins=50)
	plt.hist([math.log(h) for h in PUC if h], alpha=0.5, label="PUC, "+ str(np.mean(PUC)) + "," + str(np.var(PUC)),bins=50)
	plt.hist([math.log(h) for h in LI if h], alpha=0.5, label="LI, "+ str(np.mean(LI))+ "," + str(np.var(LI)) ,bins=50)
	plt.legend()
	plt.show()

def to_name(chromosome):
	try:
		i  	= int(chrom[3:])
	except:
		i 	= 24
	return i

def make_gene_bed(FILE,OUT):
	G 	= {}
	FHW = open(OUT, "w")
	header 	= True
	with open(FILE) as FH :
		for line in FH:
			if not header:
				start, stop 	= line.split("\t")[4:6]
				start, stop 	= int(start), int(stop)
				chrom 			= line.split("\t")[2]
				if "random" not in chrom:
					name 			= line.split("\t")[1]
					if chrom not in G:
						G[chrom] 	= {}
					if start not in G[chrom]:
						G[chrom][start]={}
					if stop not in G[chrom][start]:
						G[chrom][start][stop] 	= name
			else:
				header = False
	chromosomes 	= [( to_name(chrom), chrom)  for chrom in G]
	chromosomes.sort()
	chromosomes 	= [CHR for i, CHR in chromosomes]
	LS 				= list()
	for chrom in chromosomes:
		A 	= [(start, stop,G[chrom][start][stop] ) for start in G[chrom] for stop in G[chrom][start] if 1000<(stop -start) < 200000 ]
		LS+=[ stop-start for start, stop, chrom in A]
		A.sort()
		for start, stop, name in A:
			FHW.write(chrom + "\t" + str(start) + "\t" + str(stop) + "\t" + name+ "\n" )

	

if __name__ == "__main__":

	GENE_BED 	= True
	NOISE_BED 	= False
	if GENE_BED:
		FILE = "/Users/joazofeifa/Lab/genome_files/RefSeqHG19.txt"
		OUT  =  "/Users/joazofeifa/Lab/genome_files/RefSeqHG19.bed"
		make_gene_bed(FILE,OUT)

	if NOISE_BED:

		# FILE = "/Users/joazofeifa/Lab/genome_files/RefSeqHG18.txt"
		# refseq(FILE)
		NONE 	= "/Users/joazofeifa/Lab/genome_files/hg19_no_anntation.bed"
		NONE2 	= "/Users/joazofeifa/Lab/genome_files/hg18_no_anntation.bed"

		HCT_F 	= "/Users/joazofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/DMSO2_3.pos.BedGraph"
		DANK_F 	= "/Users/joazofeifa/Lab/gro_seq_files/Danko_2014/GSE66031_ac16.unt.all_plus.bw.bedGraph.mod.bedGraph"
		PUC_F 	= "/Users/joazofeifa/Lab/gro_seq_files/Puc2015/bedgraph_files/siControl_1h_vehicle_hg18_forward.bedGraph"
		Li_F 	= "/Users/joazofeifa/Lab/gro_seq_files/Li2013/bedgraph_files/GSM1115996_Groseq-MCF7-E2-rep1_2.forward.bedGraph"
		NONE 	= getNONE(NONE)
		HCT 	= get_distributions(HCT_F, NONE)
		DANK 	= get_distributions(DANK_F,NONE)
		NONE 	= getNONE(NONE2)
		PUC 	= get_distributions(PUC_F,NONE)
		LI 		= get_distributions(Li_F, NONE)
		compare(HCT,DANK,PUC,LI)
	