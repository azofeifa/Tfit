import node
def load_intervals(FILE):
	G 	= {}
	with open(FILE) as FH:
		for line in FH:
			chrom,start, stop, name 	= line.strip("\n").split("\t")
			if chrom not in G:
				G[chrom]=list()
			G[chrom].append((int(start),int(stop), name ))
	for chrom in G:
		G[chrom].sort()
		G[chrom] 	= node.tree(G[chrom])
	return G
def collect(FILES, G, OUT):
	A 	= {}
	for i,FILE in enumerate(FILES):
		with open(FILE) as FH:
			for line in FH:
				chrom,start, stop, cov 	= line.strip("\n").split("\t")
				if chrom in G:
					FINDS 	= G[chrom].searchInterval((int(start), int(stop)))
					for st, sp, name in FINDS:	
						if name not in A:
							A[name]=[0,0, st, sp, chrom]
						A[name][i]+=float(cov)
	FHW=open(OUT, "w")
	for name in A:
		pos, neg, st, sp, chrom 	= A[name]
		FHW.write(name+"\t" + chrom + "\t" + str(st) + "\t" + str(sp) + "\t" + str(pos)+","+str(neg)+"\n")




if __name__ == "__main__":

	HCT_forward 	= "/Users/joazofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/DMSO2_3.pos.BedGraph"
	HCT_reverse 	= "/Users/joazofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/DMSO2_3.neg.BedGraph"

	Andersson2014_forward 	= "/Users/joazofeifa/Lab/gro_seq_files/Andersson2014/bedgraph_files/SRR1596501.fastqbowtie2.sorted.pos.BedGraph"
	Andersson2014_reverse 	= "/Users/joazofeifa/Lab/gro_seq_files/Andersson2014/bedgraph_files/SRR1596501.fastqbowtie2.sorted.neg.BedGraph"

	Core2014_forward 	= "/Users/joazofeifa/Lab/gro_seq_files/Core2014/bedgraph_files/SRR1554312.fastqbowtie2.sorted.neg.BedGraph"
	Core2014_reverse 	= "/Users/joazofeifa/Lab/gro_seq_files/Core2014/bedgraph_files/SRR1554312.fastqbowtie2.sorted.pos.BedGraph"
	
	Luo2014_forward 	= "/Users/joazofeifa/Lab/gro_seq_files/Luo2014/bedgraph_files/SRR1015584.fastqbowtie2.sorted.pos.BedGraph"
	Luo2014_reverse 	= "/Users/joazofeifa/Lab/gro_seq_files/Luo2014/bedgraph_files/SRR1015584.fastqbowtie2.sorted.neg.BedGraph"


	RefSeq 		= "/Users/joazofeifa/Lab/genome_files/RefSeqHG19.bed"
	OUT 		= "/Users/joazofeifa/Lab/genome_files/Core2014Gene_Counts.tsv"
	G 			= load_intervals(RefSeq)
	collect((Core2014_forward, Core2014_reverse), G, OUT)




