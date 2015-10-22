import time, sys
sys.path.append("/Users/joazofeifa/Lab/interval_searcher/")
sys.path.append("/Users/joazofeifa/Lab/")
from interval_searcher import intervals, node
def load(FILE, bidir_bootstrap=False, test=False,FILTER=None):
	G,t 	= {},0
	with open(FILE) as FH:
		for line in FH:
			if test and t > 100:
				break
			t+=1
			if line[0]!= "#":
				chrom,start, stop 	= line.split("\t")[:3]
				start, stop 		= int(start), int(stop)
				if FILTER is None or (chrom in FILTER and not FILTER[chrom].searchInterval((int(start), int(stop)))):
					if chrom not in G:
						G[chrom] 		= list()
					if not bidir_bootstrap:
						G[chrom].append((int(start), int(stop)))
					else:
						ci 				= line.split("\t")[4].split("_")[0]
						G[chrom].append((int(start)-500 , int(stop)+500 , float(ci)))					

	for chrom in G:
		G[chrom].sort()
	return G
def get_intersection_chromosome(set_a, D, i):
	if i < len(D):
		set_b 	= set(D[i].keys())
		return get_intersection_chromosome(set_a & set_b, D, i+1)
	return set_a 
def compute_overlaps(*args):
	shared_chromosomes 	= get_intersection_chromosome(set(args[0].keys()), args, 1)
	NM,NB 				= 0.0,0.0
	MARK_ALL 			= sum([len(args[0][chrom]) for chrom in shared_chromosomes ])
	BIDIR_ALL 			= sum([len(args[1][chrom]) for chrom in shared_chromosomes ])
	for chrom in shared_chromosomes:
		LSTS 	= [a[chrom] for a in args]
		ST 		= intervals.comparison(LSTS)
		MI 		= ST.get_isolated(0)
		NI 		= ST.get_isolated(1)
		NM+=float(len(MI))
		NB+=float(len(NI))
		
	print  BIDIR_ALL-NB, 1 - (NB / BIDIR_ALL), MARK_ALL-NM, MARK_ALL, 1- (NM/MARK_ALL)  
def get_TSS_tree(FILE):
	G 	= {}
	with open(FILE) as FH:
		for line in FH:
			chrom, start, stop 	= line.strip("\n").split("\t")[:3]
			start, stop 		= int(start), int(stop)
			if chrom not in G:
				G[chrom] 	= list()
			G[chrom].append((start, stop))
	for chrom in G:
		G[chrom].sort()
		G[chrom]=node.tree(G[chrom])
	return G
if __name__ == "__main__":
	TSS_T 	= get_TSS_tree("/Users/joazofeifa/Lab/genome_files/TSS.bed")
	# me1 	= load("/Users/joazofeifa/Lab/ChIP/HCT116/HCT-116_H3K4me1_narrowPeak.bed",FILTER=TSS_T)
	# me3 	= load("/Users/joazofeifa/Lab/ChIP/HCT116/HCT116-DS16056.peaks.fdr0.01.hg19.bed",FILTER=TSS_T)
	# ac27 	= load("/Users/joazofeifa/Lab/ChIP/HCT116/HCT-116_H3K27Ac_narrowPeak.bed",FILTER=TSS_T)
	CTCF 	= load("/Users/joazofeifa/Lab/ChIP/HCT116/SL12242_Peaks.bed.broadPeak",FILTER=TSS_T)
	# Sp1 	= load("/Users/joazofeifa/Lab/ChIP/HCT116/SL12239_Peaks.bed.broadPeak",FILTER=TSS_T)
	Sin3A 	= load("/Users/joazofeifa/Lab/ChIP/HCT116/SL12244_Peaks.bed.broadPeak",FILTER=TSS_T)
	# JunD 	= load("/Users/joazofeifa/Lab/ChIP/HCT116/SL12233_Peaks.bed.broadPeak",FILTER=TSS_T)
	YY1 	= load("/Users/joazofeifa/Lab/ChIP/HCT116/ENCFF001UET.bed",FILTER=TSS_T)
	# ATF3 	= load("/Users/joazofeifa/Lab/ChIP/HCT116/ENCFF001UDL.bed",FILTER=TSS_T)
	# TEAD4 	= load("/Users/joazofeifa/Lab/ChIP/HCT116/ENCFF001UEP.bed",FILTER=TSS_T)
	# CEBPB 	= load("/Users/joazofeifa/Lab/ChIP/HCT116/ENCFF001UDP.bed",FILTER=TSS_T)
	#EGR1 	= load("/Users/joazofeifa/Lab/ChIP/HCT116/ENCFF001UDT.bed",FILTER=TSS_T)
	# ELF1 	= load("/Users/joazofeifa/Lab/ChIP/HCT116/ENCFF001UDV.bed",FILTER=TSS_T)
	# DNAse 	= load("/Users/joazofeifa/Lab/ChIP/HCT116/wgEncodeUwDnaseHct116PkRep1.narrowPeak.bed",FILTER=TSS_T)
	# DNAse2 	= load("/Users/joazofeifa/Lab/ChIP/HCT116/wgEncodeUwDnaseHct116PkRep2.narrowPeak.bed",FILTER=TSS_T)
	# TSS 	= load("/Users/joazofeifa/Lab/genome_files/TSS.bed")
	BIDIR 	= load("/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/Allen2014_boot_DMSO2_3-2_bootstrapped_bidirectional_hits_intervals.bed", bidir_bootstrap=True,FILTER=TSS_T)
	compute_overlaps(CTCF , BIDIR)



