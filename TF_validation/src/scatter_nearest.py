import nearest_neighbor_ChIP_comparison as nn
import math
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import numpy as np
import time, sys
sys.path.append("/Users/joazofeifa/Lab/EMG/TF_predictions/src/")
sys.path.append("/Users/joazofeifa/Lab/interval_searcher/")
import node

def load_genes(FILE):
	G 	= {}
	with open(FILE) as FH :
		for line in FH:
			name,chrom,start, stop, cov 	= line.strip("\n").split("\t")
			cov 							= sum([abs(float(x)) for x in cov.split(",")])
			G[name]=chrom,int(start), int(stop), math.log(cov+1 / (float(stop) - float(start)+1),10)
	return G



def scatter(A, B,genes, CS, labels, WINDOW=pow(10,2) ):
	pass



if __name__ == "__main__":
	#============================================================================================
	#TF_BIDIR Input Files
	HCT116 			= "/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_ChIP_comparisons/"
	K562 			= "/Users/joazofeifa/Lab/gro_seq_files/Core2014/EMG_out_ChIP_comparisons/"
	HeLa 			= "/Users/joazofeifa/Lab/gro_seq_files/Andersson2014/EMG_out_ChIP_comparisons/"
	A549 			= "/Users/joazofeifa/Lab/gro_seq_files/Luo2014/EMG_out_ChIP_comparisons/"
	DIRS 			= (  HeLa, A549, K562, HCT116 )
	labels 			= (  "HeLa", "A549","K562", "HCT116" )
	
	CS 				= nn.load_all(DIRS, labels=labels, test=False )
	
	genes_HCT116 	= "/Users/joazofeifa/Lab/genome_files/DMSO2_3_Gene_Counts.tsv"
	genes_K562 		= "//Users/joazofeifa/Lab/genome_files/Core2014Gene_Counts.tsv"
	genes_HeLa 		= "/Users/joazofeifa/Lab/genome_files/Andersson2014Gene_Counts.tsv"
	genes_A549 		= "/Users/joazofeifa/Lab/genome_files/Luo2014Gene_Counts.tsv"

	genes_DIR 	 	= (genes_HCT116,genes_K562,genes_A549, genes_HeLa )
	#genes_DIR 	 	= (genes_HCT116, genes_A549, genes_HeLa )
	labels 			= ("HCT116", "K562", "HeLa", "A549" )

	genes 			= nn.load_gene_count_files(genes_DIR, labels )

	scatter(HeLa, A549, "HeLa", "A549", CS)










