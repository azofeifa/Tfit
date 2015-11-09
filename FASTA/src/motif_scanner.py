import sys,re
import numpy as np
#ACGT
import time
import math as m
def read_motif_file(FILE):
	L 	= list()

	with open(FILE) as FH:
		header=True
		for line in FH:
			if not header:
				freqs 	= [float(x) for x in line.strip("\n").split("\t")]	
				L.append(freqs)
			else:
				header 	= False
	return np.array(L)
COMP 	= {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N"}
POS 	= {"A":0, "T":3, "G":2, "C":1 }
def search_seq(M, seq, TOP=100):
	cseq 	= "".join([COMP[c] for c in seq])
	N 		= M.shape[0]
	scores 	= list()
	for i in range(len(seq)-N):
		ll 	= sum([m.log(M[j,POS[c]] ) if c in POS else -np.inf for j,c in enumerate(seq[i:i+N])  ])
		scores.append((ll,i, "+"))
	for i in range(len(seq)-N  ):
		ll 	= sum([m.log(M[j,POS[c]] ) if c in POS else -np.inf for j,c in enumerate(cseq[i:i+N])  ])
		scores.append((ll,i, "-"))

	scores.sort()
	scores=scores[::-1]
	return scores[:TOP]


def search_fasta(M, FILE, OUT, TOP=100):
	seq=""
	FHW=open(OUT,"w")
	with open(FILE) as FH:
		for line in FH:
			if ">" ==line[0]:
				if seq:
					scores 	= search_seq(M,seq, TOP=TOP)
					FINDS 	= ",".join([str(ll)+":"+str(i) +":"+s for ll,i,s in scores])
					FHW.write(header+"\t" + FINDS+"\n")

				header 	= line[1:].strip("\n")
				seq 	= ""
			else:
				seq+=line.strip("\n")



if __name__ == "__main__":
	motif_file 	= "/Users/joazofeifa/opt/HOMMER/motifs/sp1.motif"
	fasta_file 	= "/Users/joazofeifa/Lab/TF_predictions/HOMMER_OUT_FILES/Allen2014_rawChIP_fasta.fasta"
	OUT_file 	= "/Users/joazofeifa/Lab/TF_predictions/HOMMER_OUT_FILES/sp1_sites.bed"
	M 			= read_motif_file(motif_file)
	search_fasta(M, fasta_file,OUT_file, TOP=500)
