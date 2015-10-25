import os
import time
def load_motif(FILE, make_bed=""):
	G 	= {}
	header=True
	FHW 	= None
	if make_bed:
		FHW=open(make_bed, "w" )
	with open(FILE) as FH:
		for line in FH:
			if not header:
				line_array 	= line.split("\t")
				chrom,start, stop 	= line_array[1:4]
				if make_bed:
					FHW.write(chrom + "\t" + start + "\t" + stop+ "\n")
				q_val 				= line_array[6]
				if chrom not in G:
					G[chrom]=list()
				
				G[chrom].append((int(start), int(stop), float(q_val)))
			else:
				header=False
	for chrom in G:
		G[chrom].sort()
	return G
def make_query(FILE ):
	G 	= {}
	with open(FILE) as FH:
		for line in FH:
			
			chrom,start,stop 	= line.split("\t")[:3]
			if chrom not in G:
				G[chrom]=list()
			G[chrom].append((int(start), int(stop), list()))

	for chrom in G:
		G[chrom].sort()

	return G
def add_motifs(FILE,ID,G, make_bed=""):
	A 	= load_motif(FILE,make_bed=make_bed )
	for chrom in A:
		if chrom in G:
			a 	= A[chrom]
			b 	= G[chrom]
			j,N = 0,len(b)
			for a_st, a_sp, qval in a:
				while j < N and b[j][1] < a_st:
					j+=1
				if j < N and b[j][0] < a_sp:
					G[chrom][j][2].append((ID, qval))
	return G
def write_out(G, out=""):
	FHW 	= open(out, "w")
	for chrom in G:
		for start,stop, info in G[chrom]:
			INFO 	= ",".join([str(ID)+":" + str(qval) for ID, qval in info])
			if INFO:
				FHW.write(chrom+"\t"+str(start) + "\t"+ str(stop) + "\t"+ INFO  + "\n")
	FHW.close()



def read_in_directory(root, Q):
	A 	= {}
	for DIR in os.listdir(root):
		if os.path.exists(root+ DIR+ "/fimo.txt" ):
			Q 	= add_motifs(root+ DIR+ "/fimo.txt", DIR, Q,make_bed=root+DIR+".bed" )
	return Q

if __name__ == "__main__":
	root 	= "/Users/joazofeifa/Lab/ENCODE/HCT116/"
	query 	= "/Users/joazofeifa/Lab/ENCODE/HCT116/DNAse/peak_files/wgEncodeUwDnaseHct116PkRep1.narrowPeak"
	OUT 	= "/Users/joazofeifa/Lab/ENCODE/HCT116/DNase_motif_matches.bed"
	Q 		= make_query(query)
	Q 		= read_in_directory(root, Q)

	write_out(Q,out=OUT)

