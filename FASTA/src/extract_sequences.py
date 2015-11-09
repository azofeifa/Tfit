import time
def load_bidir(FILE, window=0):
	G 	= {}
	with open(FILE) as FH:
		for line in FH:
			line_array 	= line.strip("\n").split("\t")
			if len(line_array) > 3:
				chrom,start, stop, params, motfs, RAW 	= line_array
				if chrom not in G:
					G[chrom]=list()
				if window:
					center 	= float(int(start) + int(stop) )/2.
					start, stop = str(int(center-window)), str(int(center+window))
				G[chrom].append([int(start), int(stop), chrom+"|"+start+"-"+stop+"|" + params+ "|" + RAW+ "\n", "" ])


	for chrom in G:
		G[chrom].sort()
	return G
def insert_sequence(FASTA_FILE,G, OUT=""):
	chrom 	= None
	with open(FASTA_FILE) as FH:
		for line in FH:
			if ">" == line[0]:
				chrom=line[1:].strip("\n")
				pos 	= 0

				if chrom in G:
					j,N = 0,len(G[chrom])
					print chrom,
				else:
					j,N = 0,0
			elif j < N:
				line 	= line.strip("\n")
				L 		= len(line)
				pos 	+=L
				while j < N and G[chrom][j][1] < pos:
					j+=1
				if j < N and G[chrom][j][0] < pos:
					for i,c in enumerate(line):
						if G[chrom][j][0]<= pos+i <=G[chrom][j][1]:
							G[chrom][j][3]+=c
	FHW=open(OUT, "w")
	for chrom in G:
		for start,stop, header, seq in G[chrom]: 
			FHW.write(">"+header)
			FHW.write(seq.upper()+"\n")
	FHW.close()
	return G
if __name__ == "__main__":
	BIDIR_CHIP 	= "/Users/joazofeifa/Lab/TF_predictions/distances_crude/Allen2014_rawCHIP_motif_distances.tsv"
	FASTA_FILE 	= "/Users/joazofeifa/opt/genome_files/hg19.fasta"
	G 			= load_bidir(BIDIR_CHIP, window=1500)
	G 			= insert_sequence(FASTA_FILE, G, OUT="/Users/joazofeifa/Lab/TF_predictions/HOMMER_OUT_FILES/Allen2014_rawChIP_fasta.fasta")


