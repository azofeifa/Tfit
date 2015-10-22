import os,sys
sys.path.append("/Users/azofeifa/Lab/")
from interval_searcher import intervals, node
import numpy as np
import math
import time
def make_query(FILE, pad=1000):
	G 	= {}
	IDS = {}
	i 	= 1
	with open(FILE) as FH:

		for line in FH:
			if line[0]!="#":
				chrom,start, stop 	= line.split("\t")[:3]
				if chrom not in G:
					G[chrom] 		= list()
				start, stop 		= int(start),int(stop)
				x 					= (stop + start) /2.
				st 					= start-((pad / 2.) - (x - start))
				sp 				 	= stop +((pad / 2.) - (stop - x)) 
				G[chrom].append((st, sp, i))
				IDS[i] 				= line.strip("\n")
				i+=1
	for chrom in G:
		G[chrom].sort()
		G[chrom] 	= node.tree(G[chrom])
	return G, IDS
def make_tree(FILE,A,G,TF, test=False):
	with open(FILE) as FH:
		header = True
		t 		= 0
		for line in FH:
			if not header:
				pattern_name, chrom, start, stop 	= line.split("\t")[:4]
				if test and t > 10000:
					break
				t+=1
				if chrom in G:
					start, stop 	= int(start), int(stop)
					FINDS 	= G[chrom].searchInterval(( start, stop))
					if len(FINDS):
						y 	= (start + stop) / 2.
						for st, sp, i in FINDS:
							x 	= (sp+st) /2.  
							if i not in A:
								A[i] 	= {}
							if TF not in A[i]:
								A[i][TF] = {}
							if pattern_name not in A[i][TF]:
								A[i][TF][pattern_name]=list()
							A[i][TF][pattern_name].append(x -y)

			else:
				header=False
	return A
def read_in_directory(root, Q):
	A 	= {}
	for DIR in os.listdir(root):
		if os.path.exists(root+ DIR+ "/fimo.txt" ):
			print DIR
			A 	= make_tree(root+ DIR+ "/fimo.txt", A, Q,DIR.split("_")[0] )
	return A
def write_out(A,IDS, OUT):
	FHW= open(out+ "motif_distances.tsv", "w")
	for i in range(1, max(IDS.keys() ) +1):
		ID 	= IDS[i]
		D 	= ID + "\t"
		if i in A:
			for TF in A[i]:
				D+=TF+":"
				for j,pattern_name in enumerate(A[i][TF]):
					if j!=0:
						D+="_" + pattern_name+ "="
					else:
						D+=pattern_name+ "="

					for d in A[i][TF][pattern_name]:
						D+=str(d)+","
					D.strip(",")
		FHW.write(D+ "\n")
	FHW.close()




if __name__ == "__main__":
<<<<<<< HEAD
	root 	= "/Users/azofeifa/FIMO_OUT/"
	query 	= "/Users/azofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/Allen2014_DMSO2_3-1_bidirectional_hits_intervals.bed"
	out 	= "/Users/azofeifa/"
=======
	if len(sys.argv)<2:
		root 	= "/Users/joazofeifa/Lab/ENCODE/HCT116/"
		query 	= "/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/Allen2014_DMSO2_3-4_bidirectional_hits_intervals.bed"
		out 	= "/Users/joazofeifa/Desktop/"
	else:
		root 	= sys.argv[1]
		query 	= sys.argv[2]
		out 	= sys.argv[3]
>>>>>>> 7925850a5cf14ff2cc17b349714cec9bf55eecc2
	Q,IDS 	= make_query(query)
	A 		= read_in_directory(root, Q)
	write_out(A,IDS, out)
