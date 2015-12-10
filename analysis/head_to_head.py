import sys
sys.path.append("/Users/joazofeifa/Lab/interval_searcher/")
import node, time
import matplotlib.pyplot as plt
import numpy as np
def load_RefSeq_get_head_to_head(FILE):
	G 		= {}
	header 	= True
	with open(FILE) as FH:
		for line in FH:
			if not header:
				i, gene, chrom, strand, start, stop 	= line.split("\t")[:6]
				if chrom not in G:
					G[chrom] 	= list()
				G[chrom].append([int(start), int(stop), strand,None])
			else:
				header 	= False				
	return G
def read_bidireciotnal(FILE):
	G 	= {}
	i 	= 0
	with open(FILE) as FH:
		for line in FH:
			if "#" != line[0]:
				chrom,start, stop 	= line.split("\t")[:3]
				if chrom not in G:
					G[chrom]=list()
				G[chrom].append((int(start), int(stop), i))
				i+=1
	for chrom in G:
		G[chrom].sort()
		G[chrom] 	= node.tree(G[chrom])
	return G
def label(R, B):
	T 	= {}
	for chrom in R:
		if chrom in B:
			r,b 	= R[chrom], B[chrom]
			u 		= list()
			for i in range(len(r)):
				start, stop  =	r[i][0],r[i][1]
				if r[i][2]=="+":
					start, stop 	= start -2000, start+2000
				else:
					start, stop 	= stop -2000, stop+2000
					

				FINDS 		= b.searchInterval((start, stop))
				if FINDS and len(FINDS)==1:
					r[i][3] 	= FINDS[0][-1]
					u.append(r[i])
			T[chrom]=u
	return T
def draw(R):
	P 	= {}
	t 	= 0

	for chrom in R:
		R[chrom].sort()
		P[chrom] 	= {}
		r 	= R[chrom]
		N 	= len(r)
		for i,(start, stop, strand, ID) in enumerate(r):
			if strand == "-":
				j 	= i
				P[chrom][i]=list()
				while j < N:
					if r[j][2]=="+" and r[j][0] > stop and (r[j][0] - stop) < 6000:
						P[chrom][i].append(r[j])
					j+=1
		if t > 100:
			break
		t+=1
	x,y 	= list(),list()
	bins 	= [ [x,list() ] for x in np.linspace(0,6000, 25)]

	for chrom in P:
		for ID in P[chrom]:
			for head in P[chrom][ID]:
				d 	= head[0] - R[chrom][ID][1] 

				i 	= 1
				while i+1 < len(bins) :
					if bins[i-1][0] <= d <= bins[i+1][0]:
						bins[i][1].append( float(int(head[3]!=R[chrom][ID][3]) ))
					i+=1
	plt.scatter([x for x,y in bins if len(y)>1   ],[np.mean(y)+1  for x,y in bins  if len(y) > 1 ] )
	plt.plot([x for x,y in bins if len(y)>1   ],[np.mean(y)+1  for x,y in bins  if len(y) > 1 ] )
	plt.fill_between([x for x,y in bins if len(y)>1   ], 
		[np.mean(y)+1 - np.std(y)*0.5 for x,y in bins if len(y)>1   ],
		[np.mean(y)+1 + np.std(y)*0.5 for x,y in bins if len(y)>1   ], color="grey", alpha=0.5 )
	plt.grid()
	plt.xlabel("Distance from Head-to-Head Promoters (base pairs)")
	plt.ylabel("Number of Bidirectional Components")
	plt.xlim(0,6000)
	plt.show()










if __name__ == "__main__":
	FILE 	= "/Users/joazofeifa/Lab/genome_files/RefSeqHG19.txt"
	BIDIR 	= "/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/DMSO2_3-2_bidirectional_hits_intervals.bed"
	R 		= load_RefSeq_get_head_to_head(FILE) 
	B 		= read_bidireciotnal(BIDIR)
	R 		= label(R,B)
	draw(R)