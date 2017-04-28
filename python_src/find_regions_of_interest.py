import time
import math
def load(FILE):
	G 	= {}
	with open(FILE) as FH:
		t 	= 0
		for line in FH:
			chrom,start,stop, cov 	= line.strip("\n").split("\t")
			start, stop 			= float(start), float(stop)
			x 						= ((stop)+(start))/2.
			cov 					= float(cov)
			if chrom not in G:
				G[chrom] =list()
				if t > 1:
					break
				t+=1
			G[chrom].append((x,cov*(stop-start+1) ))
	return G

def window_search(G, window_size=1000, thresh=1.0):
	A 	= {}
	for chrom in G:
		D 	= G[chrom]
		N 	= len(D)
		start, stop = None, None
		j,k 		= 0,0
		S 	= 0
		for i,(x,cov) in enumerate(D):
			while  j < N  and (D[j][0]-x) < -window_size:
				S-=D[j][1]
				j+=1
			
						
			while 0<=k < N and (D[k][0] - x  ) < window_size:
				S+=D[k][1]
				k+=1
			density 	=  S/ window_size
			if (density >=thresh and start is None):
				start 	= x
			elif (density >=thresh):
				stop 	= x
			elif stop is not None:
				if chrom not in A:
					A[chrom] 	= list()
				A[chrom].append((start, stop)) 
				stop,start 	= None, None
			else:
				start,stop 	= None, None
	return A

def overlaps(A,B,out):
	FHW = open(out, "w")
	for chrom in A:
		if chrom in B:
			a,b=A[chrom],B[chrom]
			j,N 	= 0,len(b)
			for a_st, a_sp in a:
				while j < N and b[j][1] < a_st:
					j+=1
				if j < N and b[j][0] < a_sp:
					o_st, o_sp 	= max(b[j][0],a_st),min(b[j][1],a_sp)
					FHW.write(chrom+"\t" + str(int(o_st)) + "\t" + str(int(o_sp)) + "\n" )
	FHW.close()


if __name__=="__main__":
	forward="/Users/joazofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/DMSO2_3.pos.BedGraph"
	reverse="/Users/joazofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/DMSO2_3.neg.BedGraph"
	OUT  	= "/Users/joazofeifa/Lab/gro_seq_files/HCT116/window.bed"
	F 		= load(forward)
	print "loaded forward"
	R 		= load(reverse)
	print "loaded reverse"
	F 		= window_search(F, window_size=1000,thresh=1.2 )
	print "window_searched forward"
	R 		= window_search(R, window_size=1000, thresh=1.2)
	print "window_searched reverse"
	overlaps(F,R,OUT)


