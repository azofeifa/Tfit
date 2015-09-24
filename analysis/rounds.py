import sys,os
import time
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
class P:
	def __init__(self, ID, params):
		self.ID 	= ID
		self.mu,self.si,self.l,self.w,self.pi,self.N, self.fp, self.ll 	= [float(p) for p in params.split("_")]



def read_file(f):
	r 	= None
	G 	= {}
	FH 	= open(f)
	for line in FH:
		if line[:8] == "#-rounds":
			r 	= int(line.strip("\n").split(" : ")[1])
		elif line[0]!= "#":
			chrom,start,stop, params 	= line.strip("\n").split("\t")
			ID 							= chrom + ":" + start + "-" + stop
			G[ID] 						= P(ID, params)
	FH.close()
	return r,G

def read_through(DIR):
	R 			= {}
	for f in os.listdir(DIR):
		r,G 	= read_file(DIR+f)
		R[r] 	= G
	return R
def get_comps(R):
	master 	= R[max(R.keys())]
	M 		= max(R.keys())
	LLS 	= {}
	for r in R:
		LLS[r] 	= list()
		for ID in R[r]:
			if ID in master:
				if master[ID].N and R[r][ID].l>55:
				#	LLS[r].append((master[ID].ll-R[r][ID].ll) / R[r][ID].N)
				#	LLS[r].append(abs(master[ID].mu-R[r][ID].mu) /1)
				
					LLS[r].append(R[r][ID].si )

							
		else:
			print r
	return LLS
def mode(xs):
	counts,edges 	= np.histogram(xs, bins=200)
	edges 			= (edges[1:] + edges[:-1])/2.
	return max(zip(counts, edges))[1]
def plot(LLS):
	F 	= plt.figure(figsize=(15,10))
	ax 	= F.add_subplot(1,1,1)
	
	rs 	= LLS.keys()
	rs.sort()
	print rs
	ms 	= [mode(LLS[r]) for r in rs]
	print ms 
	ss 	= [np.std(LLS[r]) for r in rs]

	ax.scatter(rs, ms )
	ax.plot(rs, ms )
	
	ax.grid()
	plt.show()


if __name__ == "__main__":
	DIR 	= "/Users/joazofeifa/Lab/gro_seq_files/rounds_files/"
	R 		= read_through(DIR)
	LLS 	= get_comps(R)
	plot(LLS)
