import matplotlib.pyplot as plt
import time
import numpy as np
import math
from scipy.stats import gaussian_kde
import matplotlib.ticker
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Helvetica']
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats



def load(FILE):
	G 	= {}
	with open(FILE) as FH:
		for line in FH:
			if line[0]!="#":
				chrom,start, stop, info 	= line.strip("\n").split("\t")
				start, stop 				= int(start),int(stop)
				info 						= [float(x) for x in info.split("|")[1].split(",")]
				if chrom not in G:
					G[chrom]=list()
				if 200> info[1] > 20 and 200>info[4] > 50 and 1000>info[2] > 50:
					G[chrom].append((start, stop, info))
	for chrom in G:
		G[chrom].sort()
	return G
def match_up(A,B):
	M 	= {}
	for chrom in A:
		a,b 	= A[chrom], B[chrom]
		j,N 	= 0,len(b)
		M[chrom]=list()
		for start, stop, info in a:
			while j <N and b[j][1] < start:
				j+=1
			if j < N and b[j][0] < stop:
				o_st, o_sp 	= max(b[j][0], start), min(b[j][1], stop)
				d 			= float(o_sp - o_st)
				if d / (stop - start) > 0.9:
					M[chrom].append((info, b[j][2]))
	return M
def draw(M,t, LOG=False):
	left, width = 0.1, 0.65
	bottom, height = 0.1, 0.65
	bottom_h = left_h = left+width+0.02
	rect_scatter = [left, bottom, width, height]
	rect_histx = [left, bottom_h, width, 0.2]
	rect_histy = [left_h, bottom, 0.2, height]
	F 			= plt.figure(figsize=(15,10))
	
	axScatter = plt.axes(rect_scatter)
	axHistx = plt.axes(rect_histx)
	axHisty = plt.axes(rect_histy)

	x,y 	= [i[t] if not LOG else math.log(i[t],10) for chrom in M for i,j in M[chrom]],[j[t] if not LOG else math.log(j[t],10) for chrom in M for i,j in M[chrom]]
	


	xy = np.vstack([x,y])
	
	z = gaussian_kde(xy)(xy)
	
	axScatter.scatter(x, y, c=z, s=14, edgecolor='' )
	axScatter.grid()
	axScatter.set_xlim(min(x), max(x))
	axScatter.set_ylim(min(x), max(x))
	counts,edges 	= np.histogram(x,bins=200,normed=1)
	edges 			= (edges[1:] + edges[:-1])/2.
	
	axHistx.bar(edges, counts, edgecolor="",width=0.5*(edges[-1]-edges[0])/200)
	counts,edges 	= np.histogram(y,bins=200,normed=1)
	edges 			= (edges[1:] + edges[:-1])/2.
	axHisty.barh(edges, counts, edgecolor="",height=0.5*(edges[-1]-edges[0])/200)
	axHistx.set_xlim( axScatter.get_xlim() )
	axHisty.set_ylim( axScatter.get_ylim() )
	
	axHistx.set_xticks([])
	axHisty.set_yticks([])
	plt.show()
def compute_correlations(M):
	for t in range(6):
		for u in range(t, 6):
			if t ==1 or t==2 or t==4:
				LOG1 	= True
			else:
				LOG1 	= False
			if u ==1 or u==2 or u==4:
				LOG2 	= True
			else:
				LOG2 	= False
			x,y 	= [i[t] if not LOG1 else math.log(i[t],10) for chrom in M for i,j in M[chrom]],[j[u] if not LOG2 else math.log(j[u],10) for chrom in M for i,j in M[chrom]]
			print t,"->",u,np.corrcoef(x,y=y)[0,1]



if __name__ == "__main__":
	rep1 	= "/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/Allen2014_DMSO2_3-2_divergent_classifications.bed"
	rep2 	= "/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/Allen2014_Nutlin2_3-1_divergent_classifications.bed"
	rep1 	= load(rep1)
	rep2 	= load(rep2)
	M 		= match_up(rep1, rep2)
	compute_correlations(M)
