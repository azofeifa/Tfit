import os,sys
sys.path.append("/Users/azofeifa/Lab/")
sys.path.append("/Users/joazofeifa/Lab/")
from interval_searcher import intervals, node
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import numpy as np
from scipy import stats
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import matplotlib.cm as cm


class ChIP:
	def __init__(self, start, stop, parameters):
		self.center 			= (float(stop) + float(start)) /2.
		self.start,self.stop 	= self.center-1000, self.center+1000
		self.TSS_D 			= None
		self.TSS_strand 	= None
		self.match 			= None
		self._construct(parameters)
	def _construct(self, line):
		self.mu,self.si,self.l, self.pi,self.w , self.fp, self.wl,self.wr, self.ll, self.N 	= [float(x) for x in line.split(",")]
		self.density 	= (self.N+1) / ((self.stop -self.start)+1)
		self.fp*=100
class GRO:
	def __init__(self, start, stop, parameters):
		self.center 			= (float(stop) + float(start)) /2.
		self.start,self.stop 	= self.center-1000, self.center+1000
		self.TSS_D 			= None
		self.TSS_strand 	= None
		self.match 		= None
		self._construct(parameters)
	def _construct(self, line):
		self.mu,self.si,self.l, self.w,self.pi,self.N,self.fp,self.ll 	= [float(x) for x in line.split("_")]
		self.density 	= (self.N+1) / ((self.stop -self.start)+1)

def load_file(FILE, CHIP=True):
	G 	= {}
	with open(FILE) as FH:
		for line in FH:
			if "#" != line[0]:
				chrom,start, stop, parameters 	= line.strip("\n").split("\t")
				if chrom not in G:
					G[chrom]=list()
				if CHIP:
					D 	= ChIP(start, stop, parameters)
				else:
					D 	= GRO(start, stop, parameters)
				G[chrom].append((int(start), D))
	for chrom in G:
		G[chrom].sort()
		G[chrom] 	= [y for x,y in G[chrom]]

	return G
def load_TSS(FILE):
	G 	= {}
	with open(FILE) as FH:
		header=True
		for line in FH:
			if not header:

				chrom, strand, start,stop =line.strip("\n").split("\t")[2:6]
				if strand== "-":
					start = stop
				if chrom not in G:
					G[chrom]=list()
				G[chrom].append((int(start), strand))
			else:
				header=False
	for chrom in G:
		G[chrom].sort()
	return G
def match_up(A,B):
	for chrom in A:
		a 	= A[chrom]
		b 	= B[chrom]
		j,N = 0,len(b)
		for a_i in a:
			while j < N and b[j].stop < a_i.start:
				j+=1
			k 	= j
			MIN 	= None
			MATCH 	= None
			while k < N and b[k].start < a_i.stop:
				D 	= abs(a_i.center - b[k].center)
				if D < 1000 and (MIN is None or abs(a_i.center - b[k].center) < MIN) :
					MATCH 		= b[k]
					MIN 		= abs(a_i.center - b[k].center)
				k+=1
			if MIN is not None:
				a_i.match 	= MATCH
				MATCH.match 	= a_i

def label_TSS(A, T):
	for chrom in A:
		b 		= A[chrom]
		a 		= T[chrom]
		j,N 	= 0,len(b)
		for start, strand in a:
			while j < N and b[j].stop < start:
				j+=1
			if j < N and b[j].start < start:
				if b[j].TSS_D:
					if abs(b[j].TSS_D) > abs(start - b[j].center ):
						b[j].TSS_D 	= b[j].center - start 
							
						b[j].TSS_strand 	= strand
				else:
					b[j].TSS_D 			= b[j].center - start 
						
					b[j].TSS_strand 	= strand
				
def correlation_TSS(ChIP, GRO):
	F 	= plt.figure(figsize=(15,10))
	ax 	= F.add_subplot(1,1,1)
	xy 	= [ (c.TSS_D, c.match.TSS_D) for chrom in ChIP for c in ChIP[chrom] if c.match and c.TSS_D and c.match.TSS_D and abs(c.TSS_D) < 1000 and abs(c.match.TSS_D) < 1000]
	
	x,y = [x for x,y in xy],[y for x,y in xy]
 	ro, p 	= stats.pearsonr(x,y)
 	
	xy = np.vstack([x,y])
	z = gaussian_kde(xy)(xy)
	
	ax.scatter(x, y, c=z, s=14, edgecolor='', label="pearsons: " + str(ro)[:5] +"\n GRO mean: " + str(np.mean(y))[:5] + "\n ChIP mean: " + str(np.mean(x))[:5] )

	ax.set_xlim(-1100,1110)
	ax.set_ylim(-1100,1110)
	ax.legend(loc=(0.1,0.85))	
	ax.grid()
	ax.set_xlabel("ChIP Pol II model fit " + r'$\mu$')
	ax.set_ylabel("GRO-seq model fit " + r'$\mu$')
	
	plt.show()
def make_heatmap_histogram(ax, d, tilt=False, lim_min=0, lim_max=0):
	counts,edges 	= np.histogram(d,bins=50)
	edges 			= (edges[1:] + edges[:-1]) /2.

	norm    = mpl.colors.Normalize(vmin=min(counts), vmax=max(counts))
	cmap    = cm.jet
	m       = cm.ScalarMappable(norm=norm, cmap=cmap)
	colors  = [m.to_rgba(c) for c in counts] 


	# In[46]:
	if tilt:
		ax.barh(edges,np.ones((len(edges),)), color=colors, edgecolor=colors )
	#	ax.set_ylim(lim_min, lim_max)
	else:
		ax.bar(edges,np.ones((len(edges),)), color=colors, width=(edges[-1]-edges[0])/len(edges) , edgecolor=colors )
	#	ax.set_xlim(lim_min, lim_max)

	pass

def correlations_parameters(ChIP, GRO, p="", LOG=True):
	
	# definitions for the axes
	left, width = 0.1, 0.75
	bottom, height = 0.1, 0.75
	bottom_h = left_h = left+width+0.02
	rect_scatter = [left, bottom, width, height]
	rect_histx = [left, bottom_h, width, 0.1]
	rect_histy = [left_h, bottom, 0.1, height]
	F 			= plt.figure(figsize=(15,10))
	
	axScatter = plt.axes(rect_scatter)
	axHistx = plt.axes(rect_histx)
	axHisty = plt.axes(rect_histy)
	axHistx.set_xticklabels([])
	axHistx.set_yticklabels([])
	axHisty.set_xticklabels([])
	axHisty.set_yticklabels([])
	
	if LOG:
		xy 	= [ (math.log(getattr(c, p),10), math.log(getattr(c.match, p), 10)) for chrom in ChIP for c in ChIP[chrom] if c.match and c.TSS_D is not None ]
	else:
		xy 	= [ (getattr(c, p), getattr(c.match, p)) for chrom in ChIP for c in ChIP[chrom] if c.match and c.TSS_D is not None ]
		
	x,y = [x for x,y in xy  ],[y for x,y in xy  ]
 	
	ro, p2 	= stats.pearsonr(x,y)
 	
	xy = np.vstack([x,y])
	z = gaussian_kde(xy)(xy)

	make_heatmap_histogram(axHistx, x, lim_min=np.min(xy), lim_max=np.max(xy))
 	make_heatmap_histogram(axHisty, y, tilt=True,lim_min=np.min(xy), lim_max=np.max(xy))
 	
 	if p != "fp":
		axScatter.scatter(x, y, c=z, s=14, edgecolor='', label="pearsons: " + str(ro)[:5] +"\n GRO mean: " + str(np.mean(y))[:5] + "\n ChIP mean: " + str(np.mean(x))[:5] )
	else:
		xs 		= dict([(X,1)for X in x])
		xs 		= xs.keys()
		axScatter.scatter(xs, [np.mean([ Y for X,Y in zip(x,y) if X == XX]) for XX in xs]  )	
		

	axScatter.set_xlabel("ChIP Pol II, parameter : " + p)
	axScatter.set_ylabel("GRO-seq parameter: " + p)
	axScatter.legend(loc=(0.1,0.85))	
	axScatter.grid()
	axScatter.set_xlim(np.min(x), np.max(x))
	axHistx.set_xlim(np.min(x), np.max(x))
	axHisty.set_ylim(np.min(y), np.max(y))
	
	axScatter.set_ylim(np.min(y), np.max(y))
	
	
	plt.show()	
def boxplots_strand_footprint(GRO, ChIP):
	F 	= plt.figure(figsize=(15,10))
	ax 	= F.add_subplot(2,1,1)
	ax2 = F.add_subplot(2,1,2)
	
	pos_GRO 	= [b.fp for chrom in GRO for b in GRO[chrom] if b.TSS_strand == "+"]
	neg_GRO 	= [b.fp for chrom in GRO for b in GRO[chrom] if b.TSS_strand == "-" ]
	eRNA_GRO 	= [b.fp for chrom in GRO for b in GRO[chrom] if b.TSS_strand is None ]

	pos_ChIP 	= [b.fp for chrom in ChIP for b in ChIP[chrom] if b.TSS_strand == "+"]
	neg_ChIP 	= [b.fp for chrom in ChIP for b in ChIP[chrom] if b.TSS_strand == "-" ]
	eRNA_ChIP 	= [b.fp for chrom in ChIP for b in ChIP[chrom] if b.TSS_strand is None ]
	
	ax.set_xlim(0,7)
	bp = ax.boxplot((pos_GRO, neg_GRO, eRNA_GRO),patch_artist=True, positions=(1,2,3))
	bp2 = ax2.boxplot((pos_ChIP, neg_ChIP, eRNA_ChIP),patch_artist=True, positions=(1,2,3))
	plt.show()
	

if __name__ == "__main__":
	TSS_CORR 	= True
	param_CORR 	= False
	BOX 		= False
	ChIP_file 	= "/Users/joazofeifa/Lab/ChIP/HCT116/EMG_out_files/refined_peaks.bed"
	GRO_file 	= "/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/DMSO2_3-1_bidirectional_hits_intervals.bed"
	TSS_file 	= "/Users/joazofeifa/Lab/genome_files/RefSeqHG19.txt"
	
	ChIP 		= load_file(ChIP_file, CHIP=True)	
	GRO 		= load_file(GRO_file, CHIP=False)
	match_up(GRO, ChIP)
	if BOX:
		TSS 		= load_TSS(TSS_file)
		label_TSS(ChIP, TSS)
		label_TSS(GRO, TSS)

		boxplots_strand_footprint(GRO, ChIP)
	if TSS_CORR:
		TSS 		= load_TSS(TSS_file)
		label_TSS(ChIP, TSS)
		label_TSS(GRO, TSS)

		correlation_TSS(ChIP, GRO)
	if param_CORR:
		TSS 		= load_TSS(TSS_file)
		label_TSS(ChIP, TSS)
		label_TSS(GRO, TSS)
		correlations_parameters(ChIP, GRO, p="pi", LOG=False)
	