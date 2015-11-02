import expression_correlation as ec
import time
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.cluster.hierarchy as sch
import assign_TF as at
import matplotlib as mpl
def load_ranks(FILE):
	G 	= {}
	with open(FILE) as FH:
		for line in FH:
			TF,rank 	= line.strip("\n").split(",")
			G[TF] 		= (float(rank))
	return G
def make_counts(G,ranks,ax,FILTER1=None, FILTER2=None):
	TFS 	= dict([(tf,ranks[tf]) if tf in ranks else (tf,1.) for chrom in G for b in G[chrom] for clust in b.TFS for tf in clust.split(",")])
	TFS 	= [(TFS[tf],tf) for tf in TFS]
	TFS.sort()
	TFS 	= dict([(tf,i) for i,(r,tf) in enumerate(TFS) ])
	A 		= np.zeros((len(TFS), len(TFS)))
	NS 		= np.zeros((len(TFS), len(TFS)))
	for chrom in G:
		for b in G[chrom]:
			tfs 	= [tf for clust in b.TFS for tf in clust.split(",")]
			for i in range(len(tfs)):
				for j in range(i, len(tfs)):
					tfa, tfb 	= tfs[i], tfs[j]
					u,v 		= TFS[tfa],TFS[tfb]
					if FILTER1 is None and FILTER2 is None or (FILTER1 and not FILTER1[chrom].searchInterval((b.start, b.stop))) or (FILTER2 and FILTER2[chrom].searchInterval((b.start, b.stop))) :
				
						A[u,v]+=math.log(b.get_density()+1)
						A[v,u]+=math.log(b.get_density()+1)
						#A[u,v]+=b.w
						#A[v,u]+=b.w
						
					NS[u,v]+=1
					NS[v,u]+=1
					
	for u in range(A.shape[0]):
		for v in range(A.shape[0]):
			if NS[u,v]:
				A[u,v]/=NS[u,v]
	Y = sch.linkage(NS, method='centroid')
	Z1 = sch.dendrogram(Y,no_plot=True)
	
	idx1 = Z1['leaves']
	A = A[idx1,:]
	A = A[:,idx1]
	heatmap = ax.pcolor(A, cmap=plt.cm.Blues,vmin=0, vmax=np.max(A))
	labels 	= dict([(TFS[tf], tf) for tf in TFS])
	xs 	= np.arange(0,len(labels),2)
	ys 	= np.arange(1,len(labels),2)
	ax.set_xticks(xs)
	ax.set_yticks(ys)
	
	ax.set_xticklabels([labels[x].split("_")[0][:6] for x in xs], minor=False,rotation=45, fontsize=10)
	ax.set_yticklabels([labels[x].split("_")[0][:6] for x in ys], minor=False,rotation=45, fontsize=10)
def compares_all(G,ranks, TSS_T):
	F 	= plt.figure(figsize=(15,10))
	ax1 = F.add_axes([0.05, 0.15, 0.35, 0.8])
	ax1.set_title("Bidirectionals not in Promoters (eRNAs)\nPaused Probability")
	ax2 = F.add_axes([0.5, 0.15, 0.35, 0.8])
	make_counts(G,ranks,ax1,FILTER1=None, FILTER2=None)
	ax2.set_title("Bidirectionals in Promoters (genes)\Paused Probability")
	make_counts(G,ranks,ax2,FILTER1=None, FILTER2=None)
	cmap = mpl.cm.Blues
	norm = mpl.colors.Normalize(vmin=0, vmax=1)
	ax3 = F.add_axes([0.9, 0.15, 0.05, 0.8])

	cb1 = mpl.colorbar.ColorbarBase(ax3, cmap=cmap,
                                   norm=norm,
                                   orientation='vertical')
	plt.show()


if __name__ == "__main__":
	TSS_T 	= at.get_TSS_tree("/Users/joazofeifa/Lab/genome_files/TSS.bed")
		
	TF_RANKS 	= "/Users/joazofeifa/Lab/EMG/TF_predictions/TF_Ranks.csv"
	ranks 		= load_ranks(TF_RANKS)
	eRNA_motif  = "/Users/joazofeifa/Lab/TF_predictions/assignments_refined/Allen2014_DMSO2_3-1_0.05"
	B  			= ec.load(eRNA_motif,test=True)

	compares_all(B, ranks, TSS_T)
		