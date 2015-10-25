import expression_correlation as ec
import time
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import matplotlib as mpl
import matplotlib.cm as cm
import math , os
from matplotlib.transforms import ScaledTranslation
class collapsable:
	def __init__(self, TFS):
		self.D 	= dict([(TF,1) for TF in TFS.split(",")])
		self.ID = ",".join(self.D.keys())
		pass
	def update(self, TFS):
		insert=False
		line_array 	= TFS.split(",")
		for TF in line_array :
			if TF in self.D:
				insert=True
		if insert:
			for TF in line_array:
				self.D[TF]=1
			self.ID = ",".join(self.D.keys())
			return False
		return True


def draw_legend(ax, dims, IDS, colors, cmap,circ_tform):
	print dims
	rad 	= dims[2]/3.
	delta 	= dims[3] / float(len(IDS))

	norm = mpl.colors.Normalize(vmin=min(colors), vmax=max(colors))
	m = cm.ScalarMappable(norm=norm, cmap=cmap)
	x 		= 0.3
	y 		= 0.1
	for i in range(len(colors)):
		circle1=mpl.patches.Ellipse((x-0.1,y),0.1,0.015,color=m.to_rgba(colors[i]) )
		ax.add_artist(circle1)
		ax.text(x,y,str(i))
		y+=delta+0.01
	pass
def make_adjacency(G,display=False):
	TFS 	= dict([(tf,1) for chrom in G for b in G[chrom] for TF in b.TFS for tf in TF.split(",")   ])
	TFS 	= TFS.keys()
	IDS 	= dict([(i,TF ) for i,TF in enumerate(TFS)])
	LABELS = {}
	for L in IDS:
		LABELS[L]=str(L)
	TFS 	= dict([(TF,i) for i,TF in enumerate(TFS)])
	
	A 		= np.zeros((len(TFS), len(TFS)))
	for chrom in G:
		for b in G[chrom]:
			bTFS 		= [tf for TF in b.TFS for tf in TF.split(",") ]	
			TF_clusters = b.TFS.keys()
			for i in range(len(TF_clusters)):
				for tf_a in TF_clusters[i].split(","):
					for j in range(i+1, len(TF_clusters)):
				
						for tf_b in TF_clusters[j].split(","):
							u,v 	= TFS[tf_a],TFS[tf_b]
							A[u,v]+=1
							A[v,u]+=1
	mean 	= np.mean([A[u,v] for u in range(A.shape[0]) for v in range(A.shape[0])])
	std 	= np.std([A[u,v] for u in range(A.shape[0]) for v in range(A.shape[0])])
	B 		= np.zeros((len(TFS), len(TFS)))
	for u in range(A.shape[0]):
		for v in range(A.shape[0]):
			if A[u,v] > (mean + std):
				B[u,v]=1
	G=nx.Graph()
	colors=list()
	colors2 = list()
	cmap = plt.get_cmap('Blues_r')
	norm = mpl.colors.Normalize(vmin=np.min([sum(A[u,:]) for u in range(A.shape[0]) ]), vmax=np.max([sum(A[u,:]) for u in range(A.shape[0]) ]))
	m = cm.ScalarMappable(norm=norm, cmap=cmap)
	cmap = plt.get_cmap('cool')
	
	for u in range(A.shape[0]):
		for v in range(u,A.shape[0]):
			G.add_node(u)
			G.add_node(v)
				
			if (A[u,v] ) > (mean + std ):	
				if display:
					G.add_edge(u,v,weight=0.0005 )
				else:
					G.add_edge(u,v,weight=0.0005 )
					
				colors2.append(m.to_rgba(max(sum(A[u,:]), sum(A[v,:]) )))
				
	for node in G.nodes():
		colors.append(sum(A[node,:]))
	#show only some labels
	cc 		= [(c,i) for i,c in enumerate(colors)]
	cc.sort()
	j 		= len(cc)
	for k,(c,i) in enumerate(cc):
		LABELS[i] 	=str(j)
		j-=1

	if display:
	
		F 	= plt.figure(figsize=(15,10))
		ax 	= F.add_axes([0.1,0.1,0.7,0.8])
		axH 	= F.add_axes([0.9,0.03,0.025,0.93])
		displat_netwrork(ax, axH, G, LABELS, colors, colors2, cmap, cc, IDS)
		plt.show()
		
	

	return A,B, IDS, TFS,G,LABELS, colors, colors2,cc
def displat_netwrork(ax, axH, G, LABELS, colors, colors2, cmap, cc, IDS,draw_heat=True):
	pos = nx.spring_layout(G)
	ax.set_xticks([])
	ax.set_yticks([])
	
	nx.draw_networkx_nodes(G, pos=pos, font_size=8, node_size=300,cmap=cmap, node_color=colors, ax=ax)
	nx.draw_networkx_edges(G, pos=pos, edge_color=colors2, alpha=0.7,width=1,ax=ax)
	nx.draw_networkx_labels(G,pos,LABELS,font_size=7, ax=ax)
	if draw_heat:
		H 		= np.zeros((1,len(colors)))
		
		H[0,:] 	= [c for c,i in cc[::-1]]
		heatmap = axH.pcolorfast(H, cmap=cmap)
		
		RIGHT 	= np.arange(1,len(IDS) ,2)
		LEFT 	= np.arange(0,len(IDS) ,2)
		axH.set_xticks(LEFT[::-1])
		axH.xaxis.set_tick_params(labeltop='on')
			
		axH.set_xticklabels([IDS[cc[i][1]].split("_")[0][:5]  + "_" + LABELS[cc[i][1]] for i in LEFT], fontsize=8,rotation=45)
		axH.set_yticklabels([])
		
		ax2 = axH.twiny()
		ax2.pcolorfast(H,cmap=cmap)
		
		ax2.set_xticks(RIGHT[::-1])
		ax2.set_xticklabels([IDS[cc[i][1]].split("_")[0][:5]  + "_" + LABELS[cc[i][1]]  for i in RIGHT], fontsize=8,rotation=45)
		ax2.set_yticklabels([])
	


def iterate_over_files(thresh="0.05"):
	DIR 	= "/Users/joazofeifa/Lab/TF_predictions/assignments_refined/"
	LST 	= list()
	t 		= 0
	for FILE in os.listdir(DIR):

		if FILE.split("_")[-1]==thresh and "Puc" not in FILE:
			print FILE

			G 	= ec.load(DIR+FILE, test=False)
			weighted,unweigthed, IDS, TFS,G,LABELS, colors, colors2,cc 	= make_adjacency(G, display= False)
			LST.append((FILE, (weighted,unweigthed, IDS, TFS,G,LABELS, colors, colors2, cc) ))
			if t > 1:
				break
			t+=1
	return dict(LST)

def compute_centrality(LST, display=True, write_out=True):
	# degree_nodes= [(d,ID)for ID, d in zip(ec.keys(), ec.values())] 
	# degree_nodes.sort(reverse=True)
	# print [IDS[ID] for d, ID in degree_nodes]
	centralities 	= list()
	for FILE in LST:
		weighted,unweigthed, IDS, TFS,G,LABELS, colors, colors2,cc  	= LST[FILE]
		degree 		= nx.degree_centrality(G)
		cc 			= nx.closeness_centrality(G)
		bc 			= nx.betweenness_centrality(G)
		ec 			= nx.eigenvector_centrality(G)
		centralities.append((FILE, IDS, degree, cc, bc, ec))
	average_ranks 	= dict([(tf,[0,0]) for FILE in LST for tf in LST[FILE][2].values() ])

	for FILE , IDS, degree, cc, bc, ec in centralities:

		for measure in (degree, cc, bc, ec):
			nodes 	= [(d,ID) for ID, d in zip(measure.keys(), measure.values())] 
			nodes.sort()
			N 		= sum([d for d,ID in nodes])
			nodes 	= [ (  1-(sum([nodes[j][0] for j in range(0,i+1) ])/N), ID) for i,(d,ID) in enumerate(nodes) ]
			nodes.sort()
			nodes 	= [(IDS[ID], cdf) for cdf, ID in nodes]
			for tf , cdf in nodes:
				average_ranks[tf][0]+=cdf
				average_ranks[tf][1]+=1
	for tf in average_ranks:
		average_ranks[tf] 	= average_ranks[tf][0] / average_ranks[tf][1]
	final = [(average_ranks[tf], tf)  for tf in average_ranks]
	final.sort(reverse=False)
	
	if display:
		cmap = plt.get_cmap('Blues')
		F 	= plt.figure(figsize=(15,10))

		axN = F.add_axes([0.1,0.5,0.8,0.45])
		#axH = F.add_axes([0.1,0.87,0.8,0.05])
		axH 	= None
		weighted,unweigthed, IDS, TFS,G,LABELS, colors, colors2,cc 	= LST[FILE]
		axN.set_title("Network from " + FILE.split("_")[0] + " dataset")
		
		displat_netwrork(axN, axH, G, LABELS, colors, colors2, cmap, cc, IDS,draw_heat=False)
		ax 	= F.add_axes([0.1,0.1,0.8,0.3])
		ax.plot(range(len(final)), [cdf for cdf, tf in final]  )
		LEFT, RIGHT 	= np.arange(0,len(final),2), np.arange(1,len(final), 2)
		ax.xaxis.set_tick_params(labeltop='on')
		ax.set_xticks( LEFT )
		cmap = plt.get_cmap('Blues_r')
				
		ax.set_xticklabels([final[i][1].split("_")[0][:5] + "_" + str(LABELS[TFS[final[i][1]]]) if final[i][1] in TFS else final[i][1].split("_")[0][:5] + "_" +"0" for i in  LEFT], rotation=45, fontsize=10 )
		fill_xs 	= [(i, cdf) for i,(cdf, tf) in enumerate(final) if cdf < 0.7]
		norm = mpl.colors.Normalize(vmin=0, vmax=1)
		m = cm.ScalarMappable(norm=norm, cmap=cmap)
		for i,(cdf, tf) in enumerate(final):
			if (i + 1) < len(final):
				if i+2 < len(final):
					avg_color=(final[i][0]+final[i+1][0]+final[i+2][0]) /3.
				else:
					avg_color=(final[i][0]+final[i+1][0] ) /2.
					
				ax.fill_between([i,i+1], [0,0], [ final[i][0], final[i+1][0]], color=m.to_rgba(avg_color -0.2  ) ) 
		ax.text(10,0.2, "Pioneering Factors?", color="white", fontsize=20,
			bbox=dict(facecolor='none', edgecolor='white', linewidth=2, boxstyle='round,pad=0.4'))
		ax.text(70,0.2, "Cell/Tissue/Condition Specific Factors?", color="black", fontsize=20,
			bbox=dict(facecolor='none', edgecolor='black', linewidth=2, boxstyle='round,pad=0.4'))
		ax2 = ax.twiny()
		ax2.set_xticks( RIGHT )
			
		ax2.set_xticklabels([final[i][1].split("_")[0][:5] + "_" + str(LABELS[TFS[final[i][1]]]) if final[i][1] in TFS else final[i][1].split("_")[0][:5] + "_" + "0" for i in  RIGHT], rotation=45, fontsize=10 )
		

		
		ax.set_ylabel("Average Centrality Measure Rank\n(across 20 GRO-seq datasets)")
		
		plt.show()

	if write_out:
		FHW 	= open("/Users/joazofeifa/Lab/EMG/TF_predictions/TF_Ranks.csv", "w")
		for rank, TF in final:
			FHW.write(TF+"," + str(rank) + "\n")
	return final

if __name__ == "__main__":
	LST = iterate_over_files(thresh="0.05")
	centralities 	= compute_centrality(LST, display=False, write_out=True)










