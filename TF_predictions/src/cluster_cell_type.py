import os, sys
import time
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.cluster.hierarchy as sch
import matplotlib.patches as patches
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from scipy.stats import gaussian_kde
from sklearn import datasets
import sklearn.decomposition as sd
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import assign_TF as at

from scipy.cluster.hierarchy import dendrogram, linkage
import networkx as nx
def get_color(EXP, EXPS, BINARY=False):
	cell_types 	= {}
	colors 		= {}
	B 			= {}
	cell_types["Allen2014"] 	= "HCT116"
	cell_types["Andersson2014"] = "HeLa"
	cell_types["Core2014"] 		= "K562"
	cell_types["Danko_2015"] 	= "T-Cells"
	cell_types["Danko2014"] 	= "K562"
	cell_types["Duttke2015"] 	= "HeLa"
	cell_types["Fong2014"] 		= "HEK293T"
	cell_types["Gardini2014"] 	= "HeLa"
	cell_types["Hah2013"] 		= "MCF7"
	cell_types["Jin2013"] 		= "IMR90"
	cell_types["Jin2014"] 		= "LNCaP"
	cell_types["Le2013"] 		= "MCF7"
	cell_types["Li2013"] 		= "MCF7"
	cell_types["Liu2013"] 		= "HEK293T"
	cell_types["Puc2015"] 		= "LNCaP"
	cell_types["Margazzi2012"] 	= "A549"
	cell_types["Meng2014a"] 	= "B-Cells"
	cell_types["Saponaro2014"] 	= "HEK293T" 
	cell_types["Sigova2013"] 	= "H1-hESC"
	cell_types["Wang2014"] 		= "B-Cells"
	cell_types["Yang2013"] 		= "LNCaP"
	cell_types["Luo2014"] 		= "AC16"
	cell_types["Chen2014"] 		= "H1-hESC"

	colors["HeLa"] 				= "red"
	colors["T-Cells"] 			= "blue"
	colors["B-Cells"] 			= "blue"
	colors["LNCaP"] 			= "green"
	colors["MCF7"] 				= "purple"
	colors["HEK293T"] 			= "teal"

	B["HeLa"] 				= 0
	B["T-Cells"] 			= 1
	B["B-Cells"] 			= 1
	B["LNCaP"] 				= 2
	B["MCF7"] 				= 3
	B["HEK293T"] 			= 4
	B["HCT116"] 			= 5
	B["AC16"]	 			= 6
	B["IMR90"] 				= 7
	B["K562"] 				= 8
	B["A549"] 				= 9
	B["H1-hESC"] 			= 10	

	
	
	
	
	ct 	= None
	for exp in cell_types:
		if exp in EXP:

			ct 	= cell_types[exp]
			break
	if ct is None:
		if BINARY:
			return 11
		return "black"
	else:
		if BINARY:
			if ct in B:	
				return B[ct]
		else:		
			if ct in colors:
				return colors[ct]
		if BINARY:
			return 11
		return "black"
def add_legend(ax):
	most_x, most_y 	= (275,750)
	corner 	= (most_x,most_y)
	width 	= 30
	height 	= 150
	fs 		= 14
	patch 	=   patches.Rectangle(
        (corner[0], corner[1]),   # (x,y)
        width,          # width
        height,          # height
        color="red"
    )
	rx, ry = patch.get_xy()
	cx = rx + patch.get_width()/2.0
	cy = ry + patch.get_height()/2.0
	ax.annotate("HeLa", (cx, cy), color='w', 
            fontsize=fs, ha='center', va='center')

	ax.add_patch(patch)

	corner 	= (most_x,most_y-height*1.1)
	patch 	=   patches.Rectangle(
        (corner[0], corner[1]),   # (x,y)
        width,          # width
        height,          # height
        color="blue"
    )
	rx, ry = patch.get_xy()
	cx = rx + patch.get_width()/2.0
	cy = ry + patch.get_height()/2.0
	ax.annotate("T/B cell", (cx, cy), color='w', 
            fontsize=fs, ha='center', va='center')
	
	ax.add_patch(patch)
	corner 	= (most_x+width*1.025,most_y )
	patch 	=   patches.Rectangle(
        (corner[0], corner[1]),   # (x,y)
        width,          # width
        height,          # height
        color="teal"
    )
	rx, ry = patch.get_xy()
	cx = rx + patch.get_width()/2.0
	cy = ry + patch.get_height()/2.0
	ax.annotate("HEK293T", (cx, cy), color='w', 
            fontsize=fs, ha='center', va='center')
	
	ax.add_patch(patch)

	ax.add_patch(patch)
	corner 	= (most_x+width*1.025,most_y-height*1.1)
	patch 	=   patches.Rectangle(
        (corner[0], corner[1]),   # (x,y)
        width,          # width
        height,          # height
        color="green"
    )
	rx, ry = patch.get_xy()
	cx = rx + patch.get_width()/2.0
	cy = ry + patch.get_height()/2.0
	ax.annotate("LNCaP", (cx, cy), color='w', 
            fontsize=fs, ha='center', va='center')
	
	ax.add_patch(patch)

	ax.add_patch(patch)
	corner 	= (most_x,most_y-height*2.2)
	patch 	=   patches.Rectangle(
        (corner[0], corner[1]),   # (x,y)
        width,          # width
        height,          # height
        color="purple"
    )
	rx, ry = patch.get_xy()
	cx = rx + patch.get_width()/2.0
	cy = ry + patch.get_height()/2.0
	ax.annotate("MCF7", (cx, cy), color='w', 
            fontsize=fs, ha='center', va='center')
	
	ax.add_patch(patch)

	ax.add_patch(patch)
	corner 	= (most_x+width*1.025,most_y-height*2.2)
	patch 	=   patches.Rectangle(
        (corner[0], corner[1]),   # (x,y)
        width,          # width
        height,          # height
        color="black"
    )
	rx, ry = patch.get_xy()
	cx = rx + patch.get_width()/2.0
	cy = ry + patch.get_height()/2.0
	ax.annotate("Other", (cx, cy), color='w', 
            fontsize=fs, ha='center', va='center')
	
	ax.add_patch(patch)

def load_TFS_per(FILE):
	G 	= {}
	with open(FILE) as FH:
		for line in FH:
			if "#" == line[0]:
				exp 	= line[1:].strip("\n").split("_bidirectional_hits_intervals.txt")[0]
				G[exp] 	= {}
			else :
				TF, d 					= line.strip("\n").split("\t")
				uni, center, bi, w 		= d.split(",")
				G[exp][TF] 				= float(uni), float(center), float(bi), float(w)
	return G
def display_matrix(G, BINARY = pow(10,-10), SHOW =False  ):
	TFS 	= dict([(TF,list()) for exp in G for TF in G[exp]])
	TFS 	= dict([(TF, i) for i,TF in enumerate(TFS)])
	toTFS 	= dict([(TFS[TF], TF) for TF in TFS])
	
	A 		= np.zeros((len(TFS), len(G.keys())))
	EXPS 	= dict([ (exp, i) for i,exp in enumerate(G) ])
	toEXPS 	= dict([(EXPS[EXP], EXP) for EXP in EXPS])
	for exp in G:
		i 	= EXPS[exp]
		for TF in TFS:
			j 	= TFS[TF]

			k 	= int(bool(G[exp][TF][0] < BINARY))
			k 	= int(bool(G[exp][TF][3] > 0.2))
			k 	= G[exp][TF][3]*10
			# if G[exp][TF][0]:
			# 	k 	= math.log(G[exp][TF][0],10)
			# else:
			# 	k 	= -150
			A[j,i]= k
	D2 	= np.zeros((A.shape[1], A.shape[1]))

	for i in range(A.shape[1]):
	    for j in range(A.shape[1]):
	        D2[i,j] = np.sum(np.abs(A[:,i]-A[:,j]))
	if SHOW:
		# Compute and plot first dendrogram.
		fig = plt.figure(figsize=(8,8))
		ax2 = fig.add_axes([0.3,0.94,0.45,0.02])

		ax2.set_xticks([])		
		ax2.set_yticks([])		
		
		ax1 = fig.add_axes([0.09,0.1,0.2,0.8])

		Y = sch.linkage(D2, method='complete')
		Z1 = sch.dendrogram(Y, orientation='right')
		ax1.set_xticks([])
		ax1.set_yticks([])

		# Compute and plot second dendrogram.
		Y = sch.linkage(D2, method='complete')
		Z2 = sch.dendrogram(Y, no_plot=True)

		# Plot distance matrix.
		axmatrix = fig.add_axes([0.3,0.1,0.45,0.8])

		idx1 = Z1['leaves']
		print idx1
		idx2 = Z2['leaves']
		D = D2[idx1,:]
		D = D[:,idx2]
		im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=plt.cm.jet_r, vmin=0, vmax=np.max(D) )
		axmatrix.set_yticks(range(len(EXPS)))
		axmatrix.yaxis.tick_right()
		labels = [item.get_text() for item in axmatrix.get_yticklabels()]
		for j,i in enumerate(idx2):
			labels[j] 	= toEXPS[i]
		axmatrix.set_yticklabels(labels, fontsize=10 )
		for t in axmatrix.yaxis.get_ticklabels():
			c 	= get_color(t.get_text(), {})
			t.set_color(c)
		axmatrix.set_xticks([])
		cmap = mpl.cm.jet_r
		norm = mpl.colors.Normalize(vmin=np.min(D), vmax=100)

		# ColorbarBase derives from ScalarMappable and puts a colorbar
		# in a specified axes, so it has everything needed for a
		# standalone colorbar.  There are many more kwargs, but the
		# following gives a basic continuous colorbar with ticks
		# and labels.
		cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,
		                                norm=norm,
		                                orientation='horizontal')
		ax2.xaxis.tick_top()
		cb1.set_label('Number of Differing TFS')
		plt.show()
	return D2, TFS, EXPS, A, toEXPS, toTFS


def dendrogram_only(A,D, TFS, EXPS, toEXPS, toTFS):
	F 	= plt.figure(figsize=(15,10))
	ax 	= F.add_subplot(1,1,1)
	add_legend(ax)
	linkage_matrix = linkage(D, "ward")
	labels 	= [(EXPS[exp], exp) for exp in EXPS]
	labels.sort()
	labels 	= [exp for i,exp in labels]
	ddata = dendrogram(linkage_matrix,labels= labels,
 	                   color_threshold=0,leaf_rotation= 85, leaf_font_size=10 )
	xlbls = ax.get_xmajorticklabels()
	for lbl in xlbls:
		c 	= get_color(lbl.get_text(), EXPS)
		lbl.set_color(c)

	ticks = ax.get_xticks()
	ax.set_yticklabels([])
	plt.tight_layout()
	plt.show()

def sorted_feature_matrix(D, TFS, EXPS, A , toEXPS, toTFS, SHOW=True):
	
	newA 	= list()
	for i in range(A.shape[0]):
		if sum(A[i,:]) > 0 :
			newA.append([ A[i,j] for j in range(A.shape[1])])
		else:
			del toTFS[i]
	toTFS 	= [(i,toTFS[i]) for i in toTFS]
	toTFS.sort()
	toTFS 	= dict([(j,i[1]) for j,i in enumerate(toTFS) ])
	A 		= np.array(newA)
	D2 	= np.zeros((A.shape[0], A.shape[0]))
	for i in range(A.shape[0]):
	    for j in range(A.shape[0]):
	        D2[i,j] = np.sum(np.abs(A[i,:]-A[j,:]))
	if SHOW:
		fig = plt.figure(figsize=(8,8))
		ax2 = fig.add_axes([0.3,0.94,0.45,0.02])

		ax2.set_xticks([])		
		ax2.set_yticks([])		
		
		ax1 = fig.add_axes([0.09,0.1,0.2,0.8])

		Y = sch.linkage(D2, method='average')
		Z1 = sch.dendrogram(Y, orientation='right')
		ax1.set_xticks([])
		ax1.set_yticks([])

		# Compute and plot second dendrogram.
		Y = sch.linkage(D2, method='average')
		Z2 = sch.dendrogram(Y, no_plot=True)

		# Plot distance matrix.
		axmatrix = fig.add_axes([0.3,0.1,0.35,0.8])

		idx1 = Z1['leaves']

		idx2 = Z2['leaves']
		D2 = D2[:,idx2]
		D2 = D2[idx2,:]

		im = axmatrix.matshow(D2, aspect='auto', origin='lower', cmap=plt.cm.jet, vmin=0, vmax=np.max(D2) )
		axmatrix.set_yticks(np.arange(0,len(toTFS),10))
		axmatrix.yaxis.tick_right()
		labels = [item.get_text() for item in axmatrix.get_yticklabels()]
		for i,j 	in enumerate(np.arange(0,len(toTFS),10)):
			l 			= ",".join([toTFS[idx2[k]][:4] for k in range(j,min(j+10, len(toTFS)))])
			labels[i] 	= l
		axmatrix.set_yticklabels(labels, fontsize=10 )
		for t in axmatrix.yaxis.get_ticklabels():
			c 	= get_color(t.get_text(), {})
			t.set_color(c)
		axmatrix.set_xticks([])
		cmap = mpl.cm.jet
		norm = mpl.colors.Normalize(vmin=np.min(D), vmax=np.max(D2))

		# ColorbarBase derives from ScalarMappable and puts a colorbar
		# in a specified axes, so it has everything needed for a
		# standalone colorbar.  There are many more kwargs, but the
		# following gives a basic continuous colorbar with ticks
		# and labels.
		cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,
		                                norm=norm,
		                                orientation='horizontal')
		ax2.xaxis.tick_top()
		cb1.set_label('Similarity in Experiment Profile')
		if SHOW:
			plt.show()
	return A, D2,toTFS

def get_most_likelicelltype(i, A, toEXPS, toTFS):
	average 	= sum([ get_color(toEXPS[u], toEXPS, BINARY=True)*A[i,u] for u in range(A.shape[1]) ]) / sum(A[i,:])
	return average
def logistic(x):
	return 1- (1.0 / (1 + math.exp(-x)))
def mode(X):
	G={}
	for i in X:
		if i not in G:
			G[i]=0
		G[i]+=1
	return max([(G[x], x) for x in G ])[1]

def get_cell_type_label(i):
	B 		= {}
	B["HeLa"] 				= 0
	B["T-Cells"] 			= 1
	B["B-Cells"] 			= 1
	B["LNCaP"] 				= 2
	B["MCF7"] 				= 3
	B["HEK293T"] 			= 4
	B["HCT116"] 			= 5
	B["AC16"]	 			= 6
	B["IMR90"] 				= 7
	B["K562"] 				= 8
	B["A549"] 				= 9
	B["H1-hESC"] 			= 10	
	lbl 	= ""
	for j in B:
		if i==B[j]:
			lbl+=j
	if not lbl:
		lbl="other"

	return lbl

def finally_as_a_network(A, D2, toEXPS, toTFS):
	#toTFS by specificy
	temp 	= [(sum(A[i,:]), i)    for i in toTFS]
	temp.sort()
	toTFS 	= dict([ (j, toTFS[y] ) for j,(x,y) in enumerate(temp)])
	ro 		= [i for  s,i in temp ]
	A 		= A[ro ,:]		
	D2 		= D2[ro,:]
	D2 		= D2[:,ro]


	G=nx.Graph()
	cmap 	= plt.get_cmap('Blues')
	norm 	= mpl.colors.Normalize(vmin=0, vmax=3)
	m2 		= cm.ScalarMappable(norm=norm, cmap=cmap)

	toTFS2 	= dict([(i,i)for i in toTFS])
	colors 	= list()
	colors2 = list()
	for u in range(D2.shape[0]):
		colors2.append(sum(A[u,:]))
		for v in range(u,D2.shape[0]):
			#U,V 	= get_most_likelicelltype(u,A, toEXPS, toTFS), get_most_likelicelltype(v,A, toEXPS, toTFS)

			if D2[u,v] <4  and u!=v:
				G.add_node(u)
				G.add_node(v)
				G.add_edge(u,v,weight=0.002 )
				colors.append(m2.to_rgba(1))
	F = plt.figure(figsize=(8,8))
	ax 		= F.add_axes([0.05,0.5,0.5,0.45])
	ax2 	= F.add_axes([0.65,0.5,0.015,0.45])
	
	ax3 	= F.add_axes([0.05,0.05,0.5,0.4])
	ax4 	= F.add_axes([0.65,0.05,0.015,0.4])
	
	pos = nx.spring_layout(G)
	ax.set_xticks([])
	ax.set_yticks([])
	ax2.yaxis.tick_right()
	ax2.set_xticks([])
	vmin 	= min(colors2)
	vmax 	= max(colors2)
	norm = mpl.colors.Normalize(vmin=vmin, vmax= vmax)
	cmap = plt.get_cmap('jet' )
	m = cm.ScalarMappable(norm=norm, cmap=cmap)
	bounds 	= np.arange(0,len(toTFS), 7)
	HM 		= np.zeros((A.shape[0], 1))
	HM[:,0] = [sum(A[i,:]) for i in range(A.shape[0])]
	ax2.matshow( HM, aspect='auto', origin='lower', cmap=plt.cm.jet, vmin=1, vmax=np.max(HM), alpha=1 )
	ax2.set_xticks([])
	ax2.set_yticks(bounds)
	ax2.set_yticklabels([ str(i) + "-" +str(min(i+7, len(toTFS))) + ":" +  ",".join([toTFS[j][:5] for j in range(i,min(i+7, len(toTFS)))  ])  for i in  bounds ] , fontsize=10)
	
	nx.draw_networkx_nodes(G, pos=pos, font_size=8, node_size=100, ax=ax, alpha=1, node_color=[m.to_rgba(colors2[u]) for u in G.nodes()], 
		line_color= [m.to_rgba(colors2[u]) for u in G.nodes()] )
	nx.draw_networkx_edges(G, pos=pos,  alpha=0.7,width=3,ax=ax,edge_color=colors)

	ax.set_title("Number of Times a TF is Found in an Cell Type")

	#nx.draw_networkx_labels(G,pos,toTFS2,font_size=12, ax=ax)
	ax5 = ax2.twinx()
	ax5.yaxis.tick_left()
	ax5.set_xticks([])
	ax5.set_yticks(np.arange(0, len(toTFS), 20))
	#ax3.set_title("Percentage of GRO-seq experiments (n=36)\nwhere a TF is considered 'active' ")
	ax5.set_yticklabels([str(100*np.mean([ sum(A[j,:]) for j in range(i,min(i+20, len(toTFS))) ]) / float(36.)  )[:4] + "%" for i in np.arange(0,len(toTFS), 20) ])
	ax2.yaxis.tick_right()
	t=0
	fixed_positions=list()
	colors3 		= list()
	for i in range(A.shape[0]):
		if sum(A[i,:]) > 15 and i in G.nodes():
			for j in G.edges(i):
				G.remove_edge(j[0],j[1])

			G.remove_node(i)
		elif i in G.nodes():
			fixed_positions.append((i, (pos[i][0],pos[i][1] )))


		if i in G.nodes():
			t+=1
	NS 	= G.nodes()
	for e in NS:
		if not G.edges(e):
			G.remove_node(e)
	for e in G.edges():
		u 	= mode([ get_color(toEXPS[j], {}, BINARY=True) for j in range(A.shape[1]) if A[e[0],j] ]  ) 
		v 	= mode([ get_color(toEXPS[j], {}, BINARY=True) for j in range(A.shape[1]) if A[e[1],j] ]  ) 
#		if (u-v)==0:
#			G.add_edge(e[0],e[1],weight=0.001/(abs(u-v)+1) )
		G.add_edge(e[0],e[1],weight=1 )
	fixed_positions 	= dict(fixed_positions)
	#get_color(EXP, EXPS, BINARY=False)
	new_nodes_colors 	= [ mode([ get_color(toEXPS[j], {}, BINARY=True) for j in range(A.shape[1]) if A[n,j] ]  )   for n in G.nodes() ]

	new_nodes_colors3 	= [ (mode([ get_color(toEXPS[j], {}, BINARY=True) for j in range(A.shape[1]) if A[n,j] ]  ) , toTFS[n])   for n in G.nodes() ]
	new_nodes_colors2 	= {}
	for exp, tf in new_nodes_colors3:
		if exp not in new_nodes_colors2:
			new_nodes_colors2[exp]=list()
		new_nodes_colors2[exp].append(tf)
	norm = mpl.colors.Normalize(vmin=min(new_nodes_colors), vmax=max(new_nodes_colors))
	cmap = plt.get_cmap('jet' )
	m = cm.ScalarMappable(norm=norm, cmap=cmap)
	UNI 	= dict([(i,1) for i in new_nodes_colors])
	UNI 	= [i for i in UNI]
	new_nodes_colors 	= [m.to_rgba(i) for i in new_nodes_colors]

	HM 		= np.zeros((len(UNI), 1))
	HM[:,0] = [i for i in UNI ]
	ax4.matshow( HM, aspect='auto', origin='lower', cmap=plt.cm.jet, vmin=np.min(HM), vmax=np.max(HM), alpha=1 )
	ax4.set_xticks([])
	ax4.set_yticks(range(len(UNI)))
	ax4.set_yticklabels([ get_cell_type_label(i)  for i in UNI ])
	ax6 = ax4.twinx()
	ax6.yaxis.tick_right()
	ax6.set_xticks([])
	ax6.set_yticks(range(len(UNI)) )
	ax6.set_yticklabels([ "\n".join([",".join([new_nodes_colors2[n][k] for k in range(i,min(i+8, len(new_nodes_colors2[n]))) ] )  for i in np.arange(0, len(new_nodes_colors2[n]) , 8) ])  for n in new_nodes_colors2 ],
		fontsize=10)
		



	pos 	= nx.spring_layout(G)
	ax3.set_title("Cell Type Specific TF Coloring")
	nx.draw_networkx_nodes(G, pos=fixed_positions, font_size=8, node_size=100, ax=ax3, alpha=1, node_color=new_nodes_colors )
	nx.draw_networkx_edges(G, pos=fixed_positions,  alpha=0.7,width=3,ax=ax3,edge_color=[m2.to_rgba(1) for e in G.edges()] )
	ax3.set_xticks([])
	ax3.set_yticks([])
	
	plt.show()
def get_GC_motif(D, EXP):
	X 			= np.zeros((4,))
	N 			= 0.0
	for exp in D:
		if EXP in exp:
			A 	= D[exp]
			for i in range(A.shape[0]):
				X+=A[i,:]
				N+=1
	if N:
		X/=N
		return sum(X[1:3])/sum(X)
	return None
def compare_entropy(G):
	FILE 	= "/Users/joazofeifa/Lab/TF_predictions/HOCOMOCOv9_AD_MEME.txt"
	D 		= at.load_PSSMs(FILE,test=False)
		
	xy 		= list()
	for exp in G:
		if np.random.uniform(0,1) < 0.5:
			for TF in G[exp]:
				GC 	= get_GC_motif(D, TF)
				if GC and G[exp][TF][0]:
					xy.append((G[exp][TF][0],G[exp][TF][3], GC) )
		
	F 		= plt.figure()
	ax 		= F.add_axes([0.1, 0.1, 0.35, 0.35])
	axc 	= F.add_axes([0.5, 0.1, 0.05, 0.35])
	ax3 	= F.add_axes([0.6, 0.55, 0.35, 0.35])
	ax2 	= F.add_axes([0.1, 0.55, 0.35, 0.35])
	
	cmap 	= plt.get_cmap('Blues')

	x,y  	=  [math.log(u, 10) for u,v,z in xy if  v > 0.005   ],[v for u,v,z in xy  if   v> 0.005  ]

	c 		= [z for u,v,z in xy if v > 0.005  ] 
	print len(c), len(y)
	norm 	= mpl.colors.Normalize(vmin=min(c), vmax=max(c))
	m 		= cm.ScalarMappable(norm=norm, cmap=cmap)
	cc 		= [m.to_rgba(C) for C in c]
	
	cb1 = mpl.colorbar.ColorbarBase(axc, cmap=cmap,
	                        norm=norm,
	                        orientation='vertical')
	cb1.set_label('GC Content')
	ax.scatter(x,y, color=cc, edgecolor='' )
	ax.set_xlabel("KS-Test")
	ax.set_ylabel("True Positive Rate")
	
	ax.grid()

	ax2.hist(c, bins=25)
	ax2.set_xlabel("GC Content for All HOCOMOCO Motifs")
	ax2.set_ylabel("Frequency")
	
	XY = np.vstack([c,y])
	
	z = gaussian_kde(XY)(XY)
	
	
	ax3.scatter(c,y, s=5, edgecolor='')
	ax3.set_xlabel("GC Content")
	ax3.set_ylabel("True Positive Rate")
	
	ax3.grid()


	plt.show()
	pass
def PCA_plot(D, TFS, EXPS, A , toEXPS, toTFS):
	A 	= A.T
	pca = sd.PCA(n_components=2)
	X_r = pca.fit(A).transform(A)
	F 	= plt.figure(figsize=(15,10))
	ax 	= F.add_subplot(111)
	y 	= [get_color(toEXPS[i], EXPS, BINARY=True) for i in range(X_r.shape[0]) ]
	lda = LinearDiscriminantAnalysis(n_components=2)
	X_r2 = lda.fit(D, y).transform(D)

	ax.scatter(X_r[:,0], X_r[:,1], c=[get_color(toEXPS[i], EXPS) for i in range(X_r.shape[0])], s=150 )
	plt.show()

def PCA_plot_TFS(D, TFS, EXPS, A, toEXPS, toTFS):
	
	
	pca = sd.PCA(n_components=3)
	X_r = pca.fit(A).transform(A)
	SPECS 	= list()
	for i in range(A.shape[0]):
		SPECS.append(sum(A[i,:]))
	norm = mpl.colors.Normalize(vmin=min(SPECS), vmax= max(SPECS) )
	cmap = plt.get_cmap('jet' )
	m = cm.ScalarMappable(norm=norm, cmap=cmap)
	
	F 	= plt.figure(figsize=(15,10))
	ax 	= F.add_subplot(111, projection="3d")
	
	ax.scatter(X_r[:,0], X_r[:,1], X_r[:,2], s=100, c=[m.to_rgba(SPEC) for SPEC in SPECS] )
	plt.show()
	
	pass

if __name__ == "__main__":
	FILE 	= "/Users/joazofeifa/Lab/EMG/TF_predictions/files/TFS_PER.tsv"
	G 		= load_TFS_per(FILE)
	# compare_entropy(G)
	D, TFS, EXPS, A , toEXPS, toTFS	= display_matrix(G, BINARY=pow(10,-10), SHOW=False)
	PCA_plot(D, TFS, EXPS, A , toEXPS, toTFS)
	A, D2,toTFS= sorted_feature_matrix(D, TFS, EXPS, A , toEXPS, toTFS, SHOW=False)

#	finally_as_a_network(A, D2, toEXPS, toTFS)
	dendrogram_only(A,D, TFS, EXPS, toEXPS, toTFS)
