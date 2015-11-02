import assign_TF as at
import numpy as np
from math import log
ALPHA 	= {0: "A",  1: "C", 2: "G",3: "T" }
import matplotlib as mpl
from scipy.cluster.hierarchy import dendrogram, linkage as linkage_scipy
import matplotlib.pyplot as plt
import time
import networkx as nx
import matplotlib.cm as cm

def LOG(x):
	if x <=0:
		return -np.inf
	return log(x,2)
def compare_two(a,b, prior=0.00001):
	
	for i in range(a.shape[0]):
		a[i,:]+=prior
		a[i,:]/=sum(a[i,:])
	for i in range(b.shape[0]):
		b[i,:]+=prior
		b[i,:]/=sum(b[i,:])
	
	#look at all offshifts
	if a.shape[0] < b.shape[0]:
		A 	= a
		B 	= b
	else:
		A 	= b
		B 	= a
	OFF 	= B.shape[0]-A.shape[0]
	KL 		= np.inf
	HAMMING = np.inf
	KL_OFF 	= None
	H_OFF 	= None
	for off in range(0,B.shape[0]-A.shape[0]+1):
		h,kl 	= 0,0
		for j in range(A.shape[0]):
			h+=abs(A[j,:].argmax() - B[j+off,:].argmax() )
			aa=sum([A[j,k]*LOG((A[j,k] ) / (B[j+off,k] ) ) for k in range(4) ])
			bb=sum([B[j+off,k]*LOG((B[j+off,k] ) / (A[j,k] ) ) for k in range(4) ])
			kl+=(aa+bb)/2.
		kl/=A.shape[0]
		if h < HAMMING:
			H_OFF 	= off
			HAMMING = h
		if kl < KL:
			KL_OFF 	= off
			KL 		= kl
	kl_avg 	= np.zeros(A.shape)
	h_avg 	= np.zeros(A.shape)
	for j in range(A.shape[0]):
		for k in range(4):
			kl_avg[j,k]=A[j,k]+B[j+KL_OFF,k]
			h_avg[j,k]=A[j,k]+B[j+H_OFF,k]
		kl_avg[j,:]/=sum(kl_avg[j,:])
		h_avg[j,:]/=sum(kl_avg[j,:])




	return HAMMING, KL,kl_avg, h_avg

def perform_linkage(D, OUT=""):
	H 	= np.zeros((len(D.keys()),len(D.keys())  ))

	KL 	= np.zeros(H.shape)
	toM = dict([(i,m)for i, m in enumerate(D.keys())])
	toI = dict([(toM[i],i) for i in toM])
	FHW_KL=open(OUT+ "kl_distance_matrix.csv","w")
	FHW_H=open(OUT+ "hamming_distance_matrix.csv","w")
	
	for i in range(len(D)):
		for j in range(i+1, len(D)):
			h,kl 	= compare_two(D[toM[i]], D[toM[j]])
			H[i,j] 	= h
			H[j,i] 	= h
			KL[i,j] = kl
			KL[j,i] = kl
	header 	= "#"+ ",".join([str(i)+":"+str(toM[i]) for i in toM])
	FHW_H.write(header+"\n")
	FHW_KL.write(header+"\n")
	for i in range(KL.shape[0]):
		kl_line, h_line 	="",""
		for j in range(KL.shape[0]):
			kl_line+=str(KL[i,j])+","
			h_line+=str(H[i,j])+","
		kl_line,h_line 	= kl_line.strip(","), h_line.strip(",")
		FHW_KL.write(kl_line+"\n")
		FHW_H.write(h_line+"\n")
	FHW_H.close()
	FHW_KL.close()





	
	return 	KL,H
def linkage(KL,H, toM, SHOW=True):

	linkage_matrix = linkage_scipy(KL, "complete")
	labels 	= [(i,toM[i]) for i in toM]
	labels.sort()
	labels 	= [y for x,y in labels]
	ddata = dendrogram(linkage_matrix,labels=labels,
 	                   color_threshold=0,leaf_rotation= 85, leaf_font_size=9 )
	if SHOW:
		plt.tight_layout()
		plt.show()
	return ddata
def load_matrix(FILE):
	header 	=True
	M 	 	= {}
	L 		= list()
	with open(FILE) as FH:
		for line in FH:
			if not header:
				L.append([float(x) for x in line.strip("\n").split(",")])
			else:
				M 		= dict([(int(p.split(":")[0]),p.split(":")[1]) for p in line[1:].strip("\n").split(",")])
				header 	= False
	return np.array(L),M
def try_different_thresholds(A,toM, OUT="",res=2):
	thresholds 	= np.linspace(0.001,5.5,res)
	FHW 		= open(OUT+"thresholds_nearest_assignments.tsv", "w")
	for i,threshold in enumerate(thresholds):
		print float(i)/res
		FHW.write("#" + str(threshold) + "\n")
		line 	= ""
		for i in range(A.shape[0]):
			line 	= toM[i] + "\t"
			for j in range(A.shape[0]):
				if i != j:
					if A[i,j] < threshold:
						line+=str(toM[j]) + ","
			line 	= line.strip(",") + "\n"
			FHW.write(line)
	FHW.close()

def get_counts(FILE, show_thresh=True):
	ts 	= list()
	X 	= list()
	xs 	= list()

	with open(FILE) as FH:
		for line in FH:
			if "#" == line[0]:
				if xs:
					X.append(xs)
				xs 	= list()
				t 	= float(line[1:].strip("\n"))
				ts.append(t)
			else:
				line_array 	= line.strip("\n").split("\t")
				xs.append((line_array[0].split("_")[0], [m.split("_")[0] for m in line_array[1].split(",") if m]))
		if xs:
			X.append(xs)

	if show_thresh:
		F 		= plt.figure(figsize=(15,10))
		ax1 	= F.add_subplot(1,1,1)
		ax1.grid()
		means 	= [np.mean([float(len(x)) for y,x in xs ])  for xs in X]
		stds 	= [np.std([float(len(x)) for y,x in xs ])  for xs in X]
		ax1.scatter(ts, [np.mean([float(len(x)) for y,x in xs ])  for xs in X])
		ax1.fill_between(ts, [m-stds[i] for i,m in enumerate(means)],[m+stds[i] for i,m in enumerate(means)], color="grey", alpha=0.3 )
		ax1.set_xlabel("KL Divergence Threshold")
		ax1.set_ylabel("Average Number of Connections")
		plt.show()
	return ts, X
def threshold_as_network(X):
	X 		= dict(X)
	IDS 	= dict([(x,i) for i,x in enumerate(X.keys())])
	LABELS 	= dict([(IDS[x], x) for x in IDS ])
	A 	= np.zeros((len(IDS), len(IDS)))
	for u in X:
		for v in X[u]:
			i,j 	= IDS[u], IDS[v]
			A[i,j]+=1
	G=nx.Graph()
	cmap 	= plt.get_cmap('Blues')
	norm 	= mpl.colors.Normalize(vmin=0, vmax=3)
	m 		= cm.ScalarMappable(norm=norm, cmap=cmap)
	colors 	= list()
	for u in range(A.shape[0]):
		for v in range(u,A.shape[0]):
			G.add_node(u)
			G.add_node(v)
			if A[u,v]:
				G.add_edge(u,v,weight=0.1*A[u,v] )
				colors.append(m.to_rgba(1))
	F 	= plt.figure(figsize=(15,10))
	ax 	= F.add_subplot(1,1,1)
	pos = nx.spring_layout(G)
	ax.set_xticks([])
	ax.set_yticks([])
	
	nx.draw_networkx_nodes(G, pos=pos, font_size=8, node_size=100, ax=ax, alpha=0.2)
	nx.draw_networkx_edges(G, pos=pos,  alpha=0.7,width=3,ax=ax,edge_color=colors)
	nx.draw_networkx_labels(G,pos,LABELS,font_size=7, ax=ax)
	plt.show()



if __name__ == "__main__":
	make_distance 	= False
	OUT 			= "/Users/joazofeifa/Lab/EMG/TF_predictions/files/"
	SCIPY 			= False
	CUSTOM 			= False
	ts, xs 			= get_counts(OUT+"thresholds_nearest_assignments.tsv", show_thresh=False)
	threshold_as_network(xs[5])
	if make_distance:
		FILE 	= "/Users/joazofeifa/Lab/TF_predictions/HOCOMOCOv9_AD_MEME.txt"
		D 		= at.load_PSSMs(FILE,test=False)
		perform_linkage(D, OUT=OUT)
	if CUSTOM:
		FILE 	= "/Users/joazofeifa/Lab/TF_predictions/HOCOMOCOv9_AD_MEME.txt"
		KL,M 	= load_matrix(OUT+"kl_distance_matrix.csv")
		
		try_different_thresholds(KL, M, OUT=OUT, res=50)
	
	if SCIPY:
		KL,M 	= load_matrix(OUT+"kl_distance_matrix.csv")
		H,M 	= load_matrix(OUT+"hamming_distance_matrix.csv")
		SCIPY 	= True
		if SCIPY:
			P 	 	= linkage(KL,H,M, SHOW=True)
			extract_clusters(P, threshold=None)


