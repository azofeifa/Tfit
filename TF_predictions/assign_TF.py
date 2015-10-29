import os,sys
sys.path.append("/Users/azofeifa/Lab/")
sys.path.append("/Users/joazofeifa/Lab/")
from linkage import get_TSS_tree
from interval_searcher import intervals, node
import numpy as np
import math,re
import time
try:
	import networkx as  nx
except:
	NO_NET =True
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.cluster.hierarchy import dendrogram, linkage as linkage_scipy
import matplotlib.pyplot as plt

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
def make_tree(FILE,A,M,G,TF, test=False):
	with open(FILE) as FH:
		header = True
		t 		= 0
		for line in FH:
			if not header:
				pattern_name, chrom, start, stop,strand,score,p,q 	= line.split("\t")[:8]
				if pattern_name not in M:
					M[pattern_name]=0
				M[pattern_name]+=1	
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
							A[i][TF][pattern_name].append((x -y, float(p), float(q),strand))

			else:
				header=False
	return A,M
def read_in_directory(root, Q):
	A,M 	= {},{}
	for DIR in os.listdir(root):
		if os.path.exists(root+ DIR+ "/fimo.txt" ):
			print DIR
			A,M 	= make_tree(root+ DIR+ "/fimo.txt", A, M,Q,DIR.split("_")[0] , test=False)
	return A,M
def write_out(A,M,IDS, OUT):
	FHW= open(out, "w")
	for tf_pattern in M:
		FHW.write(tf_pattern+"\t"+str(M[tf_pattern])+ "\n"  )
	FHW.write("-----------------\n")
	for i in range(1, max(IDS.keys() ) +1):
		ID 	= IDS[i]
		D 	= ID + "\t"
		if i in A:
			for TF in A[i]:
				
				for j,pattern_name in enumerate(A[i][TF]):
					D+=pattern_name+ "="
					for d,p,q,strand in A[i][TF][pattern_name]:
						D+=str(d)+"_" + str(p) + "_" + str(q) + "_" + strand +  ";"
					D=D.strip(";")
					D+=","
			D=D.strip(",")
		FHW.write(D+ "\n")
	FHW.close()
def write_out_TF_classifications(G, OUT):
	FHW 	= open(OUT, "w")
	for chrom in G:
		for b in G[chrom]:
			FHW.write(chrom+"\t" + str(b.start) + "\t" + str(b.stop) + "\t" + b.info + "\t")
			I=""
			for TF in b.TFS:
				if abs(b.TFS[TF]) < 100:
					I+=TF+":"+str(b.TFS[TF])+","
			I.strip(",")
			FHW.write(I+"\n")
	FHW.close()

class bidir:
	def __init__(self, chrom,start, stop, info, distances):
		self.chrom,self.start, self.stop, self.info 	= chrom,int(start), int(stop), info
		self.TFS = self._parse_distances(distances)
	def _parse_distances(self,distances, FILTER=500):
		G 	= {}
		for TF in re.split(":", distances):
			if "=" in TF:
				M 	= TF.split("=")[0]
				
				DS 	= [float(x) for x in TF.split("=")[1].split(",")[:-1]]
				DS 	=[ (abs(d), d) for d in DS if abs(d) < FILTER]
				if DS:
					if M not in G:
						G[M] =min(DS)[1]
					else:
						G[M] = min(G[M], min(DS)[1])
		
		return G
def load_TF_distance_out(FILE, test=True,FILTER=None):
	G 	={}
	t 	= 0
	with open(FILE) as FH:
		for line in FH:
			if test and t > 5000:
				break
			t+=1
			chrom,start, stop, info, distances 	= line.strip("\n").split("\t")
			if chrom not in G:
				G[chrom] 	= list()
			if FILTER is None or not FILTER[chrom].searchInterval((int(start), int(stop))):
				G[chrom].append(bidir(chrom,start, stop, info, distances))
	return G
def co_occurance_over_eRNAs(G, display=True, threshold=100):
	for chrom in G:
		for b in G[chrom]:
			tfs 	= b.TFS.keys()
			dels 	= list()
			for i in range(len(tfs)):
				for j in range(i+1, len(tfs)):
					tf_a, tf_b 	= tfs[i], tfs[j]
					dist 		= b.TFS[tf_a] - b.TFS[tf_b]
					if abs(dist) < 20:
						if abs(b.TFS[tf_a]) < abs(b.TFS[tf_b]):
							dels.append(tf_b)
						else:
							dels.append(tf_b)
							

			for d in dels:
				if d in b.TFS:
					del b.TFS[d]


	TFS 	= dict([(TF, 1) for chrom in G for b in G[chrom] for TF in b.TFS])
	A 		= np.zeros((len(TFS), len(TFS)))
	OTFS 	= {}
	for i,TF in enumerate(TFS.keys()):
		TFS[TF] 	= i
		OTFS[i] 	= TF

	for chrom in G:
		for B in G[chrom]:
			for tf_a in B.TFS:
				for tf_b in B.TFS:
					if abs(B.TFS[tf_a]) < threshold and abs(B.TFS[tf_b]) < threshold:
					
						A[TFS[tf_a],TFS[tf_b]]+=1
	if display:
		G=nx.Graph()
		labels 	= {}
			
		for i in range(A.shape[0]):
			for j in range(i+1, A.shape[0]):
				if A[i,j] / float( A[i,i]+ A[j,j]) > 0.06:
					labels[i] = OTFS[i]
					labels[j] = OTFS[j]


					G.add_node(i)
					G.add_node(j)
					
					G.add_edge(i,j,weight=0.1*A[i,j] / float( A[i,i]+ A[j,j]))
					G.add_edge(j,i,weight=0.1*A[j,i] / float( A[i,i]+ A[j,j]))

		F 	= plt.figure(figsize=(15,10))
		pos 	= nx.spring_layout(G)
		ax 	= F.add_axes([0.1,0.1,0.4,0.8])
		ax2 = F.add_axes([0.6,0.1,0.3,0.35])
		ax2.hist([A[i,j] / float( A[i,i]+ A[j,j])  for i in range(A.shape[0]) for j in range(i+1, A.shape[0])] , bins=50 )
		ax3 = F.add_axes([0.6,0.5,0.3,0.35])
		NORMED 	= np.zeros(A.shape)
		for i in range(A.shape[0]):
			for j in range(A.shape[0]):
				NORMED[i,j] 	= A[i,j] / float(A[i,i] + A[j,j])
		cmap = mpl.cm.cool
		norm = mpl.colors.Normalize(vmin=0, vmax=1)
		axcb = F.add_axes([0.95, 0.5, 0.025, 0.35])
		cb1 = mpl.colorbar.ColorbarBase(axcb, cmap=cmap,
                                   norm=norm,
                                   orientation='vertical')
		ax3.pcolor(NORMED,cmap=cmap, vmin=0, vmax=1.0)
		nx.draw_networkx(G, pos=pos, font_size=12, node_size=25,labels=labels,alpha=0.5, ax=ax)
		plt.show()
	return A

def look_at_occupancy(G):
	threshold 	= np.linspace(0,600,25)
	occupants_all 	= [len(b.TFS) for chrom in G for b in G[chrom]  ]
	occupants_100 	= [len([1 for tf in b.TFS if abs(b.TFS[tf]) < 100 ]) for chrom in G for b in G[chrom]  ]
	occupants_10 	= [len([1 for tf in b.TFS if abs(b.TFS[tf]) < 10] ) for chrom in G for b in G[chrom]  ]
	
	
	F 			= plt.figure(figsize=(15,10))
	ax 			= F.add_subplot(1,1,1)
	ax.hist(occupants_all,bins=50, alpha=0.3, color="blue")
	ax.hist(occupants_100,bins=25, alpha=0.3, color="red")
	ax.hist(occupants_10,bins=25, alpha=0.3, color="green")
	
	plt.show()


def EM(X,Y,mu, si, f,g):
	ws 		= [0.5,0.5]
	t,T 	= 0,500
	converged=False
	f 		= lambda x,si: (1.0 / (math.sqrt(2*math.pi)*si) )*math.exp(-pow(x-mu,2)/(2*pow(si,2)))
	prev_si = si
	prev_w 	= ws[0]-0.2
	while t < T and not converged and si > 1:
		EX 	= [0,0]
		S,S2 		= 0.0,0.0

		for x,y in zip(X,Y):
			n,u 	= f(x,si)*ws[0], g(x)*ws[1]
			norm 	= n+u
			EX[0]+=(n/norm)*y
			EX[1]+=(u/norm)*y
		N 	= EX[0] + EX[1]
		ws[0] 	= EX[0] / N
		ws[1] 	= EX[1] / N
		if abs(prev_w - ws[0]) < 0.001:
			converged=True
		prev_w 	= ws[0]
		prev_si = si
		t+=1
	return ws[0],si

def FILTER(G,display=False,WW=0.15):
	mu  	= 0 
	si 		= 50
	TFS 	=dict([ (TF,list()) for chrom in G for b in G[chrom] for TF in b.TFS])
	dels 	= list()
	for TF in TFS:
		DS 		= [b.TFS[TF] for chrom in G for b in G[chrom] if TF in b.TFS]
		NULL 	= np.linspace(min(DS), max(DS),1000)
		d, pv 	= stats.ks_2samp(DS, NULL)

		if pv < 0.01:
			bins= 50
			counts,edges 	= np.histogram(DS, bins=bins, normed=1)
			edges 			= (edges[:-1] + edges[1:])/ 2.
			g 				= lambda x: 1.0 / (max(edges) - min(edges))
			f 				= lambda x,si: (1.0 / (math.sqrt(2*math.pi)*si) )*math.exp(-pow(x-mu,2)/(2*pow(si,2)))
			w,si 			= EM(edges, counts,mu,si,f,g)
			if display:
				F 	= plt.figure(figsize=(15,10))
				ax1 = F.add_subplot(1,1,1)
				
				ax1.set_title(TF + " " + str(d) + "," + str(pv)+ "," + str(w))
				ax1.bar(edges, counts,width=(edges[-1]- edges[0])/bins, alpha=0.75)
				xs 				= np.linspace(min(DS), max(DS),1000)
				ys 				= [g(x)*(1-w) + f(x,si)*w for x in xs]
				ax1.plot(xs,ys, linewidth=4., color="black", linestyle="--")
				ax1.grid()

				plt.show()
			if w < WW:
				dels.append(TF)
		else:
			dels.append(TF)
	for d in dels:
		for chrom in G:
			for b in G[chrom]:
				if d in b.TFS:
					del b.TFS[d]
	print len(dict([(tf,1)for chrom in G for b in G[chrom] for tf in b.TFS]))

class cluster:
	def __init__(self,start, stop, bidirs):
		self.bidirs 	= bidirs
	def pick_closest(self):
		for b in self.bidirs:
			MIN,ARGMIN=None, None
			for TF in b.TFS:
				if MIN is None or abs(b.TFS[TF]) < MIN:
					MIN 	= b.TFS[TF]
					ARGMIN 	= TF
			dels 	= list()

			for TF in b.TFS:
				if abs(b.TFS[TF]) > 25:
					dels.append(TF)
			for d in dels:
				del b.TFS[d]
		pass
	def in_cluster(self, tf1, tf2):
		B 	= None
		for b in self.bidirs:
			if tf1 in b.TFS:
				B 	= b
				break
		if B is None:
			return 0
		S=0
		for b in self.bidirs:
			if tf2 in b.TFS and b!=B:
				S+=1
		return S

			
def linkage(G,threshold=10000):
	clusters 	= list()
	for chrom in G:
		LST 	= G[chrom]
		j,N 	= 0,len(LST)
		while j < N :
			start, stop = LST[j].start-threshold, LST[j].stop + threshold
			ct 	= 0
			bidirs=list()
			while j < N and start < LST[j].stop  and stop > LST[j].start :
				start,stop 	= min(start, LST[j].start -threshold), max(stop, LST[j].stop +threshold)
				bidirs.append(LST[j])
				j+=1
			clusters.append(cluster(start, stop, bidirs))
			clusters[-1].pick_closest()
			j+=1
	return clusters
def co_ocurrence_linkage(clusters):
	TFS 	= dict([(TF,0) for c in clusters for b in c.bidirs for TF in b.TFS])
	LABELS 	= {}
	for i,TF in enumerate(TFS.keys()):
		TFS[TF] 	= i
		LABELS[i] 	= TF
	A 	= np.zeros((len(TFS), len(TFS)))
	for c in clusters:

		for TFi in TFS:
			for TFj in TFS:
				S 	=c.in_cluster(TFi, TFj)
				i,j 	= TFS[TFi], TFS[TFj]
				A[i,j]+=S
				A[j,i]+=S

	G=nx.Graph()
	labels 	= {}
	colors=list()
	ALL 	= list()
	for i in range(A.shape[0]):
		for j in range(i+1, A.shape[0]):
			d 			= (A[i,j]  )  +  (A[j,i] )
			d/=2.
			ALL.append(d)
			if d> 0.03:
				labels[i] 	= LABELS[i]
				labels[j] 	= LABELS[j]
		
				G.add_node(i)
				G.add_node(j)
				
				G.add_edge(i,j,weight=d ,alpha=0.5,edge_color='r')
	for node in G.nodes():
		print node
		colors.append(sum(A[node,:]))
	F 	= plt.figure(figsize=(15,10))
	pos 	= nx.spring_layout(G)
	ax 	= F.add_axes([0.1,0.1,0.4,0.8])
	ax2 = F.add_axes([0.6,0.1,0.3,0.35])
	ax2.hist(ALL , bins=50 )
	ax3 = F.add_axes([0.6,0.5,0.3,0.35])
	NORMED 	= np.zeros(A.shape)
	for i in range(A.shape[0]):
		for j in range(A.shape[0]):
			NORMED[i,j] 	= A[i,j] 
	cmap = mpl.cm.cool
	norm = mpl.colors.Normalize(vmin=0, vmax=np.max(A))
	axcb = F.add_axes([0.95, 0.5, 0.025, 0.35])
	cb1 = mpl.colorbar.ColorbarBase(axcb, cmap=cmap,
                               norm=norm,
                               orientation='vertical')
	ax3.pcolor(NORMED,cmap=cmap, vmin=np.min(NORMED), vmax=np.max(NORMED))
	nx.draw_networkx_nodes(G, pos=pos, font_size=8, node_size=300,cmap=plt.get_cmap('cool'), node_color=colors, ax=ax)
	nx.draw_networkx_edges(G, pos=pos, color="r", alpha=0.1,width=2,ax=ax)
	nx.draw_networkx_labels(G,pos,labels,font_size=10, ax=ax)
	plt.show()			
		
def load_PSSMs(FILE):
	D 	= {}
	with open(FILE) as FH:
		C 	= False
		M,MOTIF 	= list(), None
		for line in FH:

			if "MOTIF" in line[:5]:
				if M and MOTIF is not None:
					D[MOTIF] = np.array(M)
				MOTIF 	= re.split("\s",line.strip("\n"))[1]
				M  =list()

				C 	= False
			if "letter" in line[:6]:
				C 	= True
			elif C:
				line_array 	= re.split("\s+", line[1:-1].strip("\n"))
				if line_array[0]:
					M.append([float(x) for x in line_array[:-1]] )
	return D
def cluster_PSSMS(D):
	DIMS 	= [m.shape[0] for m in D.values() ]
	for DIM in range(min(DIMS), max(DIMS)):
		X 	= list()
		IDS = list()
		for d in D:

			if D[d].shape[0]==DIM:
				IDS.append(d)
				X.append(D[d])
		if X:
			A 	= np.zeros((len(X), len(X)))
			for i in range(len(X)):
				for j in range(i+1, len(X)):
					A[i,j] 	= np.sum(abs(X[i] - X[j]))
					A[j,i] 	= A[i,j]
			linkage_matrix = linkage_scipy(A, "complete")
			
			ddata = dendrogram(linkage_matrix,
		 	                  labels=IDS , color_threshold=0,leaf_rotation= 85, leaf_font_size=14 )
			plt.tight_layout()
			plt.show()

def make_eRNA_cluster_FILE(IN, OUT):
	thresh=2
	G 	= {}
	with open(IN) as FH:
		for line in FH:
			chrom,Start, Stop, info, TFS 	= line.strip("\n").split("\t")
			#condense TFS
			T 	= list()
			if TFS:
				for pair in TFS.strip(",").split(","):
					tf, distance 	= pair.split(":")
					distance 		= float(distance)
					T.append((distance,tf))
				T.sort()
				j,N 	= 0,len(T)
				OS 		= list()
				while j < N:
					start, stop = T[j][0]-thresh, T[j][0]+thresh
					O 			= {}
					while j < N and stop > (T[j][0]-thresh) and start < (T[j][0]+thresh):
						if T[j][1] not in O:
							O[T[j][1]] 	= T[j][0]
						if abs(O[T[j][1]]) > abs(T[j][0]):
							O[T[j][1]] 	= T[j][0]
						start, stop 	= min(T[j][0]-thresh, start), max(T[j][0]+thresh, stop)
						j+=1
					OS.append(O)
					j+=1
			if chrom not in G:
				G[chrom]=list()
			G[chrom].append((int(Start), int(Stop), info, OS))
	FHW	= open(OUT, "w")
	for chrom in G:
		for start, stop, info, OS in G[chrom]:
			start 	= str(start)
			stop 	= str(stop)
			TFS 	= ":".join([ ",".join([tf for tf in O ]) + str(int(np.mean(O.values()))) for O in OS ])
			FHW.write(chrom+"\t"+start+"\t" + stop + "\t" + TFS + "\t"+info + "\n")
		
		G[chrom].sort()
	thresh=5000
	for chrom in G:
		g 		= G[chrom]
		g.sort()
		j,N 	= 0,len(G[chrom])

		while j < N:
			start, stop 	= g[j][0],g[j][1]
			ct 				= 0
			while j < N and stop > (g[j][0]-thresh) :
				stop 	= max(g[j][1],stop )
				start 	= min(start,g[j][0] )
				ct+=1
				j+=1
			if ct > 1:
				FHW.write(chrom+"\t" + str(start) + "\t" + str(stop)+  "\t" + "eRNA_cluster" + "\n")
	FHW.close()












if __name__ == "__main__":
	DIST 		= True
	load 		= False
	DENDRO 		= False
	eRNA_link	= False

	if eRNA_link:
		OUT 	= "/Users/joazofeifa/Lab/TF_predictions/eRNA_clusters.bed"
		IN 		= "/Users/joazofeifa/Lab/TF_predictions/assignments_refined/Allen2014_DMSO2_3-1_0.05"
		make_eRNA_cluster_FILE(IN,OUT)

	if DENDRO:
		FILE 	= "/Users/joazofeifa/Lab/TF_predictions/HOCOMOCOv9_AD_MEME.txt"
		D 		= load_PSSMs(FILE)
		cluster_PSSMS(D)
	if DIST:
		if len(sys.argv)<2:
			root 	= "/Users/joazofeifa/Lab/ENCODE/HCT116/"
			query 	= "/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/Allen2014_DMSO2_3-4_bidirectional_hits_intervals.bed"
			out 	= "/Users/joazofeifa/Desktop/motif_distances.tsv"
			pad 	= 2000
		else:
			root 	= sys.argv[1]
			query 	= sys.argv[2]
			out 	= sys.argv[3]
			pad 	= int(sys.argv[4])

		Q,IDS 	= make_query(query, pad=pad)
		A,M 	= read_in_directory(root, Q)
		write_out(A,M,IDS, out)
	if load:
		IN 	= "/Users/joazofeifa/Lab/TF_predictions/distances_crude/"
		OUT = "/Users/joazofeifa/Lab/TF_predictions/assignments_refined/"
		W 	= 0.05
		for IN_FILE in os.listdir(IN):
			print IN_FILE
			TSS_T 	= get_TSS_tree("/Users/joazofeifa/Lab/genome_files/TSS.bed")
		
			FILE 		= IN + IN_FILE
			G 			= load_TF_distance_out(FILE,test=False,FILTER=None)
			FILTER(G, WW=W)
			OUT_FILE 	=  OUT +  IN_FILE +  "_" + str(W)
			write_out_TF_classifications(G, OUT_FILE)
			#clusters 	= linkage(G)
			#co_ocurrence_linkage(clusters)
		#co_occurance_over_eRNAs(G)
		#look_at_occupancy(G)
