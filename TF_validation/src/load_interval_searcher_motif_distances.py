import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import matplotlib.cm as cm
import math
def load(FILE,  test=False):
	t 	= 0
	G 	= {}
	with open(FILE) as FH:
		for line in FH:
			t+=1
			if t > 1000 and test :
				break
			chrom,start, stop, info 	= line.strip("\n").split("\t")
			center 						= (float(stop) + float(start)) / 2.
			info, distances 			= info.split("|")
			distances 					= distances.split(",")
			D 							= [d.split(":") for d in distances]
			if len(D[0]) > 2:
				try:
					for st_sp, info, file_name in D:
						x 	= sum([float(i) for i in st_sp.split("-")]) / 2.
						d 	= center - x
						if info not in G:
							G[info]=list()
						G[info].append(d)
				except:
					print D, D[0]
					break
	return G	
def get_w(d, mu=0, s=100, bins=200):
	counts,edges = np.histogram(d, bins=bins)
	X = np.zeros((bins, 2))
	edges=(edges[1:]+ edges[:-1])/2.

	X[:,0]=edges
	X[:,1]=counts
	def normal(x, mu, si):
		return (1.0 / (math.sqrt(2*math.pi) * si ))* math.exp(-pow(x-mu,2)/(2*pow(si,2)))
	def uniform(x,a,b):
		return 1.0 / (b-a)
	w =	0.5
	t = 0
	a,b = min(X[:,0]), max(X[:,0])
	max_iterations=100
	while t < max_iterations:
		EN, EU = 0,0
		ll 	= 0
		for i in range(X.shape[0]):

			x,y = X[i,:]
			rn= w*normal(x, mu, s) / (w*normal(x, mu, s) + (1-w)*uniform(x, a,b))
			ru = (1-w)*uniform(x, a,b) / (w*normal(x, mu, s) + (1-w)*uniform(x, a,b))
			EN+=(rn*y)
			EU+=(ru*y)
			ll+=math.log(rn + ru)*y
		w = EN / (EU+ EN)
		t+=1
	return w


def display(G, specs):
	F 	= plt.figure()
	i 	= 0
	st 	= 0.1
	hspace 	= 0.1
	hgap 	= 0.01
	for spec in specs:
		motifs  	= [m for m in G.keys() if spec in m ]
		if not len(motifs):
			print "not found?", spec
		else:
			distances 		= list()
			for m in motifs:
				distances 	+= G[m]
			X,edges 	= np.histogram(distances, bins=100, range=(-100,100))
			norm    		= mpl.colors.Normalize(vmin=min(X), vmax=max(X) )
			cmap    		= cm.Blues
			m       		= cm.ScalarMappable(norm=norm, cmap=cmap)
			colors 			= [m.to_rgba(d) for d in X]
			ax 				= F.add_axes([st , i*hspace+hgap, 0.8  ,hspace*0.5  ])
			ax.set_xticks([])
			ax.set_yticks([])

			ax.bar(range(len(X)),np.ones( (len(X),)), color=colors, width=1.0    , edgecolor=colors )
			ax.spines['right'].set_visible(False)
			ax.spines['left'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.spines['top'].set_visible(False)			
			i+=1
	plt.show()

def load_distance_parsed(FILE):
	G 	= {}
	with open(FILE) as FH:
		for line in FH:
			motif, w, distances 	= line.strip("\n").split("\t")
			distances 				= [float(d) for d in distances.split(",")]
			w 						= float(w)
			G[motif] 				= (w, distances)
	return G

def display_matrix(G,bins=100, write_out=False):
	if write_out:
		FHW 	= open("/Users/joazofeifa/Lab/TF_predictions/distances_parsed_EM_w.tsv", "w")
		for i,m in enumerate(G.keys()):
			w 	= get_w(G[m])
			print i
			FHW.write(m+"\t" + str(w) +"\t" + ",".join([str(d) for d in G[m]]) + "\n" )
			FHW.flush()

			
	G 	= load_distance_parsed("/Users/joazofeifa/Lab/TF_predictions/distances_parsed_EM_w.tsv" )

	motifs 	= [(w,x) for x,(w, d) in zip(G.keys(), G.values() ) ]
	motifs.sort()


	N 	= len(motifs)
	X 	= np.zeros(( N,bins))
	labels 	= list()


	t 	= 0
	r 	= 5
	M 	= ""
	for i,(w,m) in enumerate(motifs):
		counts,edges 	= np.histogram( G[m][1] ,bins=bins, range=(-500,500))
		X[i,:] 	 		= counts
		X[i,:]/=float(sum(X[i,:]))
		if t > r:
			M 	= M.strip(",")
			labels.append(M)
			t=0
			M=""
		M+=m.split("_")[0]+","
		t+=1
	F 	= plt.figure(figsize=(15,10))
	ax 	= F.add_subplot(1,1,1)
	ax.set_xticks([])
	ax.set_yticks(np.linspace(0,X.shape[0], len(labels)))
	ax.set_yticklabels(labels, fontsize=8)
	heatmap = ax.pcolor(X, cmap=plt.cm.Blues, vmin=np.min(X)  , vmax=np.max(X)*0.9)
	plt.show()











if __name__ == "__main__":
	FILE 	= "/Users/joazofeifa/Lab/TF_predictions/interval_searcher_motif_distances/Allen2014_overlaps.bed"
	FILE 	= "/Users/joazofeifa/Lab/TF_predictions/interval_searcher_motif_distances/FIS_RUN_overlaps.bed"
	G 		= load(FILE)
	#G 		= None
	specs 	= ("ZBTB","ELF1", "EGR1", "MAZ", "YY1", "MAX", "TEAD4", "ATF3", "SRF", "SP1", "USF1", "JUND", "CEBPB", "FOX", "FOSL1", "REST", "CTCF", "HEY" )
	display(G, specs)
	#display_matrix(G)

