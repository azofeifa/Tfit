import time
import matplotlib.pyplot as plt
import numpy as np
import math,os
from scipy.stats import gaussian_kde,pearsonr
from scipy import stats
from pylab import setp
class bidir:
	def __init__(self,chrom,start,stop, info, TFS):
		self.chrom,self.start, self.stop= chrom,float(start), float(stop)
		self.mt 						= 1
		self.TFS 						= {}
		self.collapse 					= {"SP1_f1": "Zinc",
											"WT1_f1" : "Zinc",
											"SP4_f1" : "Zinc",
											"MAZ_f1" : "Zinc",			
											"ZN148_si":"Zinc",
											"ZBT7B" : "Zinc"
											}
		self._parse_info(info)
		self._parse_TFS(TFS)
		self.DNAse 	= False

	def _parse_info(self, info):
		line_array 	= [float(x) for x in info.split("_")]
		self.mu,self.si, self.l, self.w, self.pi,self.N,self.fp, self.ll = line_array
	def _parse_TFS(self, TFS):
		INIT_TFS =[ (float(pair.split(":")[1]),  pair.split(":")[0]  ) for pair in TFS.split(",") if pair ]
		INIT_TFS.sort()
		#concantentate
		i, N 	= 0, len(INIT_TFS)
		while i < N:
			start, stop 	= INIT_TFS[i][0]-self.mt, INIT_TFS[i][0]+ self.mt
			IDS 			= ""
			old=i
			while i < N and INIT_TFS[i][0]+self.mt > start and INIT_TFS[i][0]-self.mt < stop:
				start, stop = min(start, INIT_TFS[i][0]-self.mt), max(stop,INIT_TFS[i][0]+self.mt )	
				IDS+=INIT_TFS[i][1] + ","
				i+=1
			IDS 	= IDS.strip(",")
			NIDS 	= ""
			for ID in IDS.split(","):
				if ID in self.collapse:
					NIDS+=self.collapse[ID]+","
				else:
					NIDS+=ID+","
			NIDS=NIDS.strip(",")
			i+=1
			self.TFS[IDS] 	= (stop+start) /2.
	def get_density(self):
		return float(self.N+1) /  (self.stop - self.start)

def load(FILE,test=True):
	G 	= {}
	t 	= 0
	with open(FILE) as FH:
		for line in FH:
			line_array 	= line.strip("\n").split("\t")
			chrom,start, stop,info 	= line_array[:4]
			if test and t > 1000:
				break
			t+=1
			if len(line_array) > 4:
				TFS 	= line_array[4]
			else:
				TFS 	= ""
			if chrom not in G:
				G[chrom]= list()
			b 	= bidir(chrom,start, stop, info, TFS)
			G[chrom].append((int(start), b))
	#want to collapse TF clusters...
	for chrom in G:
		G[chrom].sort()
		G[chrom] 	= [y for x,y in G[chrom]]
	return G
def unique_motifs_per_bidir(GS,FILES):

	D 	= [[ len(b.TFS.keys()) for chrom in G for b in G[chrom]] for G in GS   ] 
	F 	= plt.figure(figsize=(15,10))
	
	ax1 = F.add_subplot(1,2,1)
	bp 	= ax1.boxplot(D,patch_artist=True)
	ax1.set_xticklabels(FILES,rotation=90)
	ax1.set_ylabel("Number of Unique Motifs\nPer Bidirectional")
	for box in bp['boxes']:
		# change outline color
		box.set( color='#7570b3', linewidth=2)
		# change fill color
		box.set( facecolor = '#1b9e77' )

	## change color and linewidth of the whiskers
	for whisker in bp['whiskers']:
		whisker.set(color='#7570b3', linewidth=2)

	## change color and linewidth of the caps
	for cap in bp['caps']:
		cap.set(color='#7570b3', linewidth=2)

	## change color and linewidth of the medians
	for median in bp['medians']:
		median.set(color='#b2df8a', linewidth=2)

	## change the style of fliers and their fill
	for flier in bp['fliers']:
		flier.set(marker='o', color='#e7298a', alpha=0.5)
	ax1.grid()
	ax1.set_ylim((-1,13))
	xys 	= [  [ (len(b.TFS.keys()),math.log(b.get_density(),10))   for chrom in G for b in G[chrom] ] for G in GS  ]
	xys 	= [(pearsonr([x for x,y in d],[y for x,y in d])[0],d,i) for i,d in enumerate(xys)]
	ro,xy,i 	= max(xys)
	
	ax2 = F.add_subplot(1,2,2)
	x,y = [x for x,y in xy],[y for x,y in xy]
	ro, pv 	= pearsonr(x,y)
	xy = np.vstack([x,y])
	
	z = gaussian_kde(xy)(xy)
	ax2.set_title("Dataset: " + str(FILES[i]))
	ax2.scatter(x, y, c=z, s=16, edgecolor='', label=r'$\rho$' + ": " + str(ro)[:4] + "\n" + r'$p_{val}$' + ": " + str(pv) )
	ax2.grid()
	ax2.set_xlabel("Number of Motifs")
	ax2.set_ylabel("Log of Read Density")
	ax2.legend()
	plt.tight_layout()

	plt.show()

def load_rank_list(FILE):
	G 	= {}
	with open(FILE) as FH:
		for line in FH:
			TF, rank 	= line.strip("\n").split(",")
			G[TF] 	= float(rank)
	return G
def setBoxColors(bp):
	for i,box in enumerate(bp['boxes']):
	    # change outline color
		if i == 0:
			box.set( color='black', linewidth=1)
			# change fill color
			box.set( facecolor = 'r' )
		else:
			box.set( color='black', linewidth=1)
			# change fill color
			box.set( facecolor = 'b' )
		

	## change color and linewidth of the whiskers
	for whisker in bp['whiskers']:
	    whisker.set(color='#7570b3', linewidth=2)

	## change color and linewidth of the caps
	for cap in bp['caps']:
	    cap.set(color='#7570b3', linewidth=2)

	## change color and linewidth of the medians
	for median in bp['medians']:
	    median.set(color='#b2df8a', linewidth=2)

	## change the style of fliers and their fill
	for flier in bp['fliers']:
	    flier.set(marker='o', color='#e7298a', alpha=0.5)

def draw_box_plots_motif_influence(ranks, T,ax,FILTER1=0, FILTER2=pow(10,-10), title="" ):
	
	delta 	= len(ranks)
	st 	= 0
	bs 	= list()
	SORTED 	= [ (np.mean([d for tf in T[TF] for d in T[TF][tf][0] if d > -1 ]) - np.mean( [d for tf in T[TF] for d in T[TF][tf][1] if d > -1] ),TF)  for TF in ranks ]
	#SORTED.sort()
	final_TF= [TF for d, TF in SORTED]
	pvals 	= list()
	Final_TF=list()
	diffs 	= list()
	ONN,OFFF = list(),list()
	II 	= 0
	for i,TF in enumerate(final_TF):
		OFF 	= [d for tf in T[TF] for d in T[TF][tf][0] if d > -1 ]
		ON 		= [d for tf in T[TF] for d in T[TF][tf][1] if d > -1]
		ONN+=ON
		OFFF+=OFF
		if II < 12:
			if ON and OFF:
				diffs.append(np.mean(ON) - np.mean(OFF))
				pval 	= stats.ttest_ind(OFF,ON)[1]
				if FILTER1  <= pval <= FILTER2:
					II+=1

					if pval:
						mag 	= math.log(pval)
					else:
						mag 	= -900
					pvals.append(mag)
					
					bp=ax.boxplot([OFF, ON], 
						positions = [st, st+1], widths=0.7,patch_artist=True)
					bs.append(bp)
					setBoxColors(bp)
					Final_TF.append(TF)
					st+=4

			else:
				pval 	= 1

		
	ax.set_xlim(-7,st+8)
	ax.set_ylim(-1,13.5)
	positions 	= np.arange(-3.5, st+8, 4)
	print positions
	h 	= 11
	# for i,p in enumerate(positions):
	# 	if h == 12:
	# 		h=12.5
	# 	elif h==12.5:
	# 		h=13
	# 	else:
	# 		h=12
	# 	ax.text(p-0.5,h, r'$p<10^{' + str(pvals[i])[:6] + r'}$' )
	bp 	= ax.boxplot([OFFF,ONN], positions=[st, st+1],widths=0.6, patch_artist=True)
	setBoxColors(bp)

	ax.set_title(title)
	ax.set_xticks(positions)
	ax.set_xticklabels([" "]+[tf.split("_")[0]  for tf in Final_TF]+ ["combined"], rotation=45)
	hB, = plt.plot([-300,-300],'b-',linewidth=3.)
	hR, = plt.plot([-300,-300],'r-',linewidth=3.)
	ax.legend((hB, hR),('Motif Presence', 'Motif Abscence') )
	ax.grid()
	return 	ONN,OFF


def compare_exp(G, Ranks, display=True):
	F 	= plt.figure(figsize=(15,10))
	thresholds 	= [(0,0.5), (1.0,1.0)]
	FILTERS 	= [(0,pow(10,-10)), (pow(10,-10), 1)]
	Diffs 		= list()
	TITLES 		= ("Top Ten Nodes by Centrality", "Bottom Ten nodes by Centrality")
	for I,(start, stop) in enumerate(thresholds):
		ranks 		= [(Ranks[tf], tf ) for tf in Ranks if start<= Ranks[tf] <= stop and "E2F" not in tf and "NFIC" not in tf ]
		FILTER 		= dict([(tf, 1 ) for tf in Ranks if  Ranks[tf] ==1.0 ])
		ranks.sort()
		ranks 		= dict([(tf, d ) for d,tf in ranks   ])
		T 			= {}
		BS 			= [b for chrom in G for b in G[chrom]]
		for TF in ranks:
			if TF not in T:
				T[TF]=dict([(tf, [list(), list()] ) for chrom in G for b in G[chrom] for TF_cluster in b.TFS for tf in TF_cluster.split(",")] )
			for b in BS:
				clusters_N 	= len(b.TFS)
				clusters 	= b.TFS.keys()
				for i in range(clusters_N):
					others 	= dict([(tf,0) for j in range(clusters_N) for tf in clusters[j].split(",") if j!=i ])
					for tf in clusters[i].split(","):
						if tf in FILTER:
							if TF not in others:
								T[TF][tf][0].append(math.log(b.get_density()))
							else:
								T[TF][tf][1].append(math.log(b.get_density()))
		if I == 1:
			ax 	= F.add_axes([0.1,0.1,0.8,0.35])
		else:
			ax 	= F.add_axes([0.1,0.6,0.8,0.35])
			
		ON,OFF=draw_box_plots_motif_influence(ranks , T,ax,FILTER1=FILTERS[I][0],FILTER2=FILTERS[I][1], title=TITLES[I])
		Diffs.append((ON,OFF))
	# ax 	=  F.add_axes([0.55,0.1,0.4,0.85])
	# bp1 = ax.boxplot([Diffs[0][1],Diffs[0][0]],positions=[1,2],patch_artist=True )
	# bp2 = ax.boxplot([Diffs[1][1],Diffs[1][0]],positions=[3,4],patch_artist=True )
	# setBoxColors(bp1)
	# setBoxColors(bp2)
	# ax.set_xticks([1.5, 3.5])
	# ax.set_title("Combined High and Low Centrality")
	# ax.set_xticklabels(["High Centrality", "Low Centrality"], rotation=45)
	# hB, = plt.plot([-300,-300],'b-',linewidth=3.)
	# hR, = plt.plot([-300,-300],'r-',linewidth=3.)
	# ax.legend((hB, hR),('Motif Presence', 'Motif Abscence') )
	# ax.grid()
	# ax.set_ylim(-1,13.5)
	
	# ax.set_xlim(0,5)
	plt.show()

def distribution_of_GRO_signal_over_bdidireciontals(G):
	F 	= plt.figure(figsize=(15,10))
	ax 	= F.add_subplot(1,1,1)
	ax.grid()
	D 	= [math.log(b.get_density()) for chrom in G for b in G[chrom]]
	ax.hist(D, bins=200, normed=1, alpha=0.5)
	xs 	= np.linspace(min(D), max(D), 1000)
	mu 	= np.mean(D)
	si 	= np.std(D)
	f 	= lambda x: 1. / (math.sqrt(2*math.pi)*si)*math.exp(-pow(x-mu,2)/(2*pow(si,2)))
	ys 	= map(f, xs)
	ax.plot(xs,ys,linewidth=3., linestyle="--", color="black")
	plt.show()
	pass





if __name__ == "__main__":
	UNIQUE 	= False
	PIONEER = False
	GRO_strength 	= False
	if GRO_strength:
		IN_FILE 	= "/Users/joazofeifa/Lab/TF_predictions/assignments_refined/Allen2014_DMSO2_3-1_0.05"
		G 			= load(IN_FILE, test=True)
		distribution_of_GRO_signal_over_bdidireciontals(G)

	if PIONEER:
		RANK_FILE 	= "/Users/joazofeifa/Lab/EMG/TF_predictions/TF_Ranks.csv"
		IN_FILE 	= "/Users/joazofeifa/Lab/TF_predictions/assignments_refined/SRR1552480_Core2014-4_0.05"
		ranks 		= load_rank_list(RANK_FILE)
		G 			= load(IN_FILE, test=True)
		compare_exp(G,ranks)

	if UNIQUE:
		GS 		= list()
		IN 		= "/Users/joazofeifa/Lab/TF_predictions/assignments_refined/"
		thresh 	= "0.05"
		t 		= 0
		DS 		= list()
		for FILE in os.listdir(IN):
			IN_FILE = IN + FILE
			if IN_FILE[-4:]==thresh and "Puc" not in IN_FILE:
				GS.append(G)
				if "Allen" in FILE:
					DS.append(FILE.split("_")[0])
					
				elif "Sap" not in FILE:
					DS.append(FILE.split("_")[1])
				else:
					DS.append(FILE.split("_")[3])
					
				
		unique_motifs_per_bidir(GS,DS)