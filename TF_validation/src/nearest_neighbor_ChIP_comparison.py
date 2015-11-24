import sys,os
sys.path.append("/Users/joazofeifa/Lab/EMG/TF_predictions/src/")
sys.path.append("/Users/joazofeifa/Lab/interval_searcher/")
import node
from ROC_from_k_model_fits import segment
import matplotlib.pyplot as plt
import numpy as np
import time
import math
class comparison:

	def __init__(self, cell_type_1, cell_type_2, TF):
		self.cell_type_1, self.cell_type_2 	= cell_type_1, cell_type_2
		self.TF 							= TF
		self.cell_type_1_model_fits 		= {}
		self.cell_type_2_model_fits 		= {}
		self.unique_to_cell_type_1 			= {}
		self.unique_to_cell_type_2 			= {}
		self.intersect_of_both_ct_1 		= {}
		self.intersect_of_both_ct_2 		= {}
	def make_searchable(self):
		GG 	= list()
		SS 	= (self.unique_to_cell_type_1,self.unique_to_cell_type_2)

		for S in SS:
			G  	= {}
			for i in S:
				s 	= S[i]
				if s.chrom not in G:
					G[s.chrom]=list()
				G[s.chrom].append((s.start, s.stop, s))
			for chrom in G:
				G[chrom].sort()
				G[chrom] 	= node.tree(G[chrom])
			GG.append(G)
		self.unique_to_cell_type_1_searchable 	= GG[0]
		self.unique_to_cell_type_2_searchable 	= GG[1]
	def search(self, cell_type_1, cell_type_2, chrom,start, stop):
		a,b,c,d 	= cell_type_1, cell_type_2,  self.cell_type_1, self.cell_type_2

	
		if (a == c or a == d) and (b == c or b == d) and chrom in self.unique_to_cell_type_1_searchable:
			return self.unique_to_cell_type_1_searchable[chrom].searchInterval((start, stop))
		return list()



	def insert_model_fit(self, FILE, cell_type_1=False):
		G,S = {},None
		with open(FILE) as FH:
			for line in FH:
				if line[0]==">":
					if S is not None :
						G[ID]=S
					S 	=  segment(line)
					ID 	= line.split("|")[0]
				elif S is not None:
					S.add_fit(line)
		if cell_type_1 == True:
			self.cell_type_1_model_fits 	= G
		else:
			self.cell_type_2_model_fits 	= G

	def get_bidir_counts(self, penality=500):
		cell_type_1_s 	= set([ID for ID in self.cell_type_1_model_fits if self.cell_type_1_model_fits[ID].get_BIC_mdoel(penality)  ])
		cell_type_2_s 	= set([ID for ID in self.cell_type_2_model_fits if self.cell_type_2_model_fits[ID].get_BIC_mdoel(penality)  ])
		
		N 				= len(self.cell_type_1_model_fits)
		I 				= cell_type_1_s.intersection(cell_type_2_s)
		
		self.intersect_of_both_ct_1 	= dict([ (i, self.cell_type_1_model_fits[i])  for i in I])
		self.intersect_of_both_ct_2 	= dict([ (i, self.cell_type_2_model_fits[i])  for i in I])
		self.unique_to_cell_type_1  	= dict([(i, self.cell_type_1_model_fits[i]) for i in cell_type_1_s.difference(cell_type_2_s) ])
		self.unique_to_cell_type_2  	= dict([(i, self.cell_type_2_model_fits[i]) for i in cell_type_2_s.difference(cell_type_1_s) ])
		return len(self.unique_to_cell_type_1), len(self.unique_to_cell_type_2), len(self.intersect_of_both_ct_1), N


			
def load_all(DIRS, test=False, labels=()):
	T 	= dict(zip(labels, DIRS))
	G 	= {}
	CS 	= list()
	SORT = dict([(l,i) for i,l in enumerate(labels)])
	for cell_type, DIR in zip(T.keys(), T.values()):
		for TF_OUT 	in os.listdir(DIR):
			TF_OUT 	= TF_OUT.split("-")[0]
			if TF_OUT not in G:
				G[TF_OUT]=list()
			G[TF_OUT].append(cell_type)
	for i,TF in enumerate(G):
		cell_type_1, cell_type_2 			= G[TF]
		cell_type_1_dir, cell_type_2_dir	= T[cell_type_1], T[cell_type_2]
		cell_type_1_file  					= [cell_type_1_dir + F for F in os.listdir(cell_type_1_dir) if TF in F ][0]
		cell_type_2_file  					= [cell_type_2_dir + F for F in os.listdir(cell_type_2_dir) if TF in F ][0]
		c 									= comparison(cell_type_1, cell_type_2, TF)
		c.insert_model_fit(cell_type_1_file, cell_type_1=True)
		c.insert_model_fit(cell_type_2_file, cell_type_1=False)
		CS.append((SORT[cell_type_1], c))

		if test and i > 10:
			break
	CS.sort()
	CS 	= [y for x,y in CS]
	for C in CS:
		C.get_bidir_counts()
		C.make_searchable()
	return CS
def graph_percents(CS):
	half 	= len(CS) / 2
	F 			= plt.figure(figsize=(15,10))
	for ii in (1,2):
		if ii == 1:
			st ,	sp 	= 0,half
			ax 			= F.add_subplot(1,2,1)
		else:
			st 	, sp 	= half, len(CS)
			ax 			= F.add_subplot(1,2,2)
		print st, sp
			


		ct1_unique 	= list()
		ct2_unique 	= list()
		ct1_unique_counts 	= list()
		ct2_unique_counts 	= list()
		
		intersect 	= list()
		intersect_counts = list()
		all_counts 	= list()
		labels 		= list()
		positions 	= list()
		PS 			= list()
		delta 		= 0.4
		height 		= 0.4
		gap 		= height*0.1
		colors 		= {"HCT116":"g", "K562":"b", "A549": "r", "HeLa": "y"}
		ct1_colors 	= list()
		ct2_colors 	= list()
		
		for i,C in enumerate(CS[st:sp]):
			ct1_u, ct2_u, I, N= C.get_bidir_counts()

			ct1_unique.append(float(ct1_u)/float(N))
			ct2_unique.append(float(ct2_u)/float(N))
			ct1_unique_counts.append(ct1_u)	
			ct2_unique_counts.append(ct2_u)	
			intersect_counts.append(I)	
			all_counts.append(N)
			ct1_colors.append(colors[C.TF.split("_")[0]])
			ct2_colors.append(colors[C.TF.split("_")[1]])
			
			intersect.append(float(I)/N )
			labels.append(C.TF)
			positions.append(i*(delta))

		positions 	= np.array(positions)
		ct1_unique, ct2_unique 	= np.array(ct1_unique), np.array(ct2_unique)
		intersect 	= np.array(intersect)
		ax.barh(positions - gap  , intersect, height=height/4., color=ct1_colors , hatch="//" ,edgecolor="w" )
		ax.barh(positions + gap*2 , intersect, height=height/4., color=ct2_colors , hatch="//" ,edgecolor="w")
		
		rects1 		= ax.barh(positions - gap  , ct1_unique, height=height/4., color=ct1_colors, left=intersect  ,edgecolor="w")
		rects2 		=ax.barh(positions + gap*2 ,ct2_unique, height=height/4., color=ct2_colors, left=intersect  ,edgecolor="w")
		ax.set_yticks(positions +(height )*0.25 )
		ax.set_yticklabels( [tf.split("_")[-1] +", N=" + str(all_counts[i]) for i,tf in enumerate(labels)]  )
		if (ii==1):
			L 	= {}
			for i,r in enumerate(rects1):
				if r._facecolor not in L:
					L[r._facecolor] = (r,labels[i].split("_")[0])
			for i,r in enumerate(rects2):
				if r._facecolor not in L:
					L[r._facecolor] = (r,labels[i].split("_")[1])

			ax.legend([L[color][0] for color in L ],[L[color][1] for color in L ], loc=(0.9,0.85))
		for i,TF in enumerate(labels):
			per 	= intersect[i]
			per2 	= intersect[i]+ct1_unique[1]

			pos 	= positions[i] - gap + (height/8.)
			# ax.text(per  , pos    , TF.split("_")[0]+ "    ",color="white" ,fontsize=10,
			# 			verticalalignment='center', horizontalalignment="right" )
			
			pos 	= positions[i] + gap*2 + (height/8.)
			# ax.text(per  , pos    , TF.split("_")[1]+ "    ",color="white" ,fontsize=10 ,
			# 	verticalalignment='center', horizontalalignment="right" )



		ax.grid()
		ax.set_xlim(0,1.2)
		ax.set_ylim(min(positions)-delta, max(positions)+delta)
	plt.tight_layout()
	plt.show()

def load_gene_counts(FILE ):
	G 	= {}
	with open(FILE) as FH:
		for line in FH:
			name,chrom,start, stop, cov 	= line.strip("\n").split("\t")
			cov 							= sum([abs(float(x)) for x in cov.split(",")])
			if chrom not in G:
				G[chrom]= 	list()
			G[chrom].append((int(start) , int(stop) , (name, cov)))
	for chrom in G:
		G[chrom].sort()
		G[chrom]= node.tree(G[chrom])
	return G
def load_gene_count_files(FILES, labels):
	G 	= {}
	for F, label in zip(FILES, labels):
		G[label] 	= load_gene_counts(F)
	return G

def format_bp(bp):
	colors 	= ("green", "red", "blue")
	for i,box in enumerate(bp['boxes']):
	    # change outline color
	    box.set( color='#7570b3', linewidth=2)
	    # change fill color
	    box.set( facecolor = colors[i] )

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
	    flier.set(marker='o', color='#e7298a', alpha=0.)
def nearest(genes, CS, W = pow(10,5)):
	L 	= list()
	for i,C in enumerate(CS ):
		C.get_bidir_counts()
		ct1, ct2 			= C.cell_type_1,C.cell_type_2
		if ct1 in genes and ct2 in genes:
			#want to compare the gene expression of cell type 1 and cell type 2 
			#in situataions where its unique to either ct1 or ct2 (we expect these to be different) and in both ct1 and ct2 (null kinda)
			HS 	= list()
			for H in (C.intersect_of_both_ct_1 , C.unique_to_cell_type_1, C.unique_to_cell_type_2):
				A,B 	= list(),list()
				for s in H:
					ct_1 				= H[s]

					chrom, start, stop 	= ct_1.chrom, ct_1.start, ct_1.stop
					FINDS 				= genes[ct1][chrom].searchInterval((start, stop))
					x 					= (float(stop) + float(start))/2.
					FINDS_ct1 			= genes[ct1][chrom].searchInterval((start-W, stop+W))
					FINDS_ct2 			= genes[ct2][chrom].searchInterval((start - W, stop + W))
					if FINDS_ct1 and FINDS_ct2:
						info_ct1 			=  min([ (abs( x - ((st + sp)/2.)), name, cov) for st, sp, (name, cov) in FINDS_ct1]) 
						info_ct2 			=  min([ (abs( x - ((st + sp)/2.)), name, cov) for st, sp, (name, cov) in FINDS_ct2]) 
						A.append(math.log(info_ct1[2]+1,10))
						B.append(math.log(info_ct2[2]+1,10))
				HS.append((A, B))

			L.append((C.cell_type_1, C.cell_type_2, C.TF, HS   ))
	delta 	= len(L)/4
	t 		= 0
	z 		= 0
	F 		= plt.figure(figsize=(15,10))
	ax 		= F.add_subplot(2,2,1)
	ax.grid()
	ax.set_ylim(-2,4)
	ticks 	= [-1]
	labels 	= [""]
	for i,(ct1, ct2, tf, HS) in enumerate(L):
		if t > delta:
			ax.set_xticks(ticks)
			ax.set_xticklabels(labels,fontsize=10)
			ticks 	= [-1]
			labels 	= [""]
			z+=1
			ax 	= F.add_subplot(2,2,z+1)
			ax.grid()
			ax.set_ylim(-2,4)
			t=0
		DIFF_INTERSECT 	= [a-b for a,b in zip(HS[0][0],HS[0][1])]
		unique_ct1 		= [a-b for a,b in zip(HS[1][0],HS[1][1])]
		unique_ct2 		= [a-b for a,b in zip(HS[2][0],HS[2][1])]
		
		bp = ax.boxplot((DIFF_INTERSECT, unique_ct1, unique_ct2),patch_artist=True, positions=((t*2) ,(t*2)+ 0.5,(t*2) + 1  ))
		ticks.append((t*2)+ 0.5)
		labels.append(ct1+"\n" + ct2 + "\n" + tf.split("_")[-1])
		format_bp(bp)
		t+=1
	ax.set_xticks(ticks)
	ax.set_xticklabels(labels,fontsize=10)
	plt.tight_layout()
	plt.show()








if __name__ == "__main__":
	SHOW_PERCENT 	= True
	#============================================================================================
	#TF_BIDIR Input Files
	HCT116 			= "/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_ChIP_comparisons/"
	K562 			= "/Users/joazofeifa/Lab/gro_seq_files/Core2014/EMG_out_ChIP_comparisons/"
	HeLa 			= "/Users/joazofeifa/Lab/gro_seq_files/Andersson2014/EMG_out_ChIP_comparisons/"
	A549 			= "/Users/joazofeifa/Lab/gro_seq_files/Luo2014/EMG_out_ChIP_comparisons/"
	#============================================================================================
	#Gene Count Files
	genes_HCT116 	= "/Users/joazofeifa/Lab/genome_files/DMSO2_3_Gene_Counts.tsv"
	genes_K562 		= "//Users/joazofeifa/Lab/genome_files/Core2014Gene_Counts.tsv"
	genes_HeLa 		= "/Users/joazofeifa/Lab/genome_files/Andersson2014Gene_Counts.tsv"
	genes_A549 		= "/Users/joazofeifa/Lab/genome_files/Luo2014Gene_Counts.tsv"

	
	DIRS 			= (HCT116, K562, HeLa, A549)
	#DIRS 			= (HCT116, HeLa, A549)
	genes_DIR 	 	= (genes_HCT116,genes_K562,genes_A549, genes_HeLa )
	#genes_DIR 	 	= (genes_HCT116, genes_A549, genes_HeLa )
	labels 			= ("HCT116", "K562", "HeLa", "A549" )
	CS 				= load_all(DIRS, labels=labels, test=False )
	if SHOW_PERCENT:
		graph_percents(CS)


















