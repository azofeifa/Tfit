import sys,os
sys.path.append("/Users/joazofeifa/Lab/EMG/TF_predictions/src/")
sys.path.append("/Users/joazofeifa/Lab/interval_searcher/")
import node
from ROC_from_k_model_fits import segment
import matplotlib.pyplot as plt
import numpy as np
def load_k_model_fits(FILE):
	S 	= None
	G 	= list()
	with open(FILE) as FH:
		for line in FH:
			if line[0]==">":
				if S is not None :
					G.append(S)
				S 	=  segment(line)
			elif S is not None:
				S.add_fit(line)
	if S is not None and S.check():
		G.append(S)
	return G
def load_directory(DIR, test=False):
	GS,t 	= list(),0
	for FILE in os.listdir(DIR):
		if ".tsv" == FILE[-4:]:
			print FILE
			G 	= load_k_model_fits(DIR+"/"+FILE)
			if test and t > 4:
				break
			GS.append((FILE.split("-")[0], G))
			t+=1	
	return GS
def load_meta_data(FILE):
	with open(FILE) as FH:
		header	= True
		M 		= {}
		for line in FH:
			if not header:
				line_array 	= line.split("\t")
				ID,TF 		= line_array[0], line_array[15].split("-")[0]
				M[ID] 		= TF
			else:
				header=False
	return M
def make_DNAse_searchable(FILE):
	G 	= {}
	with open(FILE) as FH:
		for line in FH:
			chrom,start,stop 	= line.split("\t")[:3]
			if chrom not in G:
				G[chrom]=list()
			G[chrom].append((int(start)-2000, int(stop)+2000))
	for chrom in G:
		G[chrom].sort()
		G[chrom] 	= node.tree(G[chrom])
	return G




def graph_percent(GS, M, D, H, m1, REF, penality=150, SHOW=False, OUT=""):
	R 	= {}
	U 	= {}
	Z 	= {}
	DN 	= {}
	HB 	= {}
	HD 	= {}
	for F, G in GS:
		if "POL" not in M[F]:
			G 			= [s for s in G if s.chrom in D and s.chrom in H and s.chrom in m1 ]
			TF 			= M[F]
			U[TF] 		= len(G)
			N 			= float(len(G))
			B 			= [s for s in G if s.get_BIC_mdoel(penality) ]
			DNN 		= [s for s in G if D[s.chrom].searchInterval((s.start, s.stop)) ]
			H_BIDIR 	= [s for s in B if H[s.chrom].searchInterval((s.start, s.stop))   ]
			H_DNN 		= [s for s in DNN if H[s.chrom].searchInterval((s.start, s.stop))  ]
			S 			= float(len(B))
			DNN 	 	= float(len(DNN))

			if TF not in R:
				R[TF] 	= (N, S/N)
				DN[TF] 	= DNN / N
				HB[TF] 	= float(len(H_BIDIR)) / N
				HD[TF] 	= float(len(H_DNN)) / N
			else:
				if R[TF][0] < N:
					R[TF] 	= (N, S/N)
					DN[TF] 	= DNN / N
					HB[TF] 	= float(len(H_BIDIR)) / N
					HD[TF] 	= float(len(H_DNN)) / N
	for TF in R:
		R[TF]  = [R[TF][1]]
	positions,percents,percents_2, colors,labels 	= list(),list(),list(),list(),list()

	percents_DNase,percents_DNase_2  	= list(),list()
	delta,height 		= 0.8,0.6
		
	xy 	= [(np.mean(R[r]), r, DN[r], HD[r], HB[r])  for r in R]
	xy.sort()
	for i,(mean,TF, mean_DNAse, mean_HD, mean_HB) in enumerate(xy):
		positions+=[i*(delta*1.7) ]
		percents_DNase+=[mean_HD]
		percents_DNase_2+=[mean_DNAse-mean_HD]
		percents+=[mean_HB]
		percents_2+=[mean-mean_HB]

		colors.append(i)
		labels.append(TF + "\n(N=" + str(U[TF]) + ")")

	


	if SHOW:
		F 	= plt.figure(figsize=(5,10))
		ax 	= F.add_subplot(1,1,1)
		
		ax.barh(positions, percents, height=height/2., color="b", edgecolor='white',   label=" " )
		ax.barh(positions, percents_2,height=height/2., color="b", left=percents,edgecolor='white', label=" ", alpha=0.5 )

		ax.barh([p+(delta/2.) for p in positions], percents_DNase, height=height/2., color="r", edgecolor='white' , label=" "  )
		ax.barh([p+(delta/2.) for p in positions], percents_DNase_2, height=height/2., color="r",edgecolor='white', left =percents_DNase, label=" " , alpha=0.5 )
		
		ax.set_xlim(0,1.)
		ax.set_ylim(-delta ,max(positions)+delta*2)
		positions=np.array(positions) +(( height) / 2.)

		ax.set_xticks(np.arange(0,1.1,0.1))
		ax.set_xticklabels([str(int(x*100)) for x in np.arange(0,1.1,0.1)])
		ax.set_xlabel("Percent of TF Binding Events with\nDivergent Transcription Prediction")

		ax.set_yticks(positions)
		ax.set_yticklabels(labels)

		for i,(pos, per, per2 , perDNA, perDNA2) in enumerate(zip(positions, percents, percents_2 , percents_DNase, percents_DNase_2)):
			ax.text(per  , pos -  (height/2 )   ,   str(int(100*per/(per2+per))) + "%"  ,color="black" ,fontsize=12,
				verticalalignment='top', horizontalalignment="center" )
			ax.text(perDNA  , pos  + (height/2. ) + (delta/5.)  ,     str(int(100*perDNA/(perDNA2+perDNA)  ))+ "%"  ,color="black" ,fontsize=12 ,
				verticalalignment='bottom', horizontalalignment="center" )

		plt.tight_layout()
		plt.savefig("/Users/joazofeifa/Lab/Article_drafts/EMG_paper/images/DNAse_Bidir_H3K27ac_overlap.pdf")
		plt.show()







if __name__ == "__main__":
	Allen_DIR 		= "/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_ChIP/"
	Allen_META 		= "/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_ChIP/metadata.txt"
	HCT116_DNAse 	= "/Users/joazofeifa/Lab/ChIP/HCT116/DNase/wgEncodeUwDnaseHct116PkRep1.narrowPeak.bed"
	HCT116_H3K27ac 	= "/Users/joazofeifa/Lab/ChIP/HCT116/H3K27ac/HCT-116_H3K27Ac_narrowPeak.bed"
	HCT116_H3K4me1 	= "/Users/joazofeifa/Lab/ChIP/HCT116/H3K4me1/HCT-116_H3K4me1_narrowPeak.bed"
	TSS 			= "/Users/joazofeifa/Lab/genome_files/TSS.bed"
	K562 			= "/Users/joazofeifa/Lab/ChIP/K562/K562_SAHADukeDNaseSeq.pk"
	Core_DIR 		= "/Users/joazofeifa/Lab/gro_seq_files/Core2014/EMG_out_ChIP/"
	Core_META 		= "/Users/joazofeifa/Lab/gro_seq_files/Core2014/EMG_out_ChIP/metadata.txt"
	D 				= make_DNAse_searchable(HCT116_DNAse)
	H 				= make_DNAse_searchable(HCT116_H3K27ac)
	m1 				= make_DNAse_searchable(HCT116_H3K4me1)
	REF 			= make_DNAse_searchable(TSS)
	GS 				= load_directory(Allen_DIR, test=False)
	M 				= load_meta_data(Allen_META)
	G 				= graph_percent(GS, M, D, H, m1, REF, SHOW=True)

