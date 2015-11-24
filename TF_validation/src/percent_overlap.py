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


def graph_percent(GS, M, D, penality=10, SHOW=False, OUT=""):
	R 	= {}
	U 	= {}
	Z 	= {}
	DN 	= {}
	for F, G in GS:
		if "POL" not in M[F]:
			G 	= [s for s in G if s.chrom in D]

			TF 	= M[F]
			U[TF] 	= len(G)
			N 	= float(len(G))
			B 	= [s.get_BIC_mdoel(penality) for s in G ]
			DNN 	= sum([1.0  for s in G if D[s.chrom].searchInterval((s.start, s.stop)) ])
			S 	= float(sum( B ))

			if TF not in R:
				R[TF] 	= (N, S/N)
				DN[TF] 	= DNN / N
			else:
				if R[TF][0] < N:
					R[TF] 	= (N, S/N)
					DN[TF] 	= DNN / N

	for TF in R:
		R[TF]  = [R[TF][1]]
	positions,percents,colors,labels 	= list(),list(),list(),list()
	percents_DNase  	= list()
	delta,height 		= 1.1,1.0
		
	xy 	= [(np.mean(R[r]), r, DN[r])  for r in R]
	xy.sort()
	for i,(mean,TF, mean_DNAse) in enumerate(xy):
		positions+=[i*(delta*2) ]
		percents_DNase+=[mean_DNAse]
		percents+=[mean]


		colors.append(i)
		labels.append(TF)

	


	if SHOW:
		F 	= plt.figure()
		ax 	= F.add_subplot(1,1,1)
		ax.set_title("Cell type: HCT116\nRed: % with DNAse I Signal\nPurple: % with divergent transcription signal")
		ax.barh(positions, percents, height=height/2., color="m" )
		ax.barh([p+(delta/2.) for p in positions], percents_DNase, height=height/2., color="r" )
		
		ax.set_xlim(0,1)
		ax.set_ylim(-delta ,max(positions)+delta*2)
		positions=np.array(positions) +(( height) / 2.)

		ax.set_xticks(np.arange(0,1,0.1))
		ax.set_xticklabels([str(int(x*100)) for x in np.arange(0,1,0.1)])
		ax.set_xlabel("Percent of TF Binding Events with\nDivergent Transcription Prediction")

		ax.set_yticks(positions)
		ax.set_yticklabels(labels)
		ax.grid()
		for i,(pos, per, perDNA) in enumerate(zip(positions, percents, percents_DNase)):
			ax.text(per  , pos -(height/5.)    , str(per*100)[:4]+ "%, " + "N=" + str(U[labels[i]])+ "  ",color="white" ,fontsize=10,
				verticalalignment='center', horizontalalignment="right" )
			ax.text(perDNA  , pos + (delta/2.) -(height/5.)  , str(perDNA*100)[:4]+ "%, " + "N=" + str(U[labels[i]])+ "  ",color="white" ,fontsize=10 ,
				verticalalignment='center', horizontalalignment="right" )


		plt.show()







if __name__ == "__main__":
	Allen_DIR 	= "/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_ChIP/"
	Allen_META 	= "/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_ChIP/metadata.txt"
	HCT116_DNAse= "/Users/joazofeifa/Lab/ChIP/HCT116/DNase/wgEncodeUwDnaseHct116PkRep1.narrowPeak.bed"
	K562 		= "/Users/joazofeifa/Lab/ChIP/K562/K562_SAHADukeDNaseSeq.pk"
	Core_DIR 	= "/Users/joazofeifa/Lab/gro_seq_files/Core2014/EMG_out_ChIP/"
	Core_META 	= "/Users/joazofeifa/Lab/gro_seq_files/Core2014/EMG_out_ChIP/metadata.txt"
	D 			= make_DNAse_searchable(HCT116_DNAse)
	GS 		= load_directory(Allen_DIR, test=False)
	M 		= load_meta_data(Allen_META)
	G 		= graph_percent(GS, M, D, SHOW=True)

