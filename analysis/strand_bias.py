import node, load, time
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats 
import math
def refseq(FILE):
	header 	= True
	G 		= {}
	with open(FILE) as FH:
		for line in FH:
			if not header:
				name, chrom, strand, start, stop 	= line.split("\t")[1:6]
				if chrom not in G:
					G[chrom] 	= list()
				if strand == "+":
					G[chrom].append((int(start)-1000, int (start )+1000 ,  strand ))
				else:
					G[chrom].append((int(stop)-1000, int (stop )+1000 ,  strand ))
			else:
				header=False
	for chrom in G:
		G[chrom] = node.tree(G[chrom])
	return G		


def run(L,R):
	positive 	= list()
	negative 	= list()
	eRNA 		= list()
	for chrom in L:
		if chrom in R:
			for st,sp,S  in L[chrom]:
				for model in S.models:
					if model.wEM > 0.01 :
						FINDS 	=  R[chrom].searchInterval((model.start, model.stop))
						if FINDS is None:
							eRNA.append(model)
						else:
							strand 	= FINDS[0][2]
							if strand == "+":
								positive.append(model)
							else:
								negative.append(model)

	axes,bps 	= list(), list()
	F 	= plt.figure(figsize=(10,6))
	ax 	= F.add_subplot(1,2,1)
	ax.set_title("Strand Bias Distributions")
	pos, neg, ERNA 	= [m.pi for m in positive],[m.pi for m in negative],[m.pi for m in eRNA]

	pos_neg 		= scipy.stats.ks_2samp(pos,neg)
	pos_ERNA 		= scipy.stats.ks_2samp(pos,ERNA)
	neg_ERNA 		= scipy.stats.ks_2samp(ERNA,neg)
	X 				= [1,2,3,1,3]
	Y 				= [1,1.01,1.1,1.19,1.19]
	
	props = {'arrowstyle':'|-|',\
	         'shrinkA':1,'shrinkB':1,'lw':2,"mutation_scale":5.}

	ax.annotate('', xy=(1,1.1), xytext=(2,1.1), arrowprops=props)
	ax.annotate('', xy=(2,1.2), xytext=(3,1.2), arrowprops=props)
	ax.annotate('', xy=(1,1.3), xytext=(3,1.3), arrowprops=props)
	ax.annotate(r'$p<10^{' + str(int(math.log(pos_neg[1],10))) + "}$", 
		xy=(1.35,1+0.102), zorder=300,fontsize=12)
	ax.annotate(r'$p<10^{' + str(int(math.log(pos_ERNA[1],10))) + "}$", 
		xy=(2.35,1+0.202), zorder=300,fontsize=12)
	ax.annotate(r'$p<10^{' + str(int(math.log(pos_ERNA[1],10))) + "}$", 
		xy=(1.8,1+0.302), zorder=300,fontsize=12)
	
	
	ax.set_ylim(-0.1,1.5)

	bp = ax.boxplot((pos, neg, ERNA),patch_artist=True)
	
	axes.append(ax)

	bps.append(bp)
	ax 	= F.add_subplot(1,2,2)
	ax.set_title("Foot Print Parameter Distributions")
	pos, neg, ERNA 	= [m.fp for m in positive],[m.fp for m in negative],[m.fp for m in eRNA]
	pos_neg 		= scipy.stats.ks_2samp(pos,neg)
	pos_ERNA 		= scipy.stats.ks_2samp(pos,ERNA)
	neg_ERNA 		= scipy.stats.ks_2samp(ERNA,neg)
	
	bp = ax.boxplot((pos, neg, ERNA),patch_artist=True)

	ax.annotate('', xy=(1,550), xytext=(2,550), arrowprops=props)
	ax.annotate('', xy=(2,600), xytext=(3,600), arrowprops=props)
	ax.annotate('', xy=(1,650), xytext=(3,650), arrowprops=props)
	ax.annotate(r'$p=' + str(0.1) + "$", 
		xy=(1.35,550+2), zorder=300,fontsize=12)
	ax.annotate(r'$p<10^{' + str(int(math.log(pos_ERNA[1],10))) + "}$", 
		xy=(2.35,600+2), zorder=300,fontsize=12)
	ax.annotate(r'$p<10^{' + str(int(math.log(pos_ERNA[1],10))) + "}$", 
		xy=(1.8,650+2), zorder=300,fontsize=12)
	ax.set_ylim(-100,750)

	axes.append(ax)
	bps.append(bp)
	## change outline color, fill color and linewidth of the boxes
	for i,bp in enumerate(bps):
		ax 	= axes[i]
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
		ax.set_xticklabels(["Forward Strand\nPromoter", "Reverse Strand\nPromoter", "eRNA\nUnannotated"])
		ax.grid()
	plt.tight_layout()
	plt.savefig("/Users/joazofeifa/Lab/Talks/2015/CSHL/Parameter_Comparison")
	plt.show()




if __name__ == "__main__":
	# f 	= lambda x: (1. / math.sqrt(2*math.pi*10)  )*math.exp(-pow(x-50,2)/ (2*10) )
	# g 	= lambda x: (1. / math.sqrt(2*math.pi*10)  )*math.exp(-pow(x-60,2)/ (2*10) )
	
	# print f(50)
	# print g(50)
	DIR 			="/Users/joazofeifa/Lab/gro_seq_files/HCT116/EMG_out_files/"
		
	REF 						= "/Users/joazofeifa/Lab/genome_files/RefSeqHG19.txt"
	DMSO2_3 						="Allen2014_DMSO2_3-1_bidirectional_hits_intervals.bed"
	DMSO2_3_L,DMSO2_3_G 			= load.load_model_fits_bed_file(DIR+DMSO2_3)
		
	R 					= refseq(REF)
	run(DMSO2_3_L,R)