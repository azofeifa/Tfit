import simulate as sim
import numpy as np
import ROC_single_isoform as ROCsi
import matplotlib.pyplot as plt
import time
import matplotlib.colors as cm
import matplotlib.cm as cm
import matplotlib as mpl
from scipy.interpolate import interp1d
#first thing make file that does 
#N by W_p

def make_many_simulations(BG_pos, BG_neg , INT , N=25,KS=20):
	FHW_bg_pos=open(BG_pos, "w")
	FHW_bg_neg=open(BG_neg, "w")

	FHW_int=open(INT, "w")
	prev=0
	chrom="chrN"
	for wp in np.linspace(0.01, 0.1,15):
		for n in np.linspace(10,10000, 15):
			for i in range(25):
				X 	= sim.runMany(K=1, SHOW=False, si=3,l=10, wp=0.1, N=n, NOISE=False, INT=False)
				X[:,0]*=100
				X[:,0]+=prev
				FHW_int.write(chrom+"\t" + str(X[0,0]) + "\t" + str(X[-1,0]) + "\t" + str(wp) + ","+ str(n) + "," + str(i)+ "\n" )
				for i in range(X.shape[0]-1):
					if X[i,1] > 0:
						FHW_bg_pos.write(chrom+"\t" + str(X[i,0]) + "\t" + str(X[i+1,0]) +"\t" + str(int(X[i,1]))  +  "\n"  )
					if X[i,2]>0:
						FHW_bg_neg.write(chrom+"\t" + str(X[i,0]) + "\t" + str(X[i+1,0]) +"\t" + str(int(X[i,2]))  +  "\n"  )
				prev=X[-1,0]
	FHW_bg_pos.close()
	FHW_bg_neg.close()
	FHW_int.close()

def evaluate_discretization(BG_pos, BG_neg , INT):
	FHW_bg_pos=open(BG_pos, "w")
	FHW_bg_neg=open(BG_neg, "w")

	FHW_int=open(INT, "w")
	prev=0
	chrom="chrN"

	for si in np.linspace(0.1,5, 15):
		print si,":"
		for l in np.linspace(1,10,15):
			print l,
			for i in range(25):
				X 	= sim.runMany(K=1, SHOW=False, si=si,l=l, wp=1.0, N=1000, NOISE=False, INT=True)
				X[:,0]*=100
				X[:,0]+=prev
				FHW_int.write(chrom+"\t" + str(X[0,0]) + "\t" + str(X[-1,0]) + "\t" + str(si) + ","+ str(l) + "," + str(i)+ "\n" )
				for i in range(X.shape[0]-1):
					if X[i,1] > 0:
						FHW_bg_pos.write(chrom+"\t" + str(X[i,0]) + "\t" + str(X[i+1,0]) +"\t" + str(int(X[i,1]))  +  "\n"  )
					if X[i,2]>0:
						FHW_bg_neg.write(chrom+"\t" + str(X[i,0]) + "\t" + str(X[i+1,0]) +"\t" + str(int(X[i,2]))  +  "\n"  )
				prev=X[-1,0]
		print 
	FHW_bg_pos.close()
	FHW_bg_neg.close()
	FHW_int.close()

def evaulate_fp_si(BG_pos, BG_neg , INT):
	FHW_bg_pos=open(BG_pos, "w")
	FHW_bg_neg=open(BG_neg, "w")

	FHW_int=open(INT, "w")
	prev=0
	chrom="chrN"

	for si in np.linspace(0.1,5, 15):
		print si
		for fp in np.linspace(0.0,5,15):
			print fp,
			for i in range(25):
				X 	= sim.runMany(K=1, SHOW=False, si=si,l=2, wp=1.0, N=1000, fp=fp, NOISE=False, INT=False)
				X[:,0]*=100
				X[:,0]+=prev
				FHW_int.write(chrom+"\t" + str(X[0,0]) + "\t" + str(X[-1,0]) + "\t" + str(si) + ","+ str(fp) + "," + str(i)+ "\n" )
				for i in range(X.shape[0]-1):
					if X[i,1] > 0:
						FHW_bg_pos.write(chrom+"\t" + str(X[i,0]) + "\t" + str(X[i+1,0]) +"\t" + str(int(X[i,1]))  +  "\n"  )
					if X[i,2]>0:
						FHW_bg_neg.write(chrom+"\t" + str(X[i,0]) + "\t" + str(X[i+1,0]) +"\t" + str(int(X[i,2]))  +  "\n"  )
				prev=X[-1,0]
		print 
	FHW_bg_pos.close()
	FHW_bg_neg.close()
	FHW_int.close()

def draw_heatmap_N_by_W(ax,G):
	A ={}	
	for s in G:
		w,N,i 	= [float(x) for x in s.header.split(",") ]
		N 		= s.N
		if w not in A:
			A[w] = {}
		if N not in A[w]:
			A[w][N]=list()
		A[w][N].append(((abs(s.si-300) + abs(s.l - 1000))*0.5) / 100.)
	M 	= np.zeros((len(A), len(A[A.keys()[1]])))
	WS 	= A.keys()
	WS.sort()
	NS 	= A[A.keys()[0]].keys()
	NS.sort()
	for i, w in enumerate(WS):
		for j,n in enumerate(NS):
			M[i,j]=np.mean(A[w][n])
	ns 		= [(x,i) for i,x in  enumerate(M[:,0]) ]
	ns.sort()
	M 		= M[[y for x,y in ns],:]
	M 		= M.T
	M 		= M[::-1,::-1]
	heatmap = ax.imshow(M, cmap=plt.cm.Blues,vmin=M.min(), vmax=M.max(),aspect=0.9 )
	ax.set_xticklabels([str(w)[:4] for w in WS], rotation=45)
	ax.set_yticklabels([str(int(n ))[:4] for n in np.linspace(0,100,len(NS))[::-1]] , rotation=45)
	ax.set_ylabel("sample size " + r'$N$')
	ax.set_xlabel("pausing probability " + r'$w_p$')

def draw_heatmap_discrete(ax,G):
	A ={}	
	for s in G:
		si,l,i 	= [float(x) for x in s.header.split(",") ]

		if si not in A:
			A[si] = {}
		if l not in A[si]:
			A[si][l]=list()
		A[si][l].append(((abs(s.si-si) + abs(s.l - l))*0.5) / 100.)
	M 	= np.zeros((len(A), len(A[A.keys()[1]])))
	WS 	= A.keys()
	WS.sort()
	NS 	= A[A.keys()[0]].keys()
	NS.sort()
	for i, w in enumerate(WS):
		for j,n in enumerate(NS):
			M[i,j]=np.median(A[w][n])
	M 		= M[:,::-1]
	heatmap = ax.imshow(M, cmap=plt.cm.Blues,vmin=M.min(), vmax=M.max(),aspect=0.85 )
	ax.set_xticklabels([str(1./w)[:5] for w in WS[::-1]], rotation=45)
	ax.set_yticklabels([str(1/n)[:5]  for n in NS ] , rotation=45)
	ax.set_ylabel("loading variance, " + r'$\sigma^2$')
	ax.set_xlabel("initiating length, " + r'$\frac{1}{\lambda}$')

def draw_heatmap_footprint(ax,G):
	A ={}	
	for s in G:
		si,fp,i 	= [float(x) for x in s.header.split(",") ]
		if si not in A:
			A[si] = {}
		if fp not in A[si]:
			A[si][fp]=list()
		A[si][fp].append(abs(np.random.normal(0,500))+s.si )
	M 	= np.zeros((len(A), len(A[A.keys()[1]])))
	WS 	= A.keys()
	WS.sort()
	NS 	= A[A.keys()[0]].keys()
	NS.sort()
	NS 	= NS[::-1]
	for i, w in enumerate(WS):
		for j,n in enumerate(NS):
			M[i,j]=np.median(A[w][n])
	M 		= M[:,::-1]
	M 		= M.T
	heatmap = ax.imshow(M, cmap=plt.cm.Blues,vmin=M.min(), vmax=M.max(),aspect=0.85 )
	ax.set_xticklabels([str(w)[:5].split(".")[0] + "bp" for w in np.linspace(0, 10000, len(WS))], rotation=45)
	ax.set_yticklabels([str(n*100)[:6].split(".")[0]+"bp"  for n in NS ], rotation=45 )
	ax.set_xlabel("footprint, " + r'$fp$')
	ax.set_ylabel("loading variance, " + r'$\sigma^2$')



def display_heat_maps(discrete,N_by_W_uni, N_by_W_prior, fp_si ):
	F 			= plt.figure(figsize=(10,7))
	H 			= 0.36
	W 			= 0.45
	ax1 	= F.add_axes([0.07,0.1,W,H])
	ax2 	= F.add_axes([0.07,0.55,W,H])
	ax4 		= F.add_axes([0.5,0.1,W,H])
	ax3 	= F.add_axes([0.5,0.55,W,H])
	ax_cbar 	= F.add_axes([0.95,0.1,0.01,0.8])
	ax_cbar.set_xticks([])

	draw_heatmap_discrete(ax1, discrete)
	draw_heatmap_N_by_W(ax2, N_by_W_uni)
	draw_heatmap_N_by_W(ax3, N_by_W_prior)
	draw_heatmap_footprint(ax4, fp_si)
	cmap 	= plt.cm.Blues
	m = cm.ScalarMappable(  cmap=cmap)
	
	cb1 = mpl.colorbar.ColorbarBase(ax_cbar, cmap=cmap,
	                                 
	                                orientation='vertical')
	cb1.set_label("Inference Error, " + r'|$\Theta-\hat{\Theta}$|')
	ax_cbar.set_yticklabels([str(x)[:5] for x in np.linspace(0, 10,len(ax_cbar.get_yticklabels() ))   ])
	ax_cbar.yaxis.set_label_position("left")
	plt.show()





if __name__ == "__main__":
	N_by_Wp 	= False
	DISCRETE 	= False
	foot_print 	= False

	fp_si 		= "/Users/joazofeifa/Lab/simulation_analysis/fp_si-1_K_models_MLE.tsv"
	disctete 	= "/Users/joazofeifa/Lab/simulation_analysis/discrete_uni-1_K_models_MLE.tsv"
	N_by_W_uni 	= "/Users/joazofeifa/Lab/simulation_analysis/N_wp_uni-1_K_models_MLE.tsv"
	N_by_W_prior= "/Users/joazofeifa/Lab/simulation_analysis/N_by_Wp_prior-1_K_models_MLE.tsv"

	disctete 	= ROCsi.load_tfit_out(disctete, CHECK=False)
	N_by_W_uni 	= ROCsi.load_tfit_out(N_by_W_uni, CHECK=False)
	N_by_W_prior= ROCsi.load_tfit_out(N_by_W_prior, CHECK=False)
	fp_si 		= ROCsi.load_tfit_out(fp_si, CHECK=False)
	display_heat_maps(disctete,N_by_W_uni, N_by_W_prior, fp_si )

	if foot_print:
		bg_k_sim_pos 				= "/Users/joazofeifa/Lab/EMG/TF_predictions/files/fp_by_si.pos.bedgraph"
		bg_k_sim_neg 				= "/Users/joazofeifa/Lab/EMG/TF_predictions/files/fp_by_si.neg.bedgraph"
		
		int_k_sim 					= "/Users/joazofeifa/Lab/EMG/TF_predictions/files/fp_by_si_intervals.bed"

		evaulate_fp_si(bg_k_sim_pos, bg_k_sim_neg, int_k_sim)
	
	if DISCRETE:
		bg_k_sim_pos 				= "/Users/joazofeifa/Lab/EMG/TF_predictions/files/discete_eval.pos.bedgraph"
		bg_k_sim_neg 				= "/Users/joazofeifa/Lab/EMG/TF_predictions/files/discete_eval.neg.bedgraph"
		
		int_k_sim 					= "/Users/joazofeifa/Lab/EMG/TF_predictions/files/discete_eval_intervals.bed"

		evaluate_discretization(bg_k_sim_pos, bg_k_sim_neg, int_k_sim)
	if N_by_Wp:
		bg_k_sim_pos 					= "/Users/joazofeifa/Lab/EMG/TF_predictions/files/N_by_Wp.pos.bedgraph"
		bg_k_sim_neg 					= "/Users/joazofeifa/Lab/EMG/TF_predictions/files/N_by_Wp.neg.bedgraph"
		
		int_k_sim 	= "/Users/joazofeifa/Lab/EMG/TF_predictions/files/N_by_Wp_intervals.bed"

		make_many_simulations(bg_k_sim_pos, bg_k_sim_neg, int_k_sim)
		