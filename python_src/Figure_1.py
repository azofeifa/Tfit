import load
import matplotlib.pyplot as plt
import numpy as np
import model
import matplotlib as mpl
import matplotlib.cm as cm

def posterior(x, rvs, Type):
	vl 	=  sum([rvs[k].pdf(x,s) for k in (0,1) for s in (-1,1)  ]) /sum([rv.pdf(x,s) for rv in rvs  for s in (-1,1) ])
	
	return vl

def draw(X):
	norm = mpl.colors.Normalize(vmin=0, vmax= 1)
	cmap = plt.get_cmap('PuOr' )
	m = cm.ScalarMappable(norm=norm, cmap=cmap)
		
	means 	= (97, 23)
	sigmas 	= (1.5,1)
	lambdas = (0.3,0.5)
	fps 	= (1.5,0.9)
	pis 	= (0.2,0.4)
	wps 	= (0.2,0.1)
	wlfs 	= (0.05, 0.01)
	wrfs 	= (0.15,0.)
	lfs 	= (X[-1,0],95)
	lrs 	= ( 25,X[0,0])
	rvs 	= [model.component_bidir(means[i], sigmas[i], lambdas[i], wps[i],pis[i] , None,foot_print=fps[i]) for i in range(2) ]
	rvs 	+=[model.component_elongation(lrs[i],means[i], wlfs[i], 0, None, None, None, 0  ) for i in range(2)]
	rvs 	+=[model.component_elongation(means[i], lfs[i],wrfs[i], 1.0, None, None, None, 0  ) for i in range(2)]
	for rv in rvs:
		print rv
	F 	= plt.figure(figsize=(15,10))
	ax 	= F.add_subplot(1,1,1)
	N 	= np.sum(X[:,1:])+700000
	X[:,1]/=N
	X[:,2]/=N
	colorsf 	= [ m.to_rgba(posterior(x,rvs, Type=1 )) for x in X[:,0] ]
	xs 	= np.linspace(X[0,0],X[-1,0],1000)
	ax.bar(X[:,0], X[:,1], color=colorsf, edgecolor=colorsf, alpha=1)
	ax.bar(X[:,0], -X[:,2], color=colorsf, edgecolor=colorsf, alpha=1)
	ax.plot(xs, [sum([rv.pdf(x,1) for rv in rvs]) for x in xs],linewidth=3.,linestyle="--",color="black")
	ax.plot(xs, [sum([-rv.pdf(x,-1) for rv in rvs]) for x in xs],linewidth=3.,linestyle="--",color="black")
	ax_cbar = F.add_axes([0.84,0.2,0.01,0.6])
	cb1 = mpl.colorbar.ColorbarBase(ax_cbar, cmap=cmap,
	                                norm=norm,
	                                orientation='vertical')
	cb1.set_label("\n"+r'$p(k=paused|\hat{\theta})$', fontsize=20)
	ax_cbar.yaxis.tick_right()
	ax_cbar.set_yticklabels(["0", "", "", "","",  "0.5", "","", "", "",  "1"]  )
	
	ax2 = ax_cbar.twinx()
	ax2.set_yticks([])
	ax2.yaxis.set_label_position("left")
	ax.grid()
	ax.set_yticklabels([""]+[str(abs((i ) )) for i in ax.get_yticks()[1:] ] )
	ax.set_ylabel("Density")

	ax.set_xlabel("Relative Genomic Position")

	plt.savefig("/Users/joazofeifa/Lab/Article_drafts/EMG_paper/images/example_gene_fig.svg")
	plt.show()

def write_out(X, OUT=""):
	FHW 	= open(OUT, "w")
	for i in range(X.shape[0]):
		FHW.write(str(X[i,0])+","+str(X[i,1])+","+str(X[i,2])+"\n")
	FHW.close()
def load_IN(FILE):
	L 	= list()
	with open(FILE) as FH:
		for line in FH:
			x,y,r 	= line.strip("\n").split(",")
			L.append((float(x), float(y), float(r)))
	return np.array(L)

if __name__=="__main__":
	WRITE 	= False
	OUT 	= "/Users/joazofeifa/Lab/Article_drafts/EMG_paper/files/Example_Gene.csv"
	if WRITE:
		X 		=  load.grab_specific_region("chr1",8012007, 8033978, 
			pos_file="/Users/joazofeifa//Lab/gro_seq_files/HCT116/bed_graph_files/DMSO2_3.pos.BedGraph", 
			neg_file="/Users/joazofeifa//Lab/gro_seq_files/HCT116/bed_graph_files/DMSO2_3.neg.BedGraph",
			SHOW 	=False, bins=300)
		X[:,0]-=X[0,0]
		X[:,0]/=100.
		write_out(X, OUT=OUT)
	X 		= load_IN(OUT)
	draw(X)