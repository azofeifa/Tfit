import simulate as sim
import model
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import time
def run(mu=0, si=1, l=0.1, pi=0.5 , we=0.5,wl=0.25, wr=0.25,  N=2,NS=(100,1000), OUT=""):
	FHW = open(OUT, "w")
	FHW.write("~" + str(mu) + "," + str(si) + "," + str(l) + "," + str(we) + "," + str(pi) + "," + str(wl) + "," + str(wr) + "\n" )
	for n in NS:
		FHW.write("#" + str(int(n)) + "\n")
		for b in range(N):		
			print n,b
			FHW.write("*" + str(b) + "\n")
			X 	= sim.runOne(mu=mu, s=si, l=1.0/l, lr=100, ll=-100, we=we,wl=wl, wr=wr, pie=pi, pil=0.1, 
				pir=0.9, N=int(n), SHOW=False , bins=200, noise=False, foot_print = 0 )
			clf = model.EMGU(noise=True, K=1,noise_max=0.,moveUniformSupport=0,max_it=200, cores=1, 
				seed=False,foot_print=0)
			t 	= clf.fit(X)
			for rv in clf.rvs:
				FHW.write("conv: " + str(t) + "\n")
				FHW.write(rv.__str__() + "\n")
				FHW.flush()
def load_sim(FILE):
	G 	={}
	with open(FILE) as FH:
		mu 	= None
		for line in FH:	
			if "~" == line[0]:
				gmu,gsi, gl, gw, gpi, gwe, gwl 	= [float(x) for x in line[1:].strip("\n").split(",")]
			elif "#" == line[0]:
				N 	= float(line[1:].strip("\n"))
				G[N]=list()
			elif "*"==line[0]:
				if mu is not None:
					G[N].append((mu,si,l, w2, pi2))
			elif "N:"==line[:2]:
				mu,si,l, w2, pi2 	= [float(x) for x in line.split(": ")[1].strip("\n").split(",") ]

			elif "U:"==line[:2]:
				a,b,w2,pi2 		= [float(x) for x in line.split(": ")[1].strip("\n").split(",") ]
	return G,gmu,gsi, gl, gw, gpi



def draw_errors(ax, no_prior, prior, gmu,gsi, gl, gw, gpi ):
	

	colors 	= ("red", "blue")
	labels 	= ("No Prior, MLE", "Week Prior, MAP")
	for i,NP in enumerate((no_prior, prior)):
		del NP[50.0]
		del NP[70.0]
		del NP[111.0]
		NS 	 	= NP.keys()	
		NO_MUS 	= [np.mean([ (abs(mu-gmu) + abs(si - gsi) + abs(l-gl) + abs(gpi - pi) + abs(gw-w) )  for mu,si, l, w, pi in NP[N] ]) - 1.1 for N in NP] 
		NO_STD 	= [np.std([ (abs(mu-gmu) + abs(si - gsi) + abs(l-gl) + abs(gpi - pi) + abs(gw-w) )  for mu,si, l, w, pi in NP[N] ]) for N in NP] 

		xy 		= [(n,mu, std) for n, mu,std in zip(NS, NO_MUS, NO_STD)]
		xy.sort()
		NS, NO_MUS,NO_STD 	= [x for x,y,z in xy],[y for x,y,z in xy],[z for x,y,z in xy]
				
		ax.plot(NS, NO_MUS, linewidth=3., label=labels[i])
		ax.fill_between(NS, [max(y-z,0) for x,y,z in xy],[y+z for x,y,z in xy], alpha=0.2, color=colors[i]   )

		ax.set_xlim(0,400)
	ax.set_ylabel("Error " + r'$|\hat{\theta} - \theta|$')
	ax.set_xlabel("Number of Data Points")
	ax.set_ylim(-2,10)
	ax.set_yticks(np.arange(-2,10,2))
	ax.set_yticklabels([""]+[str(x) for x in np.arange(-2,10,2)[1: ]] )
	
	ax.grid()
	ax.legend()

def draw_sim_one(ax, mu,si, l, w, pi):
	X 	= sim.runOne(mu=mu, s=2, l=1.0/l, lr=100, ll=-100, we=0.5,wl=0.25, wr=0.25, pie=0.5, pil=0.1, 
				pir=0.9, N=300, SHOW=False , bins=200, noise=False, foot_print = 1 )
	noise 	= np.random.uniform(-135,145, 50)
	noise2 	= np.random.uniform(-135,145, 50)
	
	counts,edges 	= np.histogram(noise, bins=200, normed=1)
	counts2,edges2 	= np.histogram(noise2, bins=200, normed=1)
	
	counts*=0.05
	counts2*=0.05
	
	rvs 			= [model.component_bidir( 0, 2, l, 0.5,0.5 , None,foot_print=1 )]
	rvs 			+=[model.component_elongation( 0, 100,  0.25,1.0 ,None , None, None, None)]
	rvs 			+=[model.component_elongation( -100, 0,  0.25,0.0 , None , None, None, None  )]


	X[:,1:]/=np.sum(X[:,1:])
	w 	= (X[-1,0] - X[0,0]) / X.shape[0]
	ax.bar(X[:,0], X[:,1], width=w, color="blue", edgecolor="blue",alpha=0.3)		
	ax.bar(X[:,0], -X[:,2], width=w, color="red", edgecolor="red",alpha=0.3)		
	ax.bar(edges[1:], counts,width=w, color="blue", edgecolor="blue",alpha=0.3)
	ax.bar(edges2[1:], -counts2, width=w, color="red", edgecolor="red",alpha=0.3)
	
	xs 	= np.linspace(X[0,0], X[-1,0], 1000)
	ysf = [sum([rv.pdf(x,1) for rv in rvs ]) for x in xs]
	ysr = [sum([-rv.pdf(x,-1) for rv in  rvs ]) for x in xs]
	ax.plot(xs,ysf,linewidth=2.5,linestyle="-", color="black")
	ax.plot(xs,ysr,linewidth=2.5,linestyle="-", color="black", label="Model")
	ax.set_yticks(np.linspace(min(ysr)-0.01, max(ysf)+0.01, 12))
	ax.set_yticklabels([""]+ [str(abs(x))[:5] for x in np.linspace(min(ysr)+0.01, max(ysf)-0.01, 10)] + [""] )
	ax.set_ylim(min(ysr)-0.01, max(ysf)+0.01)
	ax.grid()
	ax.set_ylabel("Density")
	ax.set_xticklabels([str(int(10*x)) for x in ax.get_xticks() ])
	ax.set_xlabel("Relative Genomic Coordinate")
	ax.set_title("Fitted Model " + r'$|\theta-\hat\theta|=0.05$')
	ax.legend()
	pass
def posterior(x, rvs, Type):
	vl 	=  sum([rvs[k].pdf(x,s) for k in (Type,) for s in (-1,1)  ]) /sum([rv.pdf(x,s) for rv in rvs[Type+1:Type+3] for s in (-1,1) ])
	
	return vl



def draw_posterior(ax):
	X1 	= sim.runOne(mu=-40, s=1, l=10, lr=40, ll=-200, we=0.5,wl=0.25, wr=0.25, pie=0.5, pil=0.1, 
				pir=0.9, N=3000, SHOW=False , bins=200, noise=False, foot_print = 1 )
	X2 	= sim.runOne(mu=40, s=1, l=10, lr=200, ll=-40, we=0.5,wl=0.25, wr=0.25, pie=0.5, pil=0.1, 
				pir=0.9, N=3000, SHOW=False , bins=200, noise=False, foot_print = 1 )

	rvs 			= [model.component_bidir( -40, 1, 0.1, 0.5,0.5 , None,foot_print=1 )]
	rvs 			+=[model.component_elongation( -40, 200,  0.25,1.0 ,None , None, None, None)]
	rvs 			+=[model.component_elongation( -200, -40,  0.25,0 , None , None, None, None  )]
	
	rvs 			+= [model.component_bidir( 40, 1, 0.1, 0.5,0.5 , None,foot_print=1 )]

	rvs 			+=[model.component_elongation( 40, 200,  0.25,1.0 ,None , None, None, None)]
	rvs 			+=[model.component_elongation( -200, 40,  0.25,0 , None , None, None, None  )]

	norm = mpl.colors.Normalize(vmin=0, vmax= 1)
	cmap = plt.get_cmap('PuOr' )
	m = cm.ScalarMappable(norm=norm, cmap=cmap)
	
	KS 	= (0,3)
	MIN, MAX 	= 0,0
	for i,X in enumerate((X1, X2)) :
		w 	= (X[-1,0] - X[0,0]) / X.shape[0]
		xs 	= np.linspace(X[0,0],X[-1,0], X.shape[0] )
		MAX = max(MAX, max(X[:,1]))
		MIN = min(MIN, min(-X[:,2]))
		colorsf 	= [ m.to_rgba(posterior(x,rvs, Type=KS[i] )) for x in xs ]
		
		ax.bar(X[:,0], X[:,1], width=w, color=colorsf, edgecolor=colorsf,alpha=0.3)		
		ax.bar(X[:,0], -X[:,2], width=w, color=colorsf, edgecolor=colorsf,alpha=0.3)		
	ax.set_ylim(MIN, MAX)
	ax.set_yticklabels([""]+[str(abs(int(i ) )) for i in ax.get_yticks()[1:] ] )
	ax.set_ylabel("Read Coverage")

	ax.set_title(r'$p(k|\hat{\Theta})$')
	ax.set_xlabel("Relative Genomic Position")
	ax.grid()
	
def draw_fig(no_prior, prior,gmu,gsi, gl, gw, gpi ):
	F 			= plt.figure(figsize=(10,7))
	ax_error 	= F.add_subplot(2,2,3)
	ax_example 	= F.add_subplot(2,2,1)
	ax_post 	= F.add_subplot(2,2,2)
	draw_sim_one(ax_example, gmu,gsi, gl, gw, gpi)	
	draw_errors(ax_error, no_prior, prior, gmu,gsi, gl, gw, gpi )
	draw_posterior(ax_post)
	plt.show()


def make_many_simulations(BG_pos="", BG_neg="", INT="", N=10,KS=3):
	FHW_bg_pos=open(BG_pos, "w")
	FHW_bg_neg=open(BG_neg, "w")

	FHW_int=open(INT, "w")
	prev=0
	chrom="chrN"
	for K in  range(1, KS):
		print K
		for n in range(N):
			X 	= sim.runMany(K=K, SHOW=False)
			X[:,0]*=100
			X[:,0]+=prev
			FHW_int.write(chrom+"\t" + str(X[0,0]) + "\t" + str(X[-1,0]) + "\t" + str(K) + ","+ str(n) + "\n" )
			for i in range(X.shape[0]-1):
				if X[i,1] > 0:
					FHW_bg_pos.write(chrom+"\t" + str(X[i,0]) + "\t" + str(X[i+1,0]) +"\t" + str(int(X[i,1]))  +  "\n"  )
				if X[i,2]>0:
					FHW_bg_neg.write(chrom+"\t" + str(X[i,0]) + "\t" + str(X[i+1,0]) +"\t" + str(int(X[i,2]))  +  "\n"  )


			prev=X[-1,0]
	FHW_bg_pos.close()
	FHW_bg_neg.close()
	FHW_int.close()


if __name__ == "__main__":
	OUT 							= "/Users/joazofeifa/Lab/EMG/simulation_figures/files/simulation_sample_size_prior.csv"
	no_prior,gmu,gsi, gl, gw, gpi 	= load_sim("/Users/joazofeifa/Lab/EMG/simulation_figures/files/simulation_sample_size_no_prior.csv")
	prior,gmu,gsi, gl, gw, gpi 		= load_sim("/Users/joazofeifa/Lab/EMG/simulation_figures/files/simulation_sample_size_prior.csv")
	bg_k_sim_pos 					= "/Users/joazofeifa/Lab/EMG/TF_predictions/files/many_k.pos.bedgraph"
	bg_k_sim_neg 					= "/Users/joazofeifa/Lab/EMG/TF_predictions/files/many_k.neg.bedgraph"
	
	int_k_sim 	= "/Users/joazofeifa/Lab/EMG/TF_predictions/files/many_k_intervals.bed"
	make_many_simulations(BG_pos=bg_k_sim_pos, BG_neg=bg_k_sim_neg, INT=int_k_sim)
	#draw_fig(no_prior, prior,gmu,gsi, gl, gw, gpi)


	#run(mu=0,si=1,l=0.1,N=25, NS=np.linspace(10,1000,50), OUT=OUT)