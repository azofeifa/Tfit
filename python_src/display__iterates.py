from model import component_bidir
from model import component_elongation
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D
def run(FILE):
	G 	= {}
	D 	= list()
	with open(FILE) as FH:
		for line in FH:
			if "#Data" in line:
				collect=True
			elif collect and "#" in line:
				collect 	= False
				t 			= int(line[1:].split(",")[0])
				G[t] 		= list()
			elif collect:
				D.append(([float(x) for x in line.strip("\n").split("\t") ]))
			elif not collect and "#" in line:
				t 			= int (line[1:].split(",")[0])
				G[t] 		= list()
			else:
				if "U" == line[0]:
					a,b,w,pi 	= line.strip("\n").split(": ")[1].split(",")
					G[t].append(component_elongation(float(a), float(b), float(w), float(pi) , None, "noise", None, None  )  )
				else:
					mu,si,l,w,pi 	= line.strip("\n").split(": ")[1].split(",")
					G[t].append(component_bidir(float(mu), float(si), float(l), float(w), float(pi), None))
	return G, np.array(D)
def draw_line(ax, G, t, xs, al):
	components 	= G[t]
	ysf 		= [sum([c.pdf(x, 1) for c in components])  for x in xs]
	ysr 		= [-sum([c.pdf(x,-1) for c in components])  for x in xs]
	ax.plot(xs, ysf, linewidth=2.5,color="black", alpha=al)
	ax.plot(xs, ysr, linewidth=2.5, color="black", alpha=al)

def get_likelihood_surface(X,Y, D, C):

	def fun(mu, si):
		def LOG(x):
			if x==0:
				return -1000
			return math.log(x)
		C[0].mu 	= mu
		C[1].b 		= mu
		C[2].a 		= mu
		C[0].si 	= si
		C[0].l 		= 2
		f 	= sum([LOG(sum([ B.pdf(D[i,0],1) for B in C ]))*D[i,1] for i in range(D.shape[0])  ])
		r 	= sum([LOG(sum([ B.pdf(D[i,0],-1) for B in C ]))*D[i,2] for i in range(D.shape[0])  ])
		return f+r


	zs = np.array([fun(x,y) for x,y in zip(np.ravel(X), np.ravel(Y))])
	Z = zs.reshape(X.shape)
	return Z
def plot(G,X):
	N 	= np.sum(X[:,1:])
	X[:,1]/=N
	X[:,1]*=0.6
	X[:,2]/=N
	X[:,2]*=0.7
	width 	= (X[0,0]-X[-1,0])/X.shape[0]
	F 		= plt.figure(figsize=(15,10))

	ax 		= F.add_subplot(2,1,1)
	ax.set_title("Covergence of Single Isoform Gene")

	ax.bar(X[:,0], X[:,1], width=width,alpha=0.3, label="Forward Strand (GRO)")
	ax.bar(X[:,0], -X[:,2], width=width, color="red", alpha=0.3, label="Reverse Strand (GRO)")
	ax.set_xlabel("Relative Gene Location TMEM50A (bp)")
	ax.set_ylabel("Read Count Normalized/Probility Density")
	xs 		= np.linspace(X[0,0], X[-1,0], 1000)
	draw_line(ax,G,0, xs,0.3)
	draw_line(ax,G,7, xs,0.3)
	draw_line(ax,G,8, xs,0.3)
	draw_line(ax,G,45, xs,1)
	labels = [item.get_text() for item in ax.get_xticklabels()]
	labels 	= [ str(int(l)*100) for l in ax.get_xticks().tolist()]
	ax.set_xticklabels(labels)
	ann = ax.annotate("T:1 (Start)",
                  xy=(160, 0.005), xycoords='data',
                  xytext=(200, 0.01), textcoords='data',
                  size=12, va="center", ha="center",
                  bbox=dict(boxstyle="round4", fc="w"),
                  arrowprops=dict(arrowstyle="-|>",
                                  connectionstyle="arc3,rad=0",
                                  fc="w"), 
                  )
	ann = ax.annotate("T:7",
                  xy=(115, 0.005), xycoords='data',
                  xytext=(190, 0.0135), textcoords='data',
                  size=12, va="center", ha="center",
                  bbox=dict(boxstyle="round4", fc="w"),
                  arrowprops=dict(arrowstyle="-|>",
                                  connectionstyle="arc3,rad=0",
                                  fc="w"), 
                  )
	ann = ax.annotate("T:8",
                  xy=(100, 0.01), xycoords='data',
                  xytext=(180, 0.0175), textcoords='data',
                  size=12, va="center", ha="center",
                  bbox=dict(boxstyle="round4", fc="w"),
                  arrowprops=dict(arrowstyle="-|>",
                                  connectionstyle="arc3,rad=0",
                                  fc="w"), 
                  )
	ann = ax.annotate("T:25 (Converged)",
                  xy=(88, 0.025), xycoords='data',
                  xytext=(170, 0.022), textcoords='data',
                  size=12, va="center", ha="center",
                  bbox=dict(boxstyle="round4", fc="w"),
                  arrowprops=dict(arrowstyle="-|>",
                                  connectionstyle="arc3,rad=0",
                                  fc="w"), 
                  )

	mus 		= np.linspace(75, 115, 35)
	sigmas 		= np.linspace(0.001,10, 35)
	XX, YY = np.meshgrid(mus, sigmas)
	Z 	= get_likelihood_surface(XX,YY, X,G[47])
	ax2 		= F.add_subplot(2,1,2, projection='3d')
	ax2.set_title("Log Likelihood Landscape")
	
	surf = ax2.plot_surface(XX, YY, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
    	linewidth=0, antialiased=False)
	ax2.set_xlabel("Loading Position, " + r'$\mu$')
	ax2.set_ylabel("Loading Variance, " + r'$\sigma^2$')
	ax2.set_zlabel("log-likelihood, " + r'$l(\Theta)$')


	labels = [item.get_text() for item in ax2.get_xticklabels()]
	labels 	= [ str(int(l)*100) for l in ax2.get_xticks().tolist()]
	ax2.set_xticklabels(labels)
	labels = [item.get_text() for item in ax2.get_yticklabels()]
	labels 	= [ str(int(l)*100) for l in ax2.get_yticks().tolist()]
	ax2.set_yticklabels(labels)
	
	ax2.grid()

	ax.grid()
	plt.tight_layout()
	
	plt.show()
if __name__ == "__main__":
	FILE 	= "/Users/joazofeifa/Lab/EMG/python_src/iterates.tsv"
	G, D 	= run(FILE)
	plot(G, D)