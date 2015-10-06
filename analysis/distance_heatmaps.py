import scipy.special
import scipy.stats 
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import matplotlib.cm as cm

def load(FILE):
	G 	= {}
	TF,fimo="",""
	with open(FILE) as FH:
		for line in FH:
			if "fimo" in line:
				fimo 	= line.strip("\n")

				G[TF][fimo] 	= list()
			elif line[0] == "[":
				line 	= line.strip("\n")
				line 	= line[1:-1]
				G[TF][fimo] 	= [float(x) for x in line.split(",")]
			else:
				TF 	= line.strip("\n")
				if TF not in G:
					G[TF] =	{}

	return G
class motif_distance:
	def __init__(self, fimo, D):
		self.ID 	= fimo
		self.D 		= D
		self.U 		= None
		self.ZERO 	= None
		self.mean 	= None
		self.std 	= None
		self.N 		= None
		self.calc_pvalues()
	def calc_pvalues(self):
		x 			= self.D
		start 		= min(x)
		stop 		= max(x)
		sigma 		= np.std(x)
		mu 			= np.mean(x)
		N 			= len(x)
		y 			= np.random.uniform(start, stop, N)
		z 			= mu/(sigma/math.sqrt(N))
		p 			= 1 - scipy.special.ndtr(z)
		k 			= scipy.stats.ks_2samp(x,y)
		self.U 		= k[1]
		self.ZERO 	= p
		self.mean 	= mu
		self.std 	= sigma
		self.N 		= N
	def bin(self, bins, cmap):
		counts,edges 	= np.histogram(self.D, bins)
		edges 			= (edges[1:] + edges[:-1])/2.
		
		norm    = mpl.colors.Normalize(vmin=min(counts), vmax=max(counts))
		m       = cm.ScalarMappable(norm=norm, cmap=cmap)
		colors  = [m.to_rgba(c) for c in counts] 

		return edges , counts, colors

class TF:
	def __init__(self, name):
		self.name 		= name
		self.distances  = list()
	def insert_disance(self, md):
		self.distances.append(md)
	def get_best(self):
		if self.distances:

			return max([(md.ZERO, md) for md in self.distances])[1]
		return []
	def get_most(self):
		return max([(md.N, md) for md in self.distances])[1]

	
def make_class(G):
	L 	= list()
	for name in G:
		tf 	= TF(name)
		for fimo in G[name]:
			md 	= motif_distance(fimo, G[name][fimo])
			tf.insert_disance(md)
		L.append(tf)
	return L
def display_Hmaps(L):
	L 	= [TF for TF in L if "H3" not in TF.name]
	cmap = cm.Blues

	bins = 50
	N    = 10000
	F    = plt.figure(figsize=(15,10))
	left, width = 0.05, 0.35
	height = 0.1
	bottom = 0.2
	rect_scatter = [left, bottom, width, height]
	axes 	= list()
	min_x,max_x 	= 0,0
	for i in range(6):
		TF 	= L[i]
		md 	= TF.get_most()
		edges, counts, colors 	= md.bin(55, cmap)
		rect_scatter[1]=(i*height + i*0.063 + 0.045)
		ax  	=  plt.axes(rect_scatter)
		if md.U:
			ax.set_title(TF.name + "   " r'$p-value\leq 10^{' + str(int( math.log(md.U,10) )) + "}$" )
		else:
			ax.set_title(TF.name + "   " r'$p-value=0$' )
			
		ax.bar(edges,np.ones((len(edges),)), color=colors, width=(edges[-1]-edges[0])/len(edges) , edgecolor=colors )
		axes.append(ax)
		min_x 	= min(min_x, min(edges) )
		max_x 	= max(max_x, max(edges) )
		ax.set_yticks([])
		
	i 	= 5
	rect_scatter[0] = 0.5
	for j in range(6):
		TF 	= L[j+i]
		md 						= TF.get_most()
		edges, counts, colors 	= md.bin(55, cmap)
		rect_scatter[1]=(j*height + j*0.063 + 0.045)
		ax  	=  plt.axes(rect_scatter)
		if md.U:
			ax.set_title(TF.name + "   " r'$p-value\leq 10^{' + str(int(math.log(md.U,10) )) + "}$" )
		else:
			ax.set_title(TF.name + "   " r'$p-value=0$' )

		ax.bar(edges,np.ones((len(edges),)), color=colors, width=(edges[-1]-edges[0])/len(edges) , edgecolor=colors )

		axes.append(ax)
	min_x,max_x=-750,750
	for ax in axes:
		ax.set_xlim(min_x,max_x)
	


	# mean,std=0,1
	# ax1.set_title("Normal with mean: " + str(mean)+ " and standard deviation: " + str(std) )
	# edges, colors=gen_data_and_colors(mean,std,cmap, bins=bins,N=N)
	# ax1.bar(edges,np.ones((len(edges),)), color=colors, width=(edges[-1]-edges[0])/len(edges) , edgecolor=colors )
	# #ax1.set_xlim(min(edges), max(edges))
	# minX,maxX=min(edges),max(edges)

	# ax2  =  plt.axes(rect_scatter)

	# mean,std=0,4
	# ax2.set_title("Normal with mean: " + str(mean) + " and standard deviation: " + str(std) )
	# edges, colors=gen_data_and_colors(mean,std,cmap, bins=bins,N=N)
	# ax2.bar(edges,np.ones((len(edges),)), color=colors, width=(edges[-1]-edges[0])/len(edges) , edgecolor=colors )
	# #ax2.set_xlim(min(edges), max(edges))
	# minX,maxX=min(min(edges),minX),max(max(edges),maxX)


	# rect_scatter[1]+=height+(height/2.)
	# ax3  =  plt.axes(rect_scatter)

	# mean,std=0,10
	# ax3.set_title("Normal with mean: " + str(mean)+ " and standard deviation: " + str(std) )
	# edges, colors=gen_data_and_colors(mean,std,cmap, bins=bins,N=N)

	# ax3.bar(edges,np.ones((len(edges),)), color=colors, width=(edges[-1]-edges[0])/len(edges) , edgecolor=colors )
	# #ax2.set_xlim(min(edges), max(edges))
	# minX,maxX=min(min(edges),minX),max(max(edges),maxX)


	# rect_scatter[1]+=height+(height/2.)
	# ax4  =  plt.axes(rect_scatter)

	# mean,std=4,2
	# ax4.set_title("Normal with mean: " + str(mean)+ " and standard deviation: " + str(std) )
	# edges, colors=gen_data_and_colors(mean,std,cmap, bins=bins,N=N)

	# ax4.bar(edges,np.ones((len(edges),)), color=colors, width=(edges[-1]-edges[0])/len(edges) , edgecolor=colors )
	# #ax2.set_xlim(min(edges), max(edges))
	# minX,maxX=min(min(edges),minX),max(max(edges),maxX)

	# for ax in (ax1,ax2,ax3,ax4):
	#     ax.set_xlim(minX, maxX)
	axC = F.add_axes([0.9, 0.1, 0.02, 0.8])
	norm = mpl.colors.Normalize(vmin=0, vmax=1)
	cb1 = mpl.colorbar.ColorbarBase(axC, cmap=cmap,
	                                   norm=norm,
	                                   orientation='vertical')
	cb1.set_label('Normalized to Max')
	axC.set_xticklabels([])
	plt.show()

	# pass




if __name__ == "__main__":
	FILE 	= "/Users/joazofeifa/Lab/gro_seq_files/EMG_analysis_out_files/TFMotifToBidirDistance.txt"
	TFS 	= make_class(load(FILE))
	display_Hmaps(TFS)
	pass