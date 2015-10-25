import copy as cp
import time
import matplotlib.pyplot as plt
import numpy as np
import math
class node:
	def __init__(self, start, stop, chrom,gene, method, ct):
		self.start,self.stop, self.chrom 			= start,stop, chrom
		self.method,self.cell_tissue, self.gene 	= method, ct, gene
		self.link 									= None
class concat_node:
	def __init__(self, start, stop):
		self.start, self.stop, self.nodes 			= start, stop, list()
	def union(self, a):
		self.start 	= min(self.start, a.start)
		self.stop 	= max(self.stop, a.stop)
		self.nodes.append(a)
	def intersect(self, a):
		self.start 	= max(self.start, a.start)
		self.stop 	= min(self.stop, a.stop)
		self.nodes.append(a)
	def get_overlaps(self):
		copy_nodes,C 	= cp.copy(self.nodes), list()
		while copy_nodes:
			cn 		= copy_nodes.pop()
			current = concat_node(cn.start, cn.stop)
			C.append(current)
			j, N 	= 0,len(copy_nodes)
			while j < N:
				if current.start < copy_nodes[j].stop and current.stop > copy_nodes[j].start:
					current.union(copy_nodes[j])
					copy_nodes 	= copy_nodes[:j]+copy_nodes[j+1:]
				N 	= len(copy_nodes)
				j+=1
		return C
	def get_link_count(self, ALL=False, intersect=False, union=False ):
		if ALL:
			return len(self.nodes)
		if union:
			links 	= [(l.link.start, l.link) for l in self.nodes ]
			links.sort()
			links 	= [l for st, l in links]
			j,N 	= 0,len(links)
			ct 	= 0
			while j < N:
				start, stop = links[j].start,links[j].stop
				while j < N and stop > links[j].start and start < links[j].stop :
					start,stop 	= min(start, links[j].start ), max(stop, links[j].stop)
					j+=1
				ct+=1
				j+=1
			return ct
		if intersect:
			pass


		return 0



def load(FILE, test=False):
	header,t,G =True,0,{}
	with open(FILE) as FH:
		for line in FH:
			if test and t > 100:
				break
			t+=1
			if not header:
				(a_chrom,a_st, a_sp, b_chrom, b_st, b_sp, Agene,Bgene, 
					organism, Cell_Tissue, method, cs1, cs2, cf, ID) 	= line.strip("\n").split("\t")
				a,b 	= (node(int(a_st), int(a_sp), a_chrom, Agene, method, Cell_Tissue), 
					node(int(b_st), int(b_sp), b_chrom, Bgene, method, Cell_Tissue))
				a.link,b.link 	= b,a
				for u in (a,b):
					if u.chrom not in G:
						G[u.chrom] 	= list()
					G[u.chrom].append((u.start, u))
			else:
				header=False
	for chrom in G:
		G[chrom].sort()
		G[chrom] 	= [u for st, u in G[chrom]]
	return G
def union(G):
	U 		= {}
	for chrom in G:
		A 		= G[chrom]
		U[chrom]= list()
		j,N = 0,len(A)
		while j < N:
			current 	= concat_node(A[j].start, A[j].stop)
			while j < N and current.start < A[j].stop and current.stop > A[j].start:
				current.union(A[j])
				j+=1
			U[chrom].append(current)
			j+=1
	return U
def intersect(U):
	I 	= {}
	for chrom in U:
		I[chrom]=list()
		for u in U[chrom]:
			C 	= u.get_overlaps()
			I[chrom]+=C
	return I
def bin(ax, data,bins, title=""):
	counts,edges 	= np.histogram(data, bins=bins)
	counts 			= [math.log(ct+1, 10) for ct in counts]
	edges 			= (edges[1:] + edges[:-1])/2.
	ax.bar(edges, counts, width = (edges[-1] - edges[0])/bins )
	ax.set_title(title)
	ax.grid()
	
	
def get_degree_distribution(U, I):
	bins 			= 50
	F 				= plt.figure(figsize=(15,10))
	ax1 			= F.add_subplot(2,2,1)
	bin(ax1, [u.get_link_count(ALL=True) for chrom in I for u in I[chrom]], bins, title="intersect to all" )
	ax2 			= F.add_subplot(2,2,2)
	bin(ax2, [u.get_link_count(union=True) for chrom in I for u in I[chrom]], bins, title="intersect to union" )
	ax3 			= F.add_subplot(2,2,3)
	bin(ax3, [u.get_link_count(ALL=True) for chrom in U for u in U[chrom]], bins, title="union to all" )
	ax4 			= F.add_subplot(2,2,4)
	bin(ax4, [u.get_link_count(union=True) for chrom in U for u in U[chrom]], bins, title="union to union" )
	
	plt.show()

if __name__ == "__main__":
	FILE 	= "/Users/joazofeifa/Lab/3D_interactions/4D_genome_downloads/4DGenome_ChIA-PET.txt"
	G 		= load(FILE, test=False)
	U 		= union(G)
	I 		= intersect(U)
	get_degree_distribution(U,I)