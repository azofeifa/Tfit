import os, node
import numpy as np 
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt

def unicode(c):
	try:
		c.decode('utf-8')
		return True
	except UnicodeError:
		return False

def load(DIR):
	G 	= {}
	t = 0
	for f in os.listdir(DIR):
		if "prelim" in f:
			with open(DIR+f) as FH:
				for line in FH:
					if "#" != line[0]:
						try:
							chrom,start, stop 	= line.strip("\n").split("\t")
							if chrom not in G and unicode(chrom) and len(chrom) < 6:
								G[chrom] 	= list()
							elif unicode(chrom) and len(chrom) < 6:
								if (int(start) >0 ):
									G[chrom].append((int(start), int(stop)))

						except:
							pass
			t+=1
			if t >3:
				break
	M 	= {}
	for chrom in G:
		G[chrom].sort()
		i 	= 0
		N 	= len(G[chrom])
		if N:
			g 	= G[chrom]
			M[chrom] 	= list()
			while i < N:
				o_st,o_sp 	= g[i]
				while i < N and g[i][0] < o_sp and g[i][1] > o_st:
					o_st,o_sp 	= min(o_st, g[i][0]), max(o_sp, g[i][1])
					i+=1
				M[chrom].append((o_st,o_sp))
	D = {}
	for chrom in M:
		if len(M[chrom]):
			D[chrom] 	= node.tree(M[chrom])
	print DIR
	return M, D
def merge_ALL(EXPS_M_D):
	G 		= {}
	LST 	= [(chrom ,st, sp) for a_m,d,e in EXPS_M_D for chrom in a_m for st,sp in a_m[chrom]]
	for chrom,st, sp in LST:
		if chrom not in G:
			G[chrom] 	= list()
		G[chrom].append((st, sp))
	M 		= {}
	for chrom in G:
		G[chrom].sort()
		i 	= 0
		N 	= len(G[chrom])
		if N:
			g 	= G[chrom]
			M[chrom] 	= list()
			while i < N:
				o_st,o_sp 	= g[i]
				while i < N and g[i][0] < o_sp and g[i][1] > o_st:
					o_st,o_sp 	= min(o_st, g[i][0]), max(o_sp, g[i][1])
					i+=1
				M[chrom].append((o_st,o_sp))
	return M




def get_overlaps(EXPS_M_D):
	M 	= merge_ALL(EXPS_M_D)
	DS 	= list()
	chroms 	= EXPS_M_D[0][1].keys()
	for i,(a_m, a_d, exp) in enumerate(EXPS_M_D):
		D=list()
		print i
		for chrom in chroms:
			if len([1 for a_m, A_d, exp in EXPS_M_D if chrom in A_d]) == len(EXPS_M_D):
				T 	= a_d[chrom]

				for start, stop in M[chrom]:
					if T.searchInterval((start,stop)):
						D.append(1)
					else:
						D.append(0)

		DS.append(D)


	FHW 	= open("/Users/joazofeifa/Lab/EMG/analysis/Distance.txt", 'w')
	for D in DS:
		for i,b in enumerate(D):
			if i + 1 < len(D):
				FHW.write(str(b) + ",")
			else:
				FHW.write(str(b) + "\n")
				
	FHW.close()
				
def display_dengraom(FILE,IDS):
	A = list()
	with open(FILE) as FH:
		for line in FH:
			D 	= [float(x) for x in line.strip("\n").split(",")]
			A.append(D)
	D = np.zeros((len(A), len(A)))
	print IDS
	print D.shape
	for i in range(len(A)):
		for j in range(len(A)):
			dist 	= sum([ abs(A[i][k] - A[j][k])  for k in range(len(A[i]))] )
			D[i,j] 	= dist+1
	linkage_matrix = linkage(D, "single")

	ddata = dendrogram(linkage_matrix,
                   color_threshold=1,
                   labels=IDS, orientation="right" )
	plt.tight_layout()
	plt.show()




def get_all(EXPS, ROOT):
	EXPS_M_D=list()
	for EXP in EXPS:
		DIR 	= ROOT+EXP+"/EMG_out_files/"
		M,D 	= load(DIR)
		EXPS_M_D.append((M,D, EXP))

	get_overlaps(EXPS_M_D)




if __name__ == "__main__":
	make=False
	EXPS2 	= ("Allen2014 (HCT116)", "Andersson2014 (HeLa)", "Core2014 (K562)", 
		"Duttke 2015 (HeLa)", "Jin2013 (IMR90)",  "Le2013 (MCF-7) ",
		 "Li2013 (IMR90)", 
		"Luo2014 (AC16)", "Puc2015 (LNCaP)", "Saponaro2014 (HeLa)" )
	EXPS 	= ("Allen2014", "Andersson2014", "Core2014", 
		"Duttke2015", "Jin2013",   "Le2013",
		 "Li2013", 
		"Luo2014", "Puc2015", "Saponaro2014" )
	if make:
		ROOT 	= "/Users/joazofeifa/Lab/gro_seq_files/"
		get_all(EXPS, ROOT)
	FILE 	= "/Users/joazofeifa/Lab/EMG/analysis/Distance.txt"
	display_dengraom(FILE, EXPS2)

