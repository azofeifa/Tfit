import matplotlib.pyplot as plt
import time
import numpy as np
import math
from scipy.stats import gaussian_kde

def match_UP(A_L, B_L):
	A,B ={}, {}

	for chrom in A_L:
		if chrom not in A:
			A[chrom] 	= list()
		for st, sp, S in A_L[chrom]:
			for m in S.models:
				st, sp 	= m.mu - m.std, m.mu+m.std
				A[chrom].append((st,sp, m))
	for chrom in B_L:
		if chrom not in B:
			B[chrom] 	= list()
		for st, sp, S in B_L[chrom]:
			for m in S.models:
				st, sp 	= m.mu - m.std, m.mu+m.std
				B[chrom].append((st,sp, m))
	overlaps 	= {}
	for chrom in A:
		overlaps[chrom] 	= list()
		if chrom in B:
			a,b 	= A[chrom],B[chrom]
			j,N 	= 0,len(b)
			for i,(a_st, a_sp, a_m) in enumerate(a):
				while j < N and b[j][1] < a_st:
					j+=1
				u 	= j
				o 	= list()
				while u < N and b[u][0] < a_sp:
					o_st,o_sp 	= max(b[u][0], a_st), min(b[u][1], a_sp)
					o.append(( o_sp-o_st, a[i],  b[u] ))
					u+=1
				if o:
					o 	= max(o)
					overlaps[chrom].append((o[1], o[2] ))
	return [(a[2], b[2]) for chrom in overlaps for a,b in overlaps[chrom] ]
				




def run(A_L,A_G,B_L, B_G):
	
	overlaps 	= match_UP(A_L, B_L)
	sigmas 		= [(x.si, y.si) for x,y in overlaps if x.wEM > 0.05 and y.wEM > 0.05]
	print len(sigmas)
	F 			= plt.figure(figsize=(15,10))
	ax1 		= F.add_subplot(121)
	ax2 		= F.add_subplot(122)
	x,y 		= [math.log(x,10) for x,y in sigmas],[math.log(y,10) for x,y in sigmas]
	ax1.scatter(x,y)
	ax1.grid()
	xy = np.vstack([x,y])
	
	z = gaussian_kde(xy)(xy)
	ax2.scatter(x, y, c=z, s=14, edgecolor='')
	ax2.grid()

	plt.show()



	pass