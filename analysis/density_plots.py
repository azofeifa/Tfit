import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import math 
def insert_bedgraph(A_L, FILES):
	for FILE in FILES:
		with open(FILE) as FH:
			print FILE
			prevchrom 	= ""
			ct 			= 0
			for line in FH:
				chrom,start, stop, cov = line.strip("\n").split("\t")
				if ct > pow(10,9):
					break
				ct+=1
						
				if chrom != prevchrom :
					if chrom in A_L:
						j,N 	= 0, len(A_L[chrom])
					else:
						j,N 	= 0,0
				start, stop 	= int(start), int(stop)
				while j<N and A_L[chrom][j][1] < start:
					j+=1
				if j < N and A_L[chrom][j][0]<stop:
					A_L[chrom][j][2].N+=(int(cov))
				prevchrom=chrom

def plot_density(overlaps):

	xy 	= [(a.N/float(a.stop - a.start),b.N/float(b.stop - b.start)) for a,b in overlaps if a.N and b.N ]
	x,y = [math.log(x,10) for x,y in xy ],[math.log(y,10) for x,y in xy ]
	p53 = [(a.N/float(a.stop - a.start),b.N/float(b.stop - b.start)) for a,b in overlaps if a.N and b.N and a.p53_site and b.p53_site ]
	print len(xy)
	xy 			= np.vstack([x,y])
	z 			= gaussian_kde(xy)(xy)
	F 			= plt.figure(figsize=(15,10))
	ax1 			= F.add_subplot(111)
	ax1.scatter(x, y, c=z, s=14, edgecolor='')
	ax1.scatter([math.log(x,10) for x,y in p53],[math.log(y,10) for x,y in p53], s=40, color="black")
	ax1.grid()
	ax1.set_xlabel("DMSO2_3")
	ax1.set_ylabel("Nutlin2_3")


	plt.show()



