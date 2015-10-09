import random as r
import matplotlib.pyplot as plt
import numpy as np
import math as m
import time
def resample(X):
	pass
def norm(x, mu, si):
	return 1.0 / (si*m.sqrt(2*m.pi))*m.exp(-pow(x-mu,2)/(2*pow(si,2)) )
def sample(edges, counts):
	N 		= float(sum(counts))
	cdf 	= [ sum(counts[:i+1]) / N  for i in range(0, len(counts)) ]
	pdf 	= [ counts[i]/ N  for i in range(0, len(counts)) ]
	nN 		= 0
	X 		= np.zeros((len(edges), 2))
	X[:,0] 	= edges
	while nN < N:
		U 	= r.uniform(0,1)
		j 	= 0
		while (j+ 1 < len(counts) and cdf[j] < U):
			j+=1
		if (counts[j]!= 0):
			ct 	= np.random.geometric(1-pdf[j])
			X[j,1]+=1
			nN+=1
	return X
def sample_geom(edges, counts):
	N 		= float(sum(counts))
	cdf 	= [ sum(counts[:i+1]) / N  for i in range(0, len(counts)) ]
	pdf 	= [ counts[i]/ N  for i in range(0, len(counts)) ]
	nN 		= 0
	X 		= np.zeros((len(edges), 2))
	X[:,0] 	= edges
	while nN < N:
		U 	= r.uniform(0,1)
		for j in range(0, len(counts)):
			if (counts[j]!= 0):
				ct 	= np.random.geometric(pdf[j])
				X[j,1]+=ct
				nN+=ct
	return X


def get_mean(X):
	N 	= sum(X[:,1])
	return sum([X[i,0]*X[i,1] for i in range(X.shape[0])]) / N
def compare_means(edges, counts):

	start 			= time.clock()
	means_geom 		= [abs(0-get_mean(sample_geom(edges, counts) )) for i in range(50)]
	print "geometric: ", time.clock() - start
	start 			= time.clock()
	means_reg 		= [abs(0 - get_mean(sample(edges, counts) )) for i in range(50)]
	print "reg: ", time.clock() - start
	F 				= plt.figure(figsize=(15,10))
	ax 				= F.add_subplot(1,1,1)
	bp 				= ax.boxplot((means_reg, means_geom),patch_artist=True)
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
	ax.set_xticklabels(["Correct Bootstrap", "Approx Bootstrap" ])
	ax.grid()
	plt.show()



if __name__ == "__main__":
	R 				= np.random.normal(0.,1, 1000)
	bins 			= 100
	counts,edges  	= np.histogram(R, bins=bins)
	width 			= (max(edges) - min(edges)) / bins
	edges 			= (edges[1:] + edges[:-1])/2.
	compare_means(edges, counts)