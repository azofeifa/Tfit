import matplotlib.pyplot as plt
import numpy as np
def display(G):
	marks 		= "DNase","Pol-II S2", "H3K27ac","H3K4me1","H3K4me3"
	F 			= plt.figure()
	ax 			= F.add_subplot(1,1,1)
	h 			= 2
	delta  		= 1.25
	positions  	= np.linspace(1,len(G)*3, len(G)  )
	positions2 	= positions+delta
	ax.set_xlim(0,100)
	ax.barh(positions, [G[mark][0] for mark in marks[::-1]], height=h/2., color="red",edgecolor="red")
	ax.barh(positions2, [G[mark][1] for mark in marks[::-1]], height=h/2., color="blue",edgecolor="blue")
	ax.set_ylim(-1,len(G)*4)
	ax.set_yticks([])
	ax.set_xticklabels([int(x) for x in ax.get_xticks()], fontsize=(25))
	plt.show()

if __name__ == "__main__":
	G 	= {}
	G["DNase"] 		= (30.9, 94.9)
	G["H3K27ac"] 	= (79.1, 50.1)
	G["H3K4me1"] 	= (44.6, 41.2)
	G["H3K4me3"] 	= (66.6, 63.2)
	G["Pol-II S2"] 	= (66.6, 95.2)


	display(G)

