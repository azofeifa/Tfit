import matplotlib.pyplot as plt
import numpy as np
def wall_vs_np(G,A):
	data 	= [(NP, np.mean(G[NP][1]) ) for NP in G]
	CPU 	= [(NP, np.mean(G[NP][0]) ) for NP in G]
	std 	= [np.std(G[NP][1]) for NP in G]
	c_std 	= [np.std(G[NP][0]) for NP in G]
	MAX 	= max([y for x,y in data])
	MIN 	= min([y for x,y in data])
	data.sort()
	F 		= plt.figure(figsize=(15,10))
	ax 		= F.add_subplot(2,2,1)
	ax.set_title("Wall Time vs Number of Processors (C++ OpenMP), pando")
	ax.scatter([x for x,y in data],[y for x,y in data],
		label="Average Wall Time\n(5 Runs Each)\n" + "Max: " + 
		str(MAX) + "\nMin: " + str(MIN) )
	ax.plot([x for x,y in data],[y for x,y in data])
	ax.fill_between([x for x,y in data],
		[y + std[i] for i,(x,y) in enumerate(data)], 
		[y - std[i] for i,(x,y) in enumerate(data)], 
		alpha=0.5, color="grey")
	ax.set_xlabel("Number of Processors")
	ax.set_ylabel("Wall Time (Seconds)")
	ax.legend()
	
	ax.grid()
	ax2 		= F.add_subplot(2,2,3)
	ax2.set_title("CPU Time vs Number of Processors (C++ OpenMP), pando")
	ax2.scatter([x for x,y in CPU ],[y for x,y in CPU ], 
		label="Average CPU Time\n(5 Runs Each)\nGrey Shading:\n(One Standard Deviation)")
	ax2.plot([x for x,y in CPU ],[y for x,y in CPU ])
	ax2.fill_between([x for x,y in CPU],
		[y + c_std[i] for i,(x,y) in enumerate(CPU)], 
		[y - c_std[i] for i,(x,y) in enumerate(CPU)], 
		alpha=0.5, color="grey")
	ax2.set_xlabel("Number of Processors")
	ax2.set_ylabel("CPU Ticks (Seconds)")
	ax2.grid()
	ax2.legend(loc=(0.02,0.69))

	ax3 		= F.add_subplot(2,2,2)

	data 	= [(NP, np.mean(A[NP][1]) ) for NP in A]
	CPU 	= [(NP, np.mean(A[NP][0]) ) for NP in A]
	std 	= [np.std(A[NP][1]) for NP in A]
	c_std 	= [np.std(A[NP][0]) for NP in A]
	MAX 	= max([y for x,y in data])
	MIN 	= min([y for x,y in data])


	data.sort()
	ax3.set_title("Wall Time vs Number of Processors (C++ OpenMP), vieques")
	ax3.scatter([x for x,y in data],[y for x,y in data],
		label="Average Wall Time\n(5 Runs Each)\n" + "Max: " + 
		str(MAX) + "\nMin: " + str(MIN) )
	ax3.plot([x for x,y in data],[y for x,y in data])
	ax3.fill_between([x for x,y in data],
		[y + std[i] for i,(x,y) in enumerate(data)], 
		[y - std[i] for i,(x,y) in enumerate(data)], 
		alpha=0.5, color="grey")
	ax3.set_xlabel("Number of Processors")
	ax3.set_ylabel("Wall Time (Seconds)")
	ax3.set_ylim(0, MAX+1000)
	ax3.grid()
	ax3.legend()
	
	ax4 		= F.add_subplot(2,2,4)
	ax4.set_title("CPU Time vs Number of Processors (C++ OpenMP), pando")
	ax4.scatter([x for x,y in CPU ],[y for x,y in CPU ], 
		label="Average CPU Time\n(5 Runs Each)\nGrey Shading:\n(One Standard Deviation)")
	ax4.plot([x for x,y in CPU ],[y for x,y in CPU ])
	ax4.fill_between([x for x,y in CPU],
		[y + c_std[i] for i,(x,y) in enumerate(CPU)], 
		[y - c_std[i] for i,(x,y) in enumerate(CPU)], 
		alpha=0.5, color="grey")
	ax4.set_xlabel("Number of Processors")
	ax4.set_ylabel("CPU Ticks (Seconds)")
	ax4.grid()
	ax4.legend(loc=(0.02,0.69))


	plt.show()
	pass

def delta_ll_vs_RI(G):
	I 	= {}
	J 	= {}
	L 	= {}
	N 	= 0.
	for gene in G:
		G[gene].calc_improvements(TYPE=1)
		data 	= G[gene].get_average_diff(TYPE=1)
		G[gene].calc_improvements(TYPE=2)
		data2 	= G[gene].get_average_diff(TYPE=2)
		G[gene].calc_improvements2()
		data3 	= G[gene].get_average_diff2()

		
		for rounds, comparisons in data:
			if rounds not in I:
				I[rounds] 	= list()
			I[rounds].append(comparisons)
		for rounds, comparisons in data2:
			if rounds not in J:
				J[rounds] 	= list()
			J[rounds].append(comparisons)
		for rounds, extra in data3:
			if rounds not in L:
				L[rounds] 	= {}
			for move, diff in extra:
				if move not in L[rounds]:
					L[rounds][move] 	= list()
				L[rounds][move].append(diff)



	
	F 	= plt.figure(figsize=(15,10))
	ax1 	= F.add_subplot(2,1,1)
	ax1.scatter(I.keys(), map(np.mean, I.values()))
	ax1.scatter(J.keys(), map(np.mean, J.values()), color="green")

	ax1.grid()

	ax2 	= F.add_subplot(2,1,2)
	ax2.scatter(L.keys(), map(np.mean, [L[r][0] for r in L ]))
	ax2.scatter(L.keys(), map(np.mean, [L[r][1] for r in L ]), color="red")
	ax2.scatter(L.keys(), map(np.mean, [L[r][2] for r in L ]), color="green")
	ax2.scatter(L.keys(), map(np.mean, [L[r][3] for r in L ]), color="black")
	
	ax2.grid()
	

	plt.show()





	pass


def hybrid(G):
	F 		= plt.figure(figsize=(15,10))
	ax 		= F.add_subplot(121)
	means  	= [np.mean([ y.total[1]/60.  for x,y in G[node] ]) for node in G]
	std 	= [np.std([ y.total[1]/60.  for x,y in G[node] ]) for node in G]
	ax.set_title("Linear Plot")
	ax.scatter(G.keys(), means)
	ax.plot(G.keys(), means)
	ax.fill_between(G.keys(), [ m-std[i] for i,m in enumerate(means) ], [ m+std[i] for i,m in enumerate(means) ], alpha=0.5, color="grey")
	ax.set_ylabel("Minutes [CPU/WT]")
	ax.set_xlabel("MPI Processes (nodes)")

	ax.grid()
	ax2 		= F.add_subplot(122)
	ax2.set_title("Log-Linear Plot")
	ax2.scatter(G.keys(), means)
	ax2.plot(G.keys(), means)
	ax2.fill_between(G.keys(), [ m-std[i] for i,m in enumerate(means) ], [ m+std[i] for i,m in enumerate(means) ], alpha=0.5, color="grey")
	ax2.grid()
	ax2.set_ylabel("Minutes [CPU/WT] (log space)")
	ax2.set_xlabel("MPI Processes (nodes)")
	ax2.set_yscale("log")
	plt.show()
	pass






