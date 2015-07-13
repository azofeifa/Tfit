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