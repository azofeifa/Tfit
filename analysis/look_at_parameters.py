import matplotlib.pyplot as plt
import time
import numpy as np
import math
def display(G, U=True, N=True):
	F 	= plt.figure(figsize=(15,10))
	ax1 = F.add_subplot(2,2,1)
	ax1.set_title("sigma; loading variability")
	ax1.hist([math.sqrt(x) for x in G["N"][1] if x < 4000] ,bins=100,
		label="Mean: " + str(np.mean([math.sqrt(x) for x in G["N"][1] if x < 2000]))[:7],
		normed=1)
	ax1.legend()
	ax1.grid()
	
	ax2 = F.add_subplot(2,2,3)
	ax2.set_title("lambda; intiating rate")
	ax2.hist([x/100. for x in G["N"][2] if x < 2]  ,bins=100,  
	 label="Mean: " + str(np.mean([x/100. for x in G["N"][2] if x < 2]))[:7],
	 normed=1)
	ax2.legend()
	ax2.grid()
	
	ax3=F.add_subplot(2,2,2)
	ax3.set_title("pi; strand probability")
	ax3.hist([x for x in G["N"][4] if x > 0.01 and x < 0.99],bins=100, label="Mean: " + str(np.mean(G["N"][4]))[:7],
		normed=1)
	ax3.legend()
	ax3.grid()
	
	ax4=F.add_subplot(2,2,4)
	ax4.set_title("1/lambda; E[Y]")
	ax4.hist([100. / (x) for x in G["N"][2] if 100. / (x) < 3000 and x!= 2 ], bins=100,
		label="Mean: " + str(np.mean([100. / (x) for x in G["N"][2] if 100. / (x) < 3000 and x!= 2 ]))[:7],
		normed=1)
	ax4.legend()
	ax4.grid()


	plt.show()


def run(fits, spec=None, 
		weight_thresh=0., retry_tresh=0.,
		converged=True
		):
	G 	= { "N": (list(), list(), list(), list(), list()), 
			"U": (list(), list()) }
	
	for I in fits:
			for k, model in zip(I.models.keys(), I.models.values()):
				if ((spec is None or k == spec) and 
					model.retries <= retry_tresh 
					and (not converged or model.converged)):
					for rv in model.rvs:
						if rv.w > weight_thresh:
							if rv.type == "U":
								G["U"][0].append(rv.w)
								G["U"][0].append(rv.pi)
							elif rv.type=="N":
								G["N"][0].append(rv.mu)
								G["N"][1].append(rv.si*100)
								G["N"][2].append(rv.l)
								G["N"][3].append(rv.w)

								G["N"][4].append(rv.pi)
						
	display(G)
