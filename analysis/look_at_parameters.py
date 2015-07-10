import matplotlib.pyplot as plt
import time
import numpy as np
import math
import scipy.stats as ss


def display(G, U=True, N=True):
	F 	= plt.figure(figsize=(15,10))
	ax1 = F.add_subplot(2,2,1)
	ax1.set_title("sigma; loading variability")
	data 	= [math.sqrt(x*100) for x in G["N"][1] if x < 20]
	fit_alpha, fit_loc, fit_beta = ss.gamma.fit(data)
	
	rv 							= ss.gamma(fit_alpha, fit_loc, fit_beta)
	xs 							= np.linspace(0, max(data) , 1000)

	ax1.hist(data ,bins=100,
		label="Mean: " + str(np.mean(data))[:7],
		normed=1)
	ax1.plot(xs, map(rv.pdf , xs), linewidth=2,label="alpha: " + str(fit_alpha)[:5] + 
		"\nbeta: " + str(fit_beta)[:5] + "\nloc: " + str(fit_loc)[:5])
	ax1.legend()
	ax1.grid()
	
	ax2 = F.add_subplot(2,2,3)
	ax2.set_title("lambda; intiating rate")

	data 						= [x/100. for x in G["N"][2] if x < 2]

	fit_alpha, fit_loc, fit_beta = ss.gamma.fit(data)
	
	rv 							= ss.gamma(fit_alpha, fit_loc, fit_beta)
	xs 							= np.linspace(0, max(data) , 1000)

	ax2.hist(data  ,bins=100,  
	 label="Mean: " + str(np.mean([x/100. for x in G["N"][2] if x < 2]))[:7],
	 normed=1)
	ax2.plot(xs, map(rv.pdf , xs), linewidth=2,label="alpha: " + str(fit_alpha)[:5] + 
		"\nbeta: " + str(fit_beta)[:5] + "\nloc: " + str(fit_loc)[:5])
	ax2.legend()
	ax2.grid()
	
	ax3=F.add_subplot(2,2,2)
	ax3.set_title("pi; strand probability")
	data 						= [x for x in G["N"][4] if x > 0.01 and x < 0.99]

	a,b,loc, scale  = ss.beta.fit(data)
	
	rv 							= ss.beta(a,b,loc, scale)
	xs 							= np.linspace(0, max(data) , 1000)



	ax3.hist(data,bins=100, label="Mean: " + str(np.mean(data))[:7],
		normed=1)
	ax3.plot(xs, map(rv.pdf, xs), 
		linewidth=3., label="alpha: " + str(a)[:5] + "\nbeta: " + str(b)[:5] + "\nloc: " + str(loc)[:5]
		+"\nscale: " + str(scale)[:5])
	ax3.legend()
	ax3.grid()
	
	ax4=F.add_subplot(2,2,4)
	ax4.set_title("Length of Single Isoform Genes")
	data 	= [x/100. for x in G["U"][2] if x < 150000]
	fit_alpha, fit_loc, fit_beta = ss.gamma.fit(data)
	
	rv 							= ss.gamma(fit_alpha, fit_loc, fit_beta)
	xs 							= np.linspace(0, max(data) , 1000)
	
	ax4.hist(data, bins=100, normed=1)
	ax4.plot(xs, map(rv.pdf , xs), linewidth=2,label="alpha: " + str(fit_alpha)[:7] + 
		"\nbeta: " + str(fit_beta)[:7] + "\nloc: " + str(fit_loc)[:7])
	ax4.grid()
	ax4.legend()
	plt.tight_layout()
	plt.show()


def run(fits, spec=None, 
		weight_thresh=0., retry_tresh=0.,
		converged=True
		):
	G 	= { "N": (list(), list(), list(), list(), list()), 
			"U": (list(), list(), list()) }
	
	for I in fits:
			for k, model in zip(I.models.keys(), I.models.values()):
				if ((spec is None or k == spec) and 
					model.retries <= retry_tresh 
					and (not converged or model.converged)):
					for rv in model.rvs:
						if rv.w > weight_thresh:
							if rv.type == "U":
								G["U"][0].append(rv.w)
								G["U"][1].append(rv.pi)
								G["U"][2].append(I.stop-I.start)
							elif rv.type=="N":
								G["N"][0].append(rv.mu)
								G["N"][1].append(rv.si)
								G["N"][2].append(rv.l)
								G["N"][3].append(rv.w)

								G["N"][4].append(rv.pi)
						
	display(G)
