import math 
import numpy as np
import matplotlib.pyplot as plt
import time
def bic_function(ll, n, K, penality):
	return -2.*ll + K*4*math.log(n)*penality

def run(fits):
	resolution 	= 100
	penalities 	= np.linspace(0.01, 500, resolution)
	errors 		= np.zeros((resolution,))
	N_ps 	= list()
	N_ss 	= list()
	N_ls 	= list()
	N_ws 	= list()
	
	for i,p in enumerate(penalities):
		E 		= 0.
		c 		= list()
		c1 		= list()
		c2 		= list()
		c3 		= list()
		for I in fits:
			curr 	= np.inf
			argCurr = None
			for K in I.models:
				val = bic_function(I.models[K].ll, I.N, K, p)
				if val < curr:
					argCurr 	= K
					curr 	 	= val
			E += 3.-argCurr
			d 	= [rv.pi for rv in I.models[argCurr].rvs if rv.type=="N"]
			d2 	= [rv.si for rv in I.models[argCurr].rvs if rv.type=="N"]
			d3 	= [rv.l for rv in I.models[argCurr].rvs if rv.type=="N"]
			d4 	= [rv.w for rv in I.models[argCurr].rvs if rv.type=="N"]
			
			if d:
				c.append(np.mean(d))
				c1.append(np.mean(d2))
				c2.append(np.mean(d3))
				c3.append(np.mean(d4))

		errors[i] 	= (E / len(fits))
		N_ps.append(np.mean(c))
		N_ss.append(np.mean(c1))
		N_ls.append(np.mean(c2))
		N_ws.append(np.mean(c3))
		
	i 	= 0
	errors.sort()
	while i < len(errors) and errors[i]<0:
		i+=1
	p 	= penalities[i]
	KS 	= list()
	for I in fits:
		curr 	= np.inf
		argCurr = None
		for K in I.models:
			val = bic_function(I.models[K].ll, I.N, K, p)
			if val < curr:
				argCurr 	= K
				curr 	 	= val
		E += 3.-argCurr
		KS.append(argCurr)
		
	F 	= plt.figure(figsize=(15,10))
	ax 	= F.add_subplot(2,2,1)
	ax.scatter(penalities, errors)
	ax.set_title("BIC Penality vs Accuracy in Picking Single Isoform Gene")
	ax.set_xlabel("BIC Penality")
	ax.set_ylabel("Erorr")
	ax.grid()
	ax2 	= F.add_subplot(2,2,2)
	ax2.scatter(penalities, N_ps)
	ax2.scatter(penalities, N_ws, color="green")

	ax2.set_title("BIC Penality vs W and PI Parameters")
	ax2.set_xlabel("BIC Penality")
	ax2.set_ylabel("EMG Average Weight\nand Strand Probability")
	
	ax2.grid()
	ax3 	= F.add_subplot(2,2,3)
	ax3.scatter(penalities, N_ss)
	ax3.set_title("BIC Penality vs Sigma (variance)")
	ax3.set_xlabel("BIC Penality")
	ax3.set_ylabel("EMG Average Loading Variance ")
	ax3.grid()

	ax4 	= F.add_subplot(2,2,4)
	ax4.scatter(penalities, N_ls)
	ax4.set_title("BIC Penality vs Lambda (Initiating)")
	ax4.set_xlabel("BIC Penality ")
	ax4.set_ylabel("EMG Average Initiating Rate ")
	
	ax4.grid()
	plt.tight_layout()
	
	plt.show()

