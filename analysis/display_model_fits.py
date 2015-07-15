import matplotlib.pyplot as plt
import numpy as np
def display(fits, bins=200):
	for fit in fits:
		fit.BIN(bins)
		F 	= plt.figure(figsize=(15,10))
		fit.X[:,0]-=fit.X[0,0]
		fit.X[:,0] /=100.
			
		for i,k in enumerate(fit.models):
			ax 	= F.add_subplot(2,2,i)
			N 	= np.sum(fit.X[:,1:])
			xs 	= np.linspace(fit.X[0,0], fit.X[-1,0], 1000)
			maxY = max(map(lambda x: fit.models[k].pdf(x, -1) , xs ))
			maxY = 0.25
			ax.bar(fit.X[:,0], (fit.X[:,1]*maxY)/N, width=(fit.X[-1,0]-fit.X[0,0])/bins, alpha=0.25)
			ax.bar(fit.X[:,0], -(fit.X[:,2]*maxY)/N, width=(fit.X[-1,0]-fit.X[0,0])/bins, alpha=0.25)
			ax.plot(xs, map(lambda x: fit.models[k].pdf(x, 1) , xs ), linewidth=2.)
			ax.plot(xs, map(lambda x: -fit.models[k].pdf(x, -1) , xs ), linewidth=2.)
			
			ax.grid()
		plt.show()
