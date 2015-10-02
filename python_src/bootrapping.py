import random as r
import matplotlib.pyplot as plt
import numpy as np
import math as m
def resample(X):
	pass
def norm(x, mu, si):
	return 1.0 / (si*m.sqrt(2*m.pi))*m.exp(-pow(x-mu,2)/(2*pow(si,2)) )


if __name__ == "__main__":
	R 	= np.random.normal(0.01,0.000001, 1000)
	xs 	= np.linspace(min(R), max(R), 100)
	mu 	= np.mean(R)
	si 	= np.std(R)
	print mu-si
	ys 	= [norm(x,mu,si) for x in xs]
	plt.hist(R, normed=1, bins=35)
	plt.plot(xs,ys)
	plt.show()
	pass