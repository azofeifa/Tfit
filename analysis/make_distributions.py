import numpy as np
import matplotlib.pyplot as plt
def loading(N):
	X 	= np.random.normal(0,1,N)
	plt.hist(X, bins=100, normed=1)
	plt.xlabel("Gene Location")
	plt.ylabel("Frequency")
	plt.grid()
	plt.show()

def loading_strand(N):
	X 	= np.random.normal(0,1,N)
	Y 	= np.random.normal(0,1,N)
	fct, fedges 	= np.histogram(X,bins=100, normed=1)
	rct, redges 	= np.histogram(Y,bins=100, normed=1)
	fedges 			= (fedges[1:] + fedges[:-1])/2.
	redges 			= (redges[1:] + redges[:-1])/2.
	
	plt.bar(fedges, fct, width=(fedges[-1]-fedges[0])/100)
	plt.bar(redges , -rct, width=(redges[-1]-redges[0])/100, color="red")
	plt.ylabel("Frequency")
	plt.grid()
	plt.show()

def loading_strand_exp(N):
	X 	= np.random.normal(0,1,N) + np.random.exponential(4,N)
	Y 	= np.random.normal(0,1,N) - np.random.exponential(4,N)
	fct, fedges 	= np.histogram(X,bins=100, normed=1)
	rct, redges 	= np.histogram(Y,bins=100, normed=1)
	fedges 			= (fedges[1:] + fedges[:-1])/2.
	redges 			= (redges[1:] + redges[:-1])/2.
	
	plt.bar(fedges, fct, width=(fedges[-1]-fedges[0])/100)
	plt.bar(redges , -rct, width=(redges[-1]-redges[0])/100, color="red")
	plt.ylabel("Frequency")
	plt.grid()
	plt.show()

def loading_strand_exp_UNI(N):
	X 	= np.random.normal(0,1,N) + np.random.exponential(4,N) 
	X  	= [x for x in X]	+ [ x for x in np.random.uniform(0,100, int(N*.95))]
	Y 	= np.random.normal(0,1,N) - np.random.exponential(4,N)
	fct, fedges 	= np.histogram(X,bins=100,normed=1 )
	rct, redges 	= np.histogram(Y,bins=100,normed=1 )
	fedges 			= (fedges[1:] + fedges[:-1])/2.
	redges 			= (redges[1:] + redges[:-1])/2.
	
	plt.bar(fedges, fct, width=(fedges[-1]-fedges[0])/100)
	plt.bar(redges , -rct, width=(redges[-1]-redges[0])/100, color="red")
	plt.ylabel("Frequency")
	plt.grid()
	plt.show()


N 	= 10000
#loading(N)
#loading_strand(N)
#loading_strand_exp(N)
loading_strand_exp_UNI(N)
