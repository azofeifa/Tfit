import load
import matplotlib.pyplot as plt
import numpy as np
import math
class NORM:
	def __init__(self,mu,si,w):
		self.mu 	= mu
		self.si 	= si
		self.w 		= w
	def pdf(self, x):
		return (self.w/(self.si*math.sqrt(2*math.pi)))*math.exp(-pow(x-self.mu,2)/(2*pow(self.si,2)))
class UNI:
	def __init__(self,a,b,w):
		self.a,self.b,self.w 	= a,b,w
	def pdf(self, x):
		return (self.w/abs(self.a-self.b))
def draw(rvs, X):
	F 			= plt.figure(figsize=(15,10))
	ax 				= F.add_subplot(111)
	counts,edges 	= np.histogram(X[:,0],weights=X[:,1], normed=1,bins=150)
	
	ax.bar(edges[1:],  counts  , color="blue", edgecolor="blue", alpha=0.5, width=0.25*(X[-1,0]-X[0,0])/150)
	xs 			= np.linspace(X[0,0], X[-1,0], 1000)
	ys_forward 	= map(lambda x: sum([rv.pdf(x) for rv in rvs     ]) , xs) 
	
	ax.plot(xs, ys_forward, linewidth=2.5,  color="black")
	ax.set_xticklabels([int(i) for i in ax.get_xticks()], fontsize=20)
	ax.set_yticklabels([i for i in ax.get_yticks()], fontsize=20)
	plt.show()

def fit(X,K):
	min_x 	= min(X[:,0])
	max_x 	= max(X[:,0])
	rvs 	= [NORM(np.random.uniform(min_x, max_x), 5.0,1.0/(k+1)   )  for k in range(K)]+[UNI(min_x, max_x, 0.05)]
	
	T 		= 100
	t 		= 0
	while t < T:

		EX 		= np.zeros((K,))
		EX2 	= np.zeros((K,))
		W 		= np.zeros((K,))
		for i in range(X.shape[0]):
			x 	= X[i,0]
			y 	= X[i,1]
			rs 	= [rv.pdf(x) for rv in rvs ]
			rs 	= [r/sum(rs) for r in rs ]
			for k in range(K):
				EX[k]+=x*rs[k]*y
				EX2[k]+=pow(rvs[k].mu-x,2)*rs[k]*y
				W[k]+=rs[k]*y
		weights 	= [W[k]/sum(W) for k in range(K)]
		for k in range(K):
			rvs[k].mu 	= EX[k]/W[k]
			rvs[k].si 	= math.sqrt(EX2[k]/W[k])
			rvs[k].w 	= weights[k]
		draw(rvs, X )
		t+=1




if __name__ == "__main__":
	WRITE 	= False
	FILE_NAME 	= "/Users/joazofeifa/DNAse_region_of_interest.bed"
	if WRITE:
	 	X 		=  load.grab_specific_region("chr1",212730537, 212742051, 
				pos_file="/Users/joazofeifa/Lab/ChIP/HCT116/DNase/DNAse.bedgraph", 
				neg_file="/Users/joazofeifa/Lab/ChIP/HCT116/DNase/DNAse.bedgraph",
				SHOW 	=False, bins=1000)
	 	X[:,0]-=X[0,0]
		X[:,0]/=100.
		FHW 	= open(FILE_NAME, "w")
		for i in range(X.shape[0]):
			FHW.write(str(X[i,0]) + "\t" + str(X[i,1]) + "\t" + str(X[i,2]) + "\n")
		FHW.close()
	X 	= list()

	with open(FILE_NAME) as FH:
		for line in FH:
			x,y,z 	= [float(x) for x in line.strip("\n").split("\t")]
			X.append([x,y,z])
	X 	= np.array(X)
	fit(X,4)
