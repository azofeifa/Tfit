import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math as m
from scipy.special import erf, erfc,gamma, gammainc
from scipy.stats import chi2
class RV:
	def __init__(self, sigma,lam,fp,w):
		self.si,self.l, self.fp, self.w 	= sigma, lam,fp, w
	def IN(self, x):
		return m.exp(-pow(x,2)*0.5)/m.sqrt(2*m.pi)
	def IC(self, x):
		return 0.5*(1+erf(x/m.sqrt(2.)))
	def R(self, x):
		if x > 5:
			return 1.0 / x
		N,D 	= self.IC(x), self.IN(x)
		if D < m.pow(10,-15): #python machine epsilon
			return 1.0 / m.pow(10,-15)
		return m.exp(m.log(1. - N)-m.log(D))

	def pdf(self,z,s,mu=0,l=100):
		U 	= (1.0-self.w) / l
		if s == 1:
			z-=self.fp
		else:
			z+=self.fp
		vl 		= (self.l /2.)* (s*2*(mu-z) + self.l*self.si**2. )
		p 	= self.l*self.IN((z-mu)/self.si)*self.R(self.l*self.si - s*((z - mu)/self.si))
		E 	= self.w*p
		return m.log(E+U)

def get_ll(D,rv,mu=0,l=100):
	return sum([ rv.pdf(D[i,0],1,mu=mu,l=l)*D[i,1] for i in range(len(D ))]) + sum([ rv.pdf(D[i,0],-1,mu=mu,l=l)*D[i,2] for i in range(len(D))])

def get_pvalue_chi_seq(x, k):
	CHI 	= chi2(k)
	vl1 	= gammainc(k,x)
	vl2 	= gamma(k)

	S 		= (x - k)/m.sqrt(2*k)

	pv 		= 0.5*(1+erf(S/m.sqrt(2.)))
	return pv
def across_segment(forward_bg, reverse_bg, chrom,start, stop, sigma=10,lam=200,fp=10,pi=0.5,w=0.6,ns=100.0):
	D 	= [list(),list()]
	for i,f in enumerate((forward_bg, reverse_bg)):
		FH 	= open(f,'r')
		for line in FH:
			schrom,st, sp, cov 	= line.strip("\n").split("\t")
			x 					= ( float(st) + float(sp) ) / 2.
			if schrom == chrom and start < x < stop:
				D[i].append((x,float(cov)))
			if len(D[i]) and x > stop:
				break
		FH.close()
	st,sp 	= min(min(D[0]),min(D[1]))[0],max(max(D[0]),max(D[1]))[0]
	bins 	= int((sp-st)/15.)
	X 		= np.zeros((bins,3))
	countsf = np.histogram([x for x,y in D[0]],range=(st,sp), bins=bins,weights=[y for x,y in D[0]])[0]
	countsr = np.histogram([x for x,y in D[1]],range=(st,sp), bins=bins,weights=[y for x,y in D[1]])[0]

	X[:,0] 	= np.linspace(st,sp,bins)
	X[:,0]-=X[0,0]
	X[:,0]/=100.0
	X[:,1] 	= countsf
	X[:,2] 	= countsr
	F 		= plt.figure()
	ax1,ax2 = F.add_subplot(3,1,1),F.add_subplot(3,1,2)
	ax3 	= F.add_subplot(3,1,3)
	ax1.bar(X[:,0],X[:,1],width=(X[-1,0]-X[0,0])/len(X))
	ax1.bar(X[:,0],-X[:,2],width=(X[-1,0]-X[0,0])/len(X))

	#================================================================
	sigma/=ns
	fp/=ns
	lam/=ns
	lam=1.0/lam

	rv 		= RV(sigma,lam,fp,w)
	rv2 	= RV(sigma,lam,fp,0.0)

	lls 	= list()
	W 		= 1.5 #100,bp
	
	W 		= (X[-1,0]-X[0,0])/W

	for i in range(0,len(X)):
		st,sp 	= int(max(0,i-W)),int(min(len(X), i+W))

		nD 		= X[st:sp,:]
		NN 		= np.sum(nD[:,1:])
		ll 		= get_ll(nD,rv,mu 	= X[i,0],l=W)
		ll2 	= get_ll(nD,rv2,mu 	= X[i,0],l=W)
		lls.append((2*(ll-ll2), NN ))

	ax2.plot(range(0,len(X)), [x[0] for x in lls])
	ax3.plot(range(0,len(X)), [ get_pvalue_chi_seq(x[0],150) for x in lls])
	plt.show()

def make_rand(N=1000,sigma=10,lam=200,fp=10,pi=0.5,w=0.6,ns=100.0):
	def make_sample(n=100,d=1000):
		bins 	= int(d/15.)
		D 		= [list(),list()]
		D[0] 	= np.random.randint(0,int(d),n)
		D[1] 	= np.random.randint(0,int(d),n)
		X 		= np.zeros((bins,3))
		countsf = np.histogram(D[0],range=(0,d), bins=bins )[0]
		countsr = np.histogram(D[1],range=(0,d), bins=bins )[0]

		X[:,0] 	= np.linspace(0,d,bins)
		X[:,0]-=X[0,0]
		X[:,0]/=100.0
		X[:,1] 	= countsf
		X[:,2] 	= countsr
		return X
	sigma/=ns
	fp/=ns
	lam/=ns
	lam=1.0/lam

	rv 		= RV(sigma,lam,fp,w)
	rv2 	= RV(sigma,lam,fp,0.0)
	lls 	= list()
	for i in range(N):
		X 		= make_sample(n=100)
		ll 		= get_ll(X,rv,mu 	= X[X.shape[0]/2,0],l=X[-1,0]-X[0,0])
		ll2 	= get_ll(X,rv2,mu 	= X[X.shape[0]/2,0],l=X[-1,0]-X[0,0])

		lls.append( 2*(ll-ll2) )
	plt.hist(lls,bins=30)
	plt.show()




if __name__ == "__main__":
	forward_bg 			= "/Users/joazofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/test_DMSO2_3.pos.BedGraph"
	reverse_bg 			= "/Users/joazofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/test_DMSO2_3.neg.BedGraph"

	chrom,start, stop 	= "chr1",440747,444636
	across_segment(forward_bg, reverse_bg,chrom,start, stop )
	#make_rand()