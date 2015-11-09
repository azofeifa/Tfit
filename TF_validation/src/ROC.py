import load
import matplotlib.pyplot as plt
import time
import math 
import numpy as np

from scipy.special import erf
def choose(n,k):
	if k > n:
		return 0
	r=1
	N=0
	for d in range(1,int(k+1)):
		r*=(n-N)
		r/=d
		N+=1
	return r

def normal(x,mu, si):
	return (1.0 / (math.sqrt(2*math.pi)*si  ))*math.exp(-pow(x-mu,2)/(2*pow(si,2)))
def normal_cdf(x,mu,si):
	return 0.5*(1+erf((x-mu)/(si*math.sqrt(2)) ))
def binomial(k,w,l,M):
	p 	= float(2*w) / float(l- 2*w)
	vl 	= choose(M,k)*pow(p,k)*pow(1.0-p,M-k  )
	return vl


def draw_background(ax,l, N):
	ks 	= range(1, 3000)
	w 	= 1500 
	p 	= float(2*w) / float(l- 2*w)
	mu 	= N*p

	si  = math.sqrt(mu*(1-p))

	pdfs 	= np.array([normal(k,mu, si) for k in ks])
	Spdf 	= sum(pdfs)
	ax.plot(ks, [ sum(pdfs[:i])/ Spdf for i in range(len(pdfs)) ] )
	ax.set_xscale("log")
class scanner:
	def __init__(self, chrom,start, stop, params, CHIP, sc):
		self.chrom,self.start, self.stop=chrom,start, stop
		self.center 	= (self.start + self.stop) / 2.
		self.CHIP, self.sc 	="",""
		if CHIP:
			self.CHIP 	= dict([ (x.split(":")[0],float(x.split(":")[1] )) for x in CHIP.split(",")])
		if sc:
			self.sc 	= [ [ float(y) if i < 2 else y for i,y in enumerate(x.split(":")) ] for x in sc.split(",")]
			for i in range(len(self.sc)):
				self.sc[i][1]=(self.sc[i][1]+self.start) - self.center
		mu,si,l,w,pi,N,fp,nll 	= [float(x) for x in params.split("_")]
		self.mu,self.si,self.l,self.w,self.pi,self.N,self.fp,self.ll 	= mu,si,l,w,pi,N,fp,nll
	def get_density(self):
		return self.N / (self.stop - self.start)
	
	def get_max_motif(self, D=500):
		results 	= [ ll for ll,i,strand in self.sc if abs(i) < D]
		if results:
			return max(results)
		return None
	def get_min_distance(self, D=1500, PSSM=-13):
		results 	= [ abs(i) for ll,i,strand in self.sc if abs(i) < D and ll > PSSM]
		if results:
			return min(results)
		return None
	def is_hit(self, ChIP, mu, si):
		if ChIP in self.CHIP:
			return int(normal_cdf(self.CHIP[ChIP], mu,si)   )
		return 0
	def get_BIC_ratio(self, penality=7, window=2000):
		A 	=   self.ll / self.N
		return  A 



def load_scanner_out(FILE,test=True):
	S 	= list()
	i 	= 0
	with open(FILE) as FH:
		for line in FH:
			i+=1
			if test and i > 2500:
				break
			header, scanner_info 		= line.strip("\n").split("\t")
			chrom,stsp, params, CHIP 	= header.split("|")
			start, stop 				= [float(x) for x in stsp.split("-")]
			s 							= scanner(chrom,start, stop, params, CHIP, scanner_info)
			S.append(s)
	return S

def ROC(S):
	l 	= 2951343224.0
	n 	= 1417281200.0
	w 	= 1500.0 
	p 	= float(2*w) / float(l- 2*w)
	mu 	= n*p

	si  = math.sqrt(mu*(1-p))
	si 	= 10
	print mu, math.log(mu,10)
	F 	= plt.figure()
	ax 	= F.add_subplot(1,2,1)
	ax2 = F.add_subplot(1,2,2)
	ax.set_title("Significant ChIP Peaks\nhave a smaller motif-i distance")
	
	ax.hist([ s.get_min_distance() for s  in S if s.is_hit("SP1",mu,si) if s.get_density()>3 and s.get_min_distance() is not None  ],alpha=0.4, bins=100, label="ChIP, " + r'$p<10^{-4}$' )
	ax.hist([ s.get_min_distance() for s  in S if not s.is_hit("SP1",mu,si) if s.get_density()>2 and s.get_min_distance() is not None ],alpha=0.4, bins=100, label="ChIP, " + r'$p>10^{-4}$' )
	
	# ax.hist([math.log(s.CHIP["SP1"]+1,10) for s in S if "SP1" in s.CHIP and s.is_hit("SP1", mu,si ) and s.CHIP and math.log(s.CHIP["SP1"]+1)  < 10], alpha=0.5, label="sig")
	# ax.hist([math.log(s.CHIP["SP1"]+1,10) for s in S if "SP1" in s.CHIP and not s.is_hit("SP1", mu,si  ) and  s.CHIP and math.log(s.CHIP["SP1"]+1)  < 10 ],alpha=0.5, label="not sig")

	#ax.set_xlim(0,1400)
	ax.legend()
	ax.grid()
	S 				= [(s,  s.is_hit("SP1", mu, si)) for s in S if s.get_min_distance(D=1500, PSSM=-9 )  ]
 	penality  		= np.linspace(0,1000)
 	P 				= float(len([ 1 for s,z in S if z    ]))
 	N 				= float(len([1 for s,z in S if not z  ]))
 	TPS,FPS 		= list(),list()
	##======
	for p in penality:
		TP,TN 	= 0.0,0.0
		for s, z in S:

			d 	= s.get_min_distance(D=p, PSSM=-9 )
			if d is None and not z:
				TN+=1
			elif d and z:
				TP+=1

		TPS.append(TP/P)
		FPS.append(1-(TN / N))
		print TPS[-1], FPS[-1]

	ax2.plot(FPS, TPS)

	#ax2.set_xlim(0,1.1)
	#ax2.set_ylim(0,1.1)


	plt.show()






if __name__ == "__main__":
	FILE = "/Users/joazofeifa/Lab/TF_predictions/HOMMER_OUT_FILES/sp1_sites.bed"
	S 	= load_scanner_out(FILE)
	ROC(S)