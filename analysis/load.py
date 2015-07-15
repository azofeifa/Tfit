import numpy as np
from scipy.special import erf, erfc
import math as m
class EMG:
	def __init__(self, mu, si, l, w, pi):
		self.mu, self.si, self.l 	= mu, si, l
		self.w, self.pi 			= w,pi
		self.type 					= "N"
	def IN(self, x):
		return m.exp(-pow(x,2)*0.5)/m.sqrt(2*m.pi)
	def IC(self, x):
		return 0.5*(1+erf(x/m.sqrt(2.)))

	def R(self, x):
		if x > 8:
			return 1.0 / x
		N,D 	= self.IC(x), self.IN(x)
		if D < m.pow(10,-15): #python machine epsilon
			return 1.0 / m.pow(10,-15)
		return m.exp(m.log(1. - N)-m.log(D))
	def pdf(self, z,s):
		vl 		= (self.l /2.)* (s*2*(self.mu-z) + self.l*self.si**2. )
		if vl > 200: #overflow error, try this one...
			p 	= self.l*self.IN((z-self.mu)/self.si)*self.R(self.l*self.si - s*((z - self.mu)/self.si))
		try:
			p 		= (self.l / 2.)*m.exp( vl )*erfc( (s*(self.mu- z) + self.l*self.si**2. ) / (m.sqrt(2)*self.si) )
		except:
			self.remove,p 	= True,0.
		return p*self.w*pow(self.pi, max(0, s))*pow(1-self.pi, max(0, -s))
	def __str__(self):
		return "N: " + str(self.mu) + "," + str(self.si) + "," + str(self.l) + "," + str(self.w)+ "," + str(self.pi)


class UNI:
	def __init__(self, a,b, w, pi):
		self.a, self.b 	= a,b
		self.w, self.pi = w,pi
		self.type 		= "U"
	def pdf(self, x, st):
		
		return int(self.a<=x<=self.b)*(self.w/(self.b-self.a))*pow(self.pi, max(0, st))*pow(1-self.pi, max(0, -st))
	def __str__(self):
		return "U: " + str(self.a) + "," + str(self.b) + "," + str(self.w) + "," + str(self.pi)
class model:
	def __init__(self, ll , converged, retries):
		self.ll 		= ll
		self.converged 	= converged
		self.retries 	= retries
		self.rvs 		= list()
	def pdf(self, x, st):
		return sum([rv.pdf(x, st) for rv in self.rvs])
	def __str__(self):
		return "\n".join([rv.__str__() for rv in self.rvs])

class info:
	def __init__(self, chrom, start, stop, N):
		self.chrom 	= chrom
		self.start 	= start
		self.stop 	= stop
		self.models = {}
		self.current= None
		self.N 		= N
		self.forward	= list()
		self.reverse 	= list()
		self.X 		= None
	def add_model(self, line):
		k,ll, converged, retries 	= line[1:].strip("\n").split(",")
		k,ll, retries 				= int(k), float(ll), float(retries)
		converged 					= bool(converged=="True")
		assert k not in self.models, "there should only be one model entry"
		self.models[k] 	= model(ll , converged, retries)
		self.current 	= k
	def add_components(self, line):
		TYPE, Info 	= line.strip("\n").split(": ")
		if TYPE=="N":
			mu, si, l, w, pi 	= Info.strip("\n").split(",")
			mu, si, l, w, pi 	= float(mu), float(si), float(l),float(w), float(pi)
			self.models[self.current].rvs.append(EMG(mu, si, l, w, pi))
		elif TYPE=="U":
			a,b,w, pi 	= Info.strip("\n").split(",")
			a,b,w, pi 	= float(a), float(b), float(w), float(pi)
			self.models[self.current].rvs.append(UNI(a, b, w, pi))
	def add_data(self, x,y, forward=False):
		if forward:
			self.forward.append((x,y))
		else:
			self.reverse.append((x,y)) 
	def BIN(self, bins=100):
		start, stop 	= min((min(self.forward)[0],min(self.reverse)[0])),max((max(self.forward)[0],max(self.reverse)[0]))
		self.X 			= np.zeros((bins, 3))
		self.X[:,0] 	= np.linspace(start, stop, bins)
		for j,f in enumerate((self.forward, self.reverse)):
			for x,c in f:
				i  	= 0
				while i <self.X.shape[0] and self.X[i,0] <= x:
					i+=1
				self.X[i-1, j+1]+=c
		
def EMG_out(FILE):
	I 		= None
	fits 	= list() 	
	with open(FILE) as FH:
		for line in FH:
			if "#"==line[0]:
				chrom, (stsp,N)  = line[1:].strip("\n").split(":")[0],line[1:].strip("\n").split(":")[1].split(",")
				start, stop 	= stsp.split("-")
				start, stop 	= float(start), float(stop)
				if I is not None:
					fits.append(I)
				I 				= info(chrom,start, stop, float(N))
			elif "~" == line[0]:
				I.add_model(line)
			elif I is not None:
				I.add_components(line)
	if I is not None:
		fits.append(I)
	return fits


def insert_data(FILE, intervals):
	with open(FILE) as FH:
		i 	= -1
		for line in FH:
			if line[0] == "#":
				i+=1
				chrom, start, stop 	= line[1:].strip('\n').split(',')
				GO,forward 					= False, False
				if i < len(intervals):
					if (chrom==intervals[i].chrom and int(float(start))==(intervals[i].start) 
						and int(float(stop))==(intervals[i].stop)):
						GO 	= True

			elif "~" == line[0] and GO:
				if "forward" in line:
					forward 	= True
				else:
					forward 	= False
			elif GO :
				x,y 	= line.strip("\n").split(",")
				x,y 	= float(x),float(y)
				intervals[i].add_data(x,y, forward=forward)

		










