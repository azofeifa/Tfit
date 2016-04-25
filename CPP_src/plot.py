import sys
import matplotlib.pyplot as plt
import numpy as np
import math as m
from scipy.special import erf, erfc
from scipy.stats import invgamma
X=list()
collect=False
class component_elongation:
	def __init__(self, a,b, w, pi, bidir_component, ty, classifier, j, foot_print=0):
		self.a, self.b 	= a,b
		self.w, self.pi = w,pi
		self.r 			= {1:0., -1:0.}
		self.ri 		= {1:0., -1:0.}
		self.type 	 	= ty
		self.bidir 		= bidir_component #bidirectional component to which this belongs
		if self.bidir is not None:
			setattr(self.bidir, self.type, self)
		self.c 			= classifier #larger classifier wrapper, which has the hyperparameters
		self.j 			= j #index of data from which b came from
		self.remove 	= False
		self.foot_print = foot_print

	def __str__(self):
		return "U: " + str(self.a) + "," + str(self.b) + "," + str(self.w) + "," + str(self.pi)


	def eval(self, z, forward_ct, reverse_ct):
		if forward_ct :
			self.ri[1] 	= self.pdf(z, 1)
		else:
			self.ri[1] 	= 0
		if reverse_ct :
			self.ri[-1] = self.pdf(z,-1)
		else:
			self.ri[-1] = 0.
		self.r[1]+=self.ri[1]
		self.r[-1]+=self.ri[-1]

		return self.ri[1], self.ri[-1]
	def get_current_r(self):
		return self.ri[1] + self.ri[-1]
	def set_new_parameters(self, N):
		r 				= (self.r[1] + self.r[-1])
		self.w = (r + self.c.alpha_0)  / (N+self.c.alpha_0*self.c.K*3 + self.c.K*3 )
		if self.type=="noise":
			self.w 		= min(self.c.noise_max, self.w)
		#====================================================
		self.r[1], self.r[-1] 			= 0,0
		self.ri[1], self.ri[-1] 		= 0,0

		if self.type=="forward" :
			self.a 			= self.bidir.mu +self.foot_print
			if abs(self.b - self.a) < 2:
				self.b 		= self.a + min(np.random.gamma(self.c.uniform_rate, 1), self.c.maxX)
			self.pi 		= 1.0
		elif self.type=="reverse" :
			self.b 			= self.bidir.mu - self.foot_print
			if abs(self.b - self.a) < 2:
				self.a 			= self.b - min(np.random.gamma(self.c.uniform_rate, 1), self.c.maxX)
			self.reverse 	= 0.0

	def check(self):
		if math.isnan(self.w) or math.isnan(self.a) or math.isnan(self.b):
			return True
		return False
	def reset(self):
		self.remove 		= False


	def add_stats(self,z, forward_ct, reverse_ct, norm_forward, norm_reverse):
		if forward_ct and norm_forward and (self.type=="forward" or self.type=="noise"):
			r 	= forward_ct*self.ri[1]/norm_forward
			self.r[1]+=r
		if reverse_ct and norm_reverse and (self.type=="reverse" or self.type=="noise"):
			r 	= reverse_ct*(self.ri[-1]/norm_reverse)
			self.r[-1]+=r

	def pdf(self, z,s,move_a=0, move_b=0):
		if s== 1:
			PI 	= self.pi
		else:
			PI 	= (1-self.pi)


		vl 	= int(self.a+move_a<=z<=self.b+move_b)*(self.w / abs((self.b+move_b) - (self.a + move_a) ))*PI


		return vl

class component_bidir:
	def __init__(self, mu, si, l, w,pi , classifier,foot_print=0):
		#============================
		self.mu, self.si, self.l 	= mu, si, l
		self.w, self.pi 			= w,  pi
		#============================
		self.r 						= {1:0., -1:0.} #running total
		self.ri 					= {1:0., -1:0.} #current
		#============================
		self.EX_f,self.EX_r 		= 0.,0. #running total
		self.EY_f,self.EY_r 		= 0.,0. #running total
		self.E_X2 					= 0. #running total
		self.C 						= 0.
		self.Cf 					= 0.
		self.Cr 					= 0.


		self.type 					= "EMGU"
		self.c 						= classifier #larger classifier wrapper, which has the hyperparameters
		self.remove 				= False
		self.forward  				= None #elongation component to the right
		self.reverse 		 		= None #elongation component to the left
		self.foot_print 			= foot_print
		self.prev_mu 				= self.mu
		self.move_fp 				= False
	def __str__(self):
		return "N: " + str(self.mu) + "," +str(self.si) + "," +str(self.l) + "," + str(self.w) + "," + str(self.pi)
		return ("N: " + str(self.mu) + "," +str(self.si) + "," +str(self.l) + "," + str(self.w) + "," + str(self.pi) + "\n" +
			self.forward.__str__() + "\n" +
			self.reverse.__str__())

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

	def pdf(self, z,s,move_a=0, move_b=0):
		if s == 1:
			z-=self.foot_print
		else:
			z+=self.foot_print
		vl 		= (self.l /2.)* (s*2*(self.mu-z) + self.l*self.si**2. )
		if vl > 200: #overflow error, try this one...
			p 	= self.l*self.IN((z-self.mu)/self.si)*self.R(self.l*self.si - s*((z - self.mu)/self.si))
		else:
			p 		= (self.l / 2.)*m.exp( vl )*erfc( (s*(self.mu- z) + self.l*self.si**2. ) / (m.sqrt(2)*self.si) )
		step 	= 0

		return p*self.w*pow(self.pi, max(0, s) )*pow(1-self.pi, max(0, -s) )


for line in open("/Users/joazofeifa/Lab/EMG/CPP_src/test"):
    if line[:4]=="####":
        collect=True
    elif collect and line[:2]!="lo":
        x,y,z=line.strip("\n").split(",")
        X.append((float(x), float(y),float(z)))
    elif line[:2]=="lo":
        collect=False
parameters = 20.340554,1.106444,0.266913,0.522734,1.033316,0.487863,0.253398,0.256040
rv_b 		= component_bidir(parameters[0],parameters[1],parameters[2], parameters[-3], parameters[3], None, foot_print=parameters[4])
counts,edges 	= np.histogram([x for x,y,z in X], weights=[y for x,y,z in X],normed=1, bins=100)
counts2,edges2 	= np.histogram([x for x,y,z in X], weights=[z for x,y,z in X],normed=1, bins=100)
counts*=0.5
counts2*=0.5

rv_e 		= component_elongation(parameters[0]+parameters[1], edges[-1], parameters[-2],1.0, rv_b, "forward", None, 0, foot_print=0)
rv_r 		= component_elongation(edges[0], parameters[0]-parameters[1], parameters[-1], 0.0, rv_b, "reverse", None, 0, foot_print=0)
rv_n 		= component_elongation(edges[0], edges[-1], 0.1, 0.5, rv_b, "noise", None, 0.5, foot_print=0)

xs 				= np.linspace(min(edges),max(edges),1000)
rvs 			= [rv_b, rv_e, rv_r, rv_n]
print rv_b.pi
print rv_b.foot_print

ys 				=[ sum([rv.pdf(x,1) for rv in rvs ]) for x in xs]
ysr 				=[ -sum([rv.pdf(x,-1) for rv in rvs ]) for x in xs]
plt.bar(edges[1:], counts ,width=0.25, edgecolor="blue", color="blue",alpha=0.5)
plt.bar(edges[1:],-counts2, width=0.25, edgecolor="red", color="red",alpha=0.5)
plt.plot(xs, ys,lw=2)
plt.plot(xs, ysr,lw=2)
plt.show()
