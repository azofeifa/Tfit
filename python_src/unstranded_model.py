import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf, erfc
from scipy.stats import invgamma
import math as m
import time
def simulate(N=1000, mu=0,si=1,l=0.5, pi=0.5, fp=2, SHOW=True):
	D1 	= np.random.normal(mu,si, int(N*pi) ) + np.random.exponential(1.0/l, int(N*pi)) +  fp
	D2 	= np.random.normal(mu,si, int(N*(1-pi)) ) - np.random.exponential(1.0/l, int(N*(1-pi))) - fp
	D 	= [d for d in D1]+[d for d in D2]
	if SHOW:
		F 	= plt.figure()
		ax 	= F.add_subplot(1,1,1)
		ax.hist(D, bins=50)
		plt.show()
	counts,edges 	= np.histogram(D,bins=100)
	edges 			= (edges[1:] + edges[:-1]) /2.
	X 				= np.zeros((len(counts),2))
	X[:,0] 			= edges
	X[:,1] 			= counts
	return X

class NL:
	def __init__(self, mu,si, l, pi, w, fp):
		self.mu 		= mu
		self.si 		= si
		self.l 			= l
		self.pi 		= pi
		self.foot_print = fp
		self.w 	= w
	def __str__(self):
		return "NL " + str(self.mu)[:5] + "," + str(self.si)[:5] + ","+ str(self.l)[:5] + ","+ str(self.pi)[:5] 
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
	def EY(self, z, s):
		if s==1:
			z-=self.foot_print
		else:
			z+=self.foot_print

		return max((s*(z-self.mu) -self.l*pow(self.si,2) + self.si / (self.R(self.l*self.si - s*((z-self.mu)/ self.si))),0.))



	def EY2(self, z, s):
		if s==1:
			z-=self.foot_print
		else:
			z+=self.foot_print

		return pow(self.l,2)*pow(self.si,4) + pow(self.si,2)*(2*self.l*s*(self.mu-z)+1) + pow(self.mu-z,2)  - ((self.si*(self.l*pow(self.si,2) + s*(self.mu -z)) ) /  self.R(self.l*self.si - s*((z-self.mu)/ self.si)) )	


	def pdf(self, z,s,move_a=0, move_b=0):
		if s == 1:
			z-=self.foot_print
		else:
			z+=self.foot_print
		vl 		= (self.l /2.)* (s*2*(self.mu-z) + self.l*self.si**2. )
		if vl > 200: #overflow error, try this one...
			p 	= self.l*self.IN((z-self.mu)/self.si)*self.R(self.l*self.si - s*((z - self.mu)/self.si))
		try:
			p 		= (self.l / 2.)*m.exp( vl )*erfc( (s*(self.mu- z) + self.l*self.si**2. ) / (m.sqrt(2)*self.si) )
		except:
			self.remove,p 	= True,0.
		step 	= 0

		return p*self.w*pow(self.pi, max(0, s) )*pow(1-self.pi, max(0, -s) )
def estimate(X, fp=2):
	t=0
	T=500
	converged=False

	nl 	= NL(np.random.uniform(X[0,0], X[-1,0]),  2, 0.1,0.5 ,1,fp )
	N 	= float(np.sum(X[:,1]))

	#priors
	alpha_1 	= 10
	beta_1 		= 10
	alpha_3 	= 5
	while t < T and not converged:
		EX, EY, EX2, EPI 	= 0,0,0,0
		#e-step
		xs 	= np.linspace(X[0,0], X[-1,0],1000)
		plt.hist(X[:,0], weights=X[:,1], bins=50, normed=1, alpha=0.3)
		ys 	= [nl.pdf(x,1)  for x in xs]
		ysr = [nl.pdf(x,-1)  for x in xs]
		plt.plot(xs,ys, color="blue")
		plt.plot(xs,ysr, color="red")
		
		plt.show()
		ll 	= 0
		for i in range(X.shape[0]):
			x,y 	= X[i,:]
			epi 	= (nl.pdf(x,1)  ) / ((nl.pdf(x,1) )+ (nl.pdf(x,-1) ) )
			ll+=(m.log(nl.pdf(x,1)) + m.log(nl.pdf(x,-1)))*y
			ey 		= max(epi*nl.EY(x,1)+(1.-epi)*nl.EY(x,-1),0)

			ey2 	= epi*nl.EY2(x,1)+(1.-epi)*nl.EY2(x,-1)
			ex 		= (x-(ey*epi)  + (ey *(1.-epi)	))
			
			ex 		+= (nl.foot_print*(1-epi) -nl.foot_print*epi  )
			
			EPI 	+= (epi*y)
			EY 		+= (ey*y)
			EX 		+= ex*y
			EX2 	+= max((pow(ex,2) + ey2 - pow(ey,2))*y,0)
		print ll
		#m-step
		# self.pi, self.w 				=(self.r[1] + self.c.beta_0) / (r+ self.c.beta_0*2), (r + self.c.alpha_0)  / (N+self.c.alpha_0*self.c.K*3 + self.c.K*3 )
		# self.mu 						= self.E_X  / (r+ 0.1)
		# self.si 						= pow(abs((1./ (r + 3 + self.c.alpha_1)  )*(self.E_X2 - 2*self.mu*self.E_X + r*pow(self.mu, 2) + 2*self.c.beta_1 + self.c.tau*pow(self.mu-self.c.m_0, 2)   )),0.5)
		
		mu 	= EX / N
		si 	= m.sqrt((EX2 - 2*mu*EX + (pow(mu,2)*N) + 2*beta_1 ) / (N  + 3 + alpha_1)  )
		
		l 	= N /EY
		pi 	= (EPI +alpha_3)  / (N + 2*alpha_3)

		nl 	= NL(mu,si,l,pi,1, fp)
			
		t+=1
		pass



if __name__ == "__main__":
	X 				= simulate(N=1000, mu=0,si=1,l=0.5, pi=0.15,SHOW=False, fp=2                   )
	estimate(X)


