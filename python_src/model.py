import numpy as np
from scipy.special import erf, erfc
from scipy.stats import invgamma
import simulate
import math as m
import math
import matplotlib.pyplot as plt
import time
import multiprocessing as mp
import matplotlib as mpl
import matplotlib.cm as cm
import load 
import template_window_matching as twm

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
		self.pi, self.w =(self.r[1] + self.c.beta_0) / (r+ self.c.beta_0*2), (r + self.c.alpha_0)  / (N+self.c.alpha_0*self.c.K*3 + self.c.K*3 )
		if self.type=="noise":
			self.w 		= min(self.c.noise_max, self.w)
			print "HERE", self.w
		#====================================================
		self.r[1], self.r[-1] 			= 0,0
		self.ri[1], self.ri[-1] 		= 0,0

		if self.type=="forward" :
			self.a 			= self.bidir.mu +self.foot_print
			if abs(self.b - self.a) < 2:
				self.b 		= self.a + min(np.random.gamma(self.c.uniform_rate, 1), self.c.maxX)
		elif self.type=="reverse" :
			self.b 			= self.bidir.mu - self.foot_print
			if abs(self.b - self.a) < 2:
				self.a 			= self.b - min(np.random.gamma(self.c.uniform_rate, 1), self.c.maxX)

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


		self.type 					= "EMGU"
		self.c 						= classifier #larger classifier wrapper, which has the hyperparameters
		self.remove 				= False
		self.forward  				= None #elongation component to the right
		self.reverse 		 		= None #elongation component to the left
		self.foot_print 			= foot_print


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
	def check(self):
		if self.remove:
			self.forward.remove, self.reverse.remove 	= True, True
			return True
		if math.isnan(self.mu) or math.isnan(self.l) or math.isnan(self.si):
			self.forward.remove, self.reverse.remove 	= True, True
			return True
		return False 

	
	def add_stats(self, z, forward_ct, reverse_ct, norm_forward, norm_reverse):
				
		if forward_ct and norm_forward:
			r 	= forward_ct*self.ri[1]/norm_forward
			E_Y = self.EY(z,1)
			E_X = z-E_Y -self.foot_print 

			self.EX_f+=E_X *r
			self.EY_f+=E_Y*r
			self.E_X2+=(pow(E_X, 2) + self.EY2(z,1) - pow(E_Y, 2))*r
			self.r[1]+=r
			self.C+=max( ((z-self.mu) -E_Y) *r,0)
		if reverse_ct and norm_reverse:
			r 	= reverse_ct*(self.ri[-1]/norm_reverse)
			E_Y = self.EY(z,-1)
			E_X = z+E_Y + self.foot_print

			self.EX_r+=E_X*r
			self.EY_r+=E_Y*r
			self.E_X2+=(pow(E_X, 2) + self.EY2(z,-1) - pow(E_Y, 2))*r
			self.r[-1]+=r

			self.C+=max((-(z-self.mu) -E_Y)   *r ,0)
			

	def eval(self, z,forward_ct, reverse_ct):
		if forward_ct:
			self.ri[1] 	= self.pdf(z, 1)
		else:
			self.ri[1] 	= 0.
		if reverse_ct:
			self.ri[-1] = self.pdf(z,-1)
		else:
			self.ri[-1] = 0.
		return self.ri[1], self.ri[-1]
	def get_current(self):
		return self.ri[1],self.ri[-1]



		
	def set_new_parameters(self, N): #M-step
		print self.c.alpha_1

		r  								= self.r[1] + self.r[-1]
		self.pi, self.w 				=(self.r[1] + self.c.beta_0) / (r+ self.c.beta_0*2), (r + self.c.alpha_0)  / (N+self.c.alpha_0*self.c.K*3 + self.c.K*3 )
		self.mu 						= (self.EX_f+self.EX_r)  / (r+ 0.1)
		self.si 						= pow(abs((1./ (r + 3 + self.c.alpha_1)  )*(self.E_X2 - 2*self.mu*(self.EX_f+self.EX_r) + r*pow(self.mu, 2) + 2*self.c.beta_1 + self.c.tau*pow(self.mu-self.c.m_0, 2)   )),0.5)
		self.l 							= 1.0 /(((self.EY_f+self.EY_r) + self.c.beta_2) / (r + self.c.alpha_2))
		self.l 							= min(2,self.l)
		self.foot_print 				= min((self.C / (r+0.1)), 20)

		#====================================================
		self.r[1], self.r[-1] 			= 0,0
		self.ri[1], self.ri[-1] 		= 0,0
		#====================================================
		self.EX_f,self.EX_r, self.EY_f,self.EY_r, self.E_X2,self.C 	= 0.,0.,0.,0.0,0.0,0.0
		
	def reset(self):
		#=======================================
		#new parameters for the Bidirectional
		self.mu 		= np.random.uniform(self.c.minX, self.c.maxX)
		self.l 			= 1.0/np.random.gamma((self.c.maxX-self.c.minX)/(10*self.c.K), 1)
		self.si 		= np.random.gamma((self.c.maxX-self.c.minX)/(10*self.c.K), 1)
		self.pi,self.w 	= 0.5, 1.0 / (self.c.K*3) 
		#=======================================
		#new parameters for the forward uniform
		self.forward.a 	= self.mu
		self.forward.b 	= self.mu + min(np.random.gamma(self.c.uniform_rate, 1), self.c.maxX)
		self.w 			= 1.0 / self.c.K
		#=======================================
		#new parameters for the reverse uniform
		self.reverse.b 	= self.mu
		self.reverse.a 	= self.mu - max(np.random.gamma(self.c.uniform_rate, 1), self.c.minX)
		self.reverse.w 	= 1.0 / self.c.K
		
		

		self.remove 	= False
		



class EMGU:
	def __init__(self, max_ct=pow(10, -4), max_it=200, K=2, bayes=False, noise=True, 
		noise_max=0.1, moveUniformSupport=5, cores=1,
		seed=False, foot_print=0):
		self.max_ct 	= max_ct
		self.max_it 	= max_it
		self.K 			= K
		self.bayes 		= bayes
		self.noise 		= noise
		self.noise_max 	= noise_max
		self.move 		= moveUniformSupport
		self.cores 		= cores
		self.converged 	= False
		self.resets 	= 0.
		self.foot_print = foot_print
		#=================================
		#prior parameters
		self.alpha_0 				= 1. #symmetric prior for mixing weights
		self.beta_0 				= 1. #symmetric prior for strand probabilities
		self.m_0, self.tau 			= 0, 1 #priors for component mus
		self.alpha_1, self.beta_1 	= 1, 1 #priors for component sigmas
		self.alpha_2, self.beta_2 	= 1, 1 #priors for component 
		self.peaks 					= None #prior on where the bidirectionals are from our template/bayes factor analysis
		self.seed 					= seed
		self.uniform_rate 			= None
		#=================================
		self.rvs,self.ll 			= None, None #final group of components
	def pdf(self, x,s):
		assert self.rvs is not None, "need to running fit before evalauting mixture pdf"
		vl 	= sum([rv.pdf(x, s) for rv in self.rvs ])
		return vl


	def LOG(self, x):
		if x == 0:
			return -np.inf
		return m.log(x)



	def fit(self, X):



		#=======================================
		#randomally initialize parameters
		minX, maxX 	= X[0,0], X[-1,0]
		self.minX, self.maxX 	= minX, maxX


		ws 			= np.random.dirichlet([self.alpha_0]*self.K*3).reshape(self.K, 3)
		pis 		= np.random.beta(self.beta_0, self.beta_0, self.K*3).reshape(self.K,3)
		if self.peaks is None:
			mus 		=  np.random.uniform(minX, maxX, self.K)
		else:
			mus 		= [x for x  in self.peaks]
		mus 		= [30]
		sigmas 		= np.random.gamma((maxX-minX)/(35*self.K), 1, self.K)
		lambdas 	= 1.0/np.random.gamma((maxX-minX)/(25*self.K), 1, self.K)
		#=======================================
		#assign to components

		self.uniform_rate= (maxX-minX)/(1*self.K)
		fp 			= np.random.uniform(0,2)
		print "*********", fp
		bidirs 		= [component_bidir(mus[k], sigmas[k], lambdas[k], ws[k][0], pis[k][0],self, foot_print=fp) for k in range(self.K)] 
	
		uniforms    = [component_elongation(minX, mus[k], ws[k][1], 0.5, bidirs[k], "reverse",self ,0 , foot_print=fp) for k in range(self.K)]
		uniforms   += [component_elongation(mus[k],maxX, ws[k][2], 0.5, bidirs[k], "forward",self, X.shape[0] , foot_print=fp) for k in range(self.K)]
		
		if self.noise:
			print "HERE?"
			uniforms+=[component_elongation(minX, maxX, 0.1, 0.5, bidirs[0], "noise", self, X.shape[0]) ]
		components 			= bidirs + uniforms
		N_f, N_r 			= sum(X[:,1]), sum(X[:,2])
		t,converged 		= 0 , False
		ll, prevll 			= 0., -np.inf
		st 					= time.clock()
		iter_parameters 	= list()
		
		while t < self.max_it and not converged:
			self.rvs 		= [c for c in components ]
			self.draw(X)
			
			# for rv in self.rvs:
			# 	if rv.type=="EMGU":
			# 		print rv
			
			#######
			#E-step
			#####
			ll 	= 0.
			st 	= time.clock()
			for i in range(X.shape[0]):
				#compute component responsibilities
				norm_forward, norm_reverse 	= 0.,0.
				for c in components:
					forward_eval, reverse_eval 	 = c.eval(X[i,0],X[i,1],X[i,2])
					norm_forward 				+=forward_eval
					norm_reverse 				+=reverse_eval
				#add sufficient stats
				for c in components: #add sufficient stats
					c.add_stats(X[i,0],X[i,1],X[i,2], norm_forward, norm_reverse)
				if norm_forward:
					ll+=math.log(norm_forward)*X[i,1]
				if norm_reverse:
					ll+=math.log(norm_reverse)*X[i,2]
			print ll
			#######
			#M-step
			#######
			
			N 	= sum([sum(c.r.values()) for c in components])
			for c in components:
				c.set_new_parameters(N)
				
			#check which components blew up and reset
			i  	= 0
			
			st 	= time.clock()
			#check for convergence
			if abs(ll-prevll) < self.max_ct:
				converged=True
				self.converged=True
			#move LLs...unfortunately this is brute - force...
			
			prevll 	= ll
			# for c in components:
			# 	print c
			# print "------------"
			t+=1
		self.rvs,self.ll 	= [c for c in components if c.type!="noise"], ll
	def draw(self, X):
		assert self.rvs is not None, "need to run fit before drawing"
		F 			= plt.figure(figsize=(15,10))
		for c in self.rvs:
			print c
		ax 			= F.add_subplot(111)
		ax.bar(X[:,0],  X[:,1] / float(np.sum(X[:,1:]) ), color="blue", alpha=0.25, width=(X[-1,0]-X[0,0])/X.shape[0])
		ax.bar(X[:,0], -X[:,2] / float(np.sum(X[:,1:])), color="red", alpha=0.25, width=(X[-1,0]-X[0,0])/X.shape[0])
		xs 			= np.linspace(X[0,0], X[-1,0], 1000)
		ys_forward 	= map(lambda x: sum([rv.pdf(x, 1) for rv in self.rvs if rv.type == "EMGU"]) , xs) 
		ys_reverse 	= map(lambda x: -sum([rv.pdf(x, -1) for rv in self.rvs if rv.type == "EMGU"]) , xs)
		
		ax.plot(xs, ys_forward, linewidth=3.5,  color="blue")
		ax.plot(xs, ys_reverse, linewidth=3.5,  color="red")
		ax.grid()
		plt.show()
	


if __name__ == "__main__":
	#==================================
	#testing MAP-EM procedure

	#X 	= simulate.runOne(mu=0, s=1, l=10, lr=100, ll=-50, we=0.5,wl=0.25, wr=0.25, pie=0.5, pil=0.1, pir=0.9, 
#		N=1000, SHOW=False, bins=300, noise=False, foot_print=10 )
	# chr1:20,984,647-20,991,448
	#chr1:836,835-843,549
 	X 		=  load.grab_specific_region("chr1",836835, 843549, 
			pos_file="/Users/joazofeifa//Lab/gro_seq_files/HCT116/bed_graph_files/DMSO2_3.pos.BedGraph", 
			neg_file="/Users/joazofeifa//Lab/gro_seq_files/HCT116/bed_graph_files/DMSO2_3.neg.BedGraph",
			SHOW 	=False, bins=300)
 	X[:,0]-=X[0,0]
 	X[:,0]/=100.
	clf = EMGU(noise=True, K=1,noise_max=0.1,moveUniformSupport=0,max_it=200, cores=1, 
		seed=True )
	clf.fit(X)


	# # #make test_file
	# FHW_f 	= open("three_prime_forward.bedgraph", "w")
	# FHW_r 	= open("three_prime_reverse.bedgraph", "w")
	# for i in range(X.shape[0]):
	# 	x 		= int(X[i,0])

	# 	FHW_f.write("chr1\t" + str(x)  + "\t" + str(x) + "\t" + str(int(X[i,1])) + "\n" )
	# 	FHW_r.write("chr1\t" + str(x)  + "\t" + str(x) + "\t" + str(int(X[i,2])) + "\n" )






	#2,518,131-2,523,183
	#chr1:25,656,111-25,693,262
	#chr1:28,569,727-28,580,179
	#chr1:1,403,824-1,426,126
	#chr1:1,644,778-1,661,469
	#chr1:3,757,936-3,793,447
	
	# print max(X[:,0])



	# clf.draw(X)
	#==================================
	









