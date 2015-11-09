import numpy as np
import matplotlib.pyplot as plt
import math
import time
class segment:
	def __init__(self,line):
		self.line 	= line
		self.type 	= int(line[1:].split("|")[0]!="NOISE")
		self.N 		= sum([float(x) for x in line[1:].split("|")[2].split(",")])
		self.K 		= {}
		self.w 	= None
		self.remove = False
	def add_model(self, line):
		k,ll 		= line[1:].split("\t")[0].split(",")
		k,ll 		= float(k), float(ll)
		if k == 1:
			self.w 	= float(line.split("\t")[6].split(",")[0])
		
		self.K[k] 	= ll
	def check(self):
		if self.K[0] > self.K[1]:
			return False
		if self.K[1] == -np.inf or self.K[0] == -np.inf:
			return False
		if self.type == 0 and self.w < 0.3:
			return True
		if self.type == 1 and self.w > 0.6:
			return True
		
		return False
	def get_BIC_model(self, penality):
		noise 	= -2*self.K[0] + math.log(self.N)
		alt 	= -2*self.K[1] + penality*7*math.log(self.N)
		if alt < noise:
			return 1
		return 0

def load_tfit_out(FILE):
	G 	= list()
	S 	= None
	with open(FILE) as FH:
		for line in FH:
			if   ">" == line[0]:
				if S is not None  and S.check():
					G.append(S)
				S 	= segment(line)
			elif "~" == line[0]:
				S.add_model(line)
	return G

def ROC(G):

	F 			= plt.figure(figsize=(15,10))
	ax 			= F.add_subplot(1,1,1)
	penalities  = np.linspace(0,5.5, 200)
	P,N 		= float(len([1 for s in G if s.type ])),float(len([1 for s in G if not s.type ]))
	TPS, FPS 	= list(),list()
	for p in penalities:
		TN,TP 	= 0.0,0.0
		for s in G:
			z 	= s.type
			
			c 	= s.get_BIC_model(pow(10,p))
			if (z==c and not c):
				TN+=1
			elif (z==c and c):
				TP+=1
		TPS.append(TP/P)
		FPS.append(1.0-(TN/N))
	AUC 	= sum([ (FPS[i-1] - FPS[i])*TPS[i]  for i in range(1, len(FPS)) ])
	ax.plot(FPS, TPS, linewidth=3., label="AUC: "+str(AUC)[:5])
	ax.set_xlim(0,1)
	ax.set_ylim(0,1)
	ax.grid()
	ax.legend(loc=(0.8,0.1))

	plt.show()






if __name__ == "__main__":
	FILE 	="/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/ISO_STAT-1_K_models_MLE.tsv"
	G 	 	= load_tfit_out(FILE)
	ROC(G)
	pass