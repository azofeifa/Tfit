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
		self.w 		= None
		self.remove = False
		self.header = line[1:].split("|")[0]

	def add_model(self, line):
		line_array 	= line[1:].split("\t")
		k,ll 		= line_array[0].split(",")
		k,ll 		= float(k), float(ll)
		if k == 1:
			self.w 	= float(line_array[6].split(",")[0])
			self.si = float(line_array[2])
			self.l 	= float(line_array[3])
			self.fp = float(line_array[4])
		self.K[k] 	= ll
	def check(self):
		if self.K[0] > self.K[1] or not self.N or self.K[0]==-np.inf or self.K[1] ==-np.inf:
			return False
		if self.type and self.N < 3000:
			return False
		if self.K[1] == -np.inf or self.K[0] == -np.inf:
			return False
		if self.type == 0:
			return True
		if self.type == 1:
			return True
		
		return False
	def get_BIC_model(self, penality):
		noise 	= -2*self.K[0] + 3*math.log(self.N)
		alt 	= -2*self.K[1] + penality*7*math.log(self.N)
		if alt < noise:
			return 1
		return 0
	def get_BIC_model2(self, penality):
		MIN,ARG 	= None, None
		for k in self.K:
			BIC 	= -2*self.K[k] + (k )*penality*math.log(self.N)
			if MIN is None or BIC < MIN:
				ARG = k
				MIN = BIC
		return ARG


def load_tfit_out(FILE, CHECK=True):
	G 	= list()
	S 	= None
	with open(FILE) as FH:
		for line in FH:
			if   ">" == line[0]:
				if S is not None  and ((S.check() and CHECK)  or not CHECK ) :
					G.append(S)
				S 	= segment(line)
			elif "~" == line[0]:
				S.add_model(line)
	return G

def ROC(tple):
	F 			= plt.figure(figsize=(15,10))
	for label, G in tple:
		ax 			= F.add_subplot(1,1,1)
		penalities  = np.linspace(-1,5.5, 300)
		P,N 		= float(len([1 for s in G if s.type ])),float(len([1 for s in G if not s.type ]))
		TPS, FPS 	= list(),list()
		OPT 		= None
		prev 		= 0-7
		MAX 		= None
		ARGMAX 		= None
		for p in penalities:
			TN,TP 	= 0.0,0.0
			for s in G:
				z 	= s.type
				
				c 	= s.get_BIC_model(pow(10,p))
				if (z==c and not c):
					TN+=1
				elif (z==c and c):
					TP+=1
			
			if pow(10,p)-30 > 0 and OPT is None:
				OPT 	= 1.0-(TN/N),TP/P
			if MAX is None or math.sqrt(pow((TP/P) -(1-(TN/N)),2) +  pow((1-(TN/N)) -(TP/P),2)) > MAX:
				MAX  	= math.sqrt(pow((TP/P) -0.5,2) +  pow((1-TN/N) -0.5,2))
				ARGMAX 	= 1.0-(TN/N),TP/P

			TPS.append(TP/P)
			FPS.append(1.0-(TN/N))
			prev=pow(10,p)-7
		
		L 		= sum([ (FPS[i-1] - FPS[i] )   for i in range(1, len(FPS)) ])
		AUC 	= sum([ (FPS[i-1] - FPS[i])*TPS[i]  for i in range(1, len(FPS)) ])+(1-L)
		ax.plot(FPS, TPS, linewidth=1., label= label+ " AUC: "+str(AUC)[:5],alpha=0.8)
	ax.set_xlim(-0.1,1.1)
	ax.set_ylim(-0.1,1.1)
	ax.plot([0,1],[0,1], linewidth=3., linestyle="--" ,label="random")
	ax.set_xlabel("False Positive Rate")
	ax.set_ylabel("True Positive Rate")
	ax.grid()
	ax.legend(loc=(0.8,0.1))


	plt.show()






if __name__ == "__main__":
	FILE1 	="/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/ISO_STAT_UNI_10_500-1_K_models_MLE.tsv"
	FILE2 	="/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/ISO_STAT_UNI-1_K_models_MLE.tsv"
	FILE3 	="/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/ISO_STAT_UNI_10_5000-1_K_models_MLE.tsv"

	G1 	 	= load_tfit_out(FILE1)
	G2 	 	= load_tfit_out(FILE2)
	G3 	 	= load_tfit_out(FILE3)
	
	tple 	= [("MLE", G2),("MAP (light prior)", G1), ("MAP (heavy prior)", G3)]
	ROC(tple)
	pass