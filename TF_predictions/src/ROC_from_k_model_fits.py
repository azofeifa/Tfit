import time
import matplotlib.pyplot as plt
import numpy as np
import math
class segment:
	def __init__(self,line):
		self.z 		= None
		self.DNase 	= False
		self.nearest_gene 	= None
		self.chrom,self.start,self.stop, self.fN,self.rN, self.N 	= [None for i in range(6)]
		self.K 	= {}
		self._parse_header(line)
	def _parse_header(self, line):
		#>CM_4,1|chrX:17429566-17433572|1260.000000,208.000000
		line_array 	= line[1:].strip("\n").split("|")
		self.chrom,stsp 		= line_array[1].split(":")
		self.start, self.stop 	= [int(x) for x in stsp.split("-")]
		self.fN, self.rN 		= [float(x) for x in line_array[-1].split(",")]
		self.N 					= self.fN+self.rN
		if "MO" 	== line[1:3]:
			self.z 	= 0
		else:
			self.z 	= 1
	def add_fit(self, line):
		# 		~0,-6016.446363
		# ~1,-4285.457147 17430712.388941 29.646927       20.000000       0.609224        292.240000      0.512216,0.409243,0.076397      17433572.000000 17429566.000000
		line_array 	= line[1:].strip("\n").split("\t")
		if len(line_array)<3:
			self.K[0] 	= float(line_array[0].split(",")[-1])
		else:
			self.K[1] 	= float(line_array[0].split(",")[-1])
			self.mu,self.si, self.l, self.fp 	= [float(x.split(",")[0]) for x in line_array[1:5]]
			self.wp 							= float(line_array[6].split(",")[0])
	def get_BIC_mdoel(self, penality):
		if not self.N:
			return 0;
		NULL 	= -2*self.K[0] + math.log(self.N)
		MODEL 	= -2*self.K[1] + penality*5*math.log(self.N)
		if NULL < MODEL:
			return 0
		return 1
	def check(self):
		if self.N < 10:
			return False
		if self.N > 100 and not self.z:
			return False
		if self.K[1]==-np.inf or self.K[0]==-np.inf:
			return False
		if self.K[1] < self.K[0]:
			return False
		return True

def load_k_model_fits(FILE):
	S 	= None
	G 	= list()
	with open(FILE) as FH:
		for line in FH:
			if line[0]==">":
				if S is not None and S.check():
					G.append(S)
				S 	= segment(line)
			elif S is not None:
				S.add_fit(line)
	if S is not None and S.check():
		G.append(S)
	return G

def compute_draw_ROC(ax, G,label):
	penalities 	= np.linspace(-3,8., 100) 
	TPN, TNN 	= float(len([1 for s in G if s.z ])), float(len([1 for s in G if not s.z ]))
	TPS,TNS 	= list(),list()
	print TNN, TPN, label
	vert 		= ()
	for p in penalities:
		TN,TP 	= 0.0,0.0
		for s in G:
			z 	= s.get_BIC_mdoel(pow(10,p) )
			if z == s.z and not s.z:
				TN+=1
			if z == s.z and s.z:
				TP+=1
		TNS.append(1-(TN / TNN))
		TPS.append(TP / TPN)
	ax.plot( TNS, TPS, label=label+", AUC: " + str(sum([ (TNS[i-1]-TNS[i])*TPS[i]  for i in range(1, len(TPS)  )]))[:5]  ,linewidth=2)

def ROC(GS, labels):
	F 	= plt.figure()
	ax	= F.add_subplot(1,1,1)
	for i,G in enumerate(GS):
		compute_draw_ROC(ax, G, labels[i])
	ax.plot([0,1], [0,1], linewidth=2, linestyle="--")
	ax.legend(loc=(0.7,0.1))
	
	ax.set_xlim(-0.1,1.1)
	ax.set_ylim(0,1)
	ax.grid()
	ax.set_xlabel("False Positive Rate")
	ax.set_ylabel("True Positive Rate")
	
	plt.show()



if __name__ == "__main__":
	DIR 		= "/Users/joazofeifa/Downloads/"
	G_JunD 		= load_k_model_fits(DIR + "JunD_motif-1_K_models_MLE.tsv")
	G_ATF 		= load_k_model_fits(DIR+"ATF_motif-1_K_models_MLE.tsv")
	G_CBX3 		= load_k_model_fits(DIR+"CBX3_motif-1_K_models_MLE.tsv")
	G_ELF1 		= load_k_model_fits(DIR + "ELF1_motif-1_K_models_MLE.tsv")
	G_SRF 		= load_k_model_fits(DIR + "SRF_motif-1_K_models_MLE.tsv")
	G_USF1 		= load_k_model_fits(DIR +"USF1_motif-1_K_models_MLE.tsv")

	ROC((G_ATF,G_CBX3,G_ELF1 ,G_SRF,G_JunD,G_USF1), ("ATF", "CBX3", "ELF1", "SRF", "JunD", "USF1" ))






