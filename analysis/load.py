import merge_data_types as mdt
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
sys.path.append( "/".join(os.getcwd().split("/")[:-1]) + "/python_src/")
sys.path.append( os.getcwd() + "/python_src/")
sys.path.append("/Users/azofeifa/Lab/EMG/python_src/")
import model
import simulate
import math 
def merge_data_out(FILE,just_params=False):
	G=list()
	with open(FILE) as FH:
		for line in FH:
			if "#" == line[0]:
				chrom,info 			= line[1:].strip("\n").split(":")
				start_stop, N,aN 	= info.split(",")
				start,stop 			= start_stop.split("-")
				G.append(mdt.segment(chrom,int(start),int(stop),float(N), annotation_N=int(aN)))
			elif "~" == line[0]:
				G[-1].insert_model_info(line)
			elif "N:"==line[:2] or "U:"==line[:2]:
				G[-1].insert_component(line)
			elif not just_params:
				line_array 				= line.strip("\n").split(",")
				data_type,peak, data 	= line_array[0], line_array[1],",".join(line_array[2:])
				if data_type != "dbSNP":
					data 					= [(float(d.split(",")[0]),float(d.split(",")[1])) for d in data.split(":") ]
				else:
					data 					= [(float(d.split(",")[0]), d.split(",")) for d in data.split(":")  ]
				setattr(G[-1], data_type, data)
				setattr(G[-1], data_type+"_peak", bool(peak=="True"))
				if not hasattr(G[-1], "data_types"):
					setattr(G[-1], "data_types", list())
				G[-1].data_types.append(data_type)
	return G
import time
class model:
	def __init__(self, mu, si, lam, pi, wEM, wF, wR, bF, bR):
		self.mu 	= mu
		self.si 	= si
		self.lam 	= lam
		self.pi 	= pi
		self.wEM 	= wEM
		self.wF 	= wF
		self.wR 	= wR
		self.N 		= 0
		self.bF 	= bF
		self.bR 	= bR
		self.OKAY 	= True
		self.std 	= self.si + self.lam
		self.start 	= self.mu-self.std
		self.stop 	= self.mu+self.std
		self.chrom 	= None
		self.check()
		self.annotated = False
		self.p53_site 	 = False
		self.dist 	= None
		self.fp 	= 0
		self.overlap 	= False
	def check(self):
		parameters 	= ("mu","si","lam","pi","wEM","wF", "wR" ,"bF","bR" )
		for p in parameters:
			if math.isnan(getattr(self, p)):
				self.OKAY=False
		if self.lam == 25:
			self.OKAY=False
class segment:
	def __init__(self, line, param_txt = False):
		self.models 		= list()
		self.annotated 		= False
		self.p53_site 	 	= False
		self.dist 			= None
		self.N 				= 0

		if param_txt:
			line_array 			= line.strip("\n").split("\t")
			self.chrom,st_sp 	= line_array[0].split(":")
			self.start, self.stop = [int(x) for x in st_sp.split("-")]
			self.ID 			= line_array[1]
			self.ll, self.N 	= [float(x) for x in line_array[2].split(",")]
			self.models 		= list()
			self.OKAY 			= True
			self.add_models(line_array)
			self.K 				= len(self.models)

		else:
			chrom,start, stop, parameters 	= line.strip("\n").split("\t")
			self.chrom 				= chrom 
			self.start, self.stop 	= int(start), int(stop)
			mu,si,l,w, pi, N, fp,ll 	= parameters.split("_")
			self.models.append(model(float(mu),float(si),float(l),float(pi),float(w) , 0,0,0,0   )  )
			self.models[-1].N 		= float(N)
			self.models[-1].fp 		= float(fp)
			self.models[-1].chrom 	= self.chrom
			self.models[-1].density = float(N) / (self.stop - self.start)
		self.K 				= len(self.models)
			
	def add_models(self, line_array):
		mus 	= [float(m) for m in line_array[3].split(",")]
		sis 	= [float(m) for m in line_array[4].split(",")]
		lams 	= [float(m) for m in line_array[5].split(",")]
		pis 	= [float(m) for m in line_array[6].split(",")]
		wsEMG 	= [float(m) for m in line_array[7].split(",")]
		wsF 	= [float(m) for m in line_array[8].split(",")]
		wsR 	= [float(m) for m in line_array[9].split(",")]
		bF 		= [float(m) for m in line_array[10].split(",")]
		bR 		= [float(m) for m in line_array[11].split(",")]
		for i, mu in enumerate(mus):
			self.models.append(model(mu, sis[i], lams[i], pis[i], wsEMG[i], wsF[i], wsR[i], bF[i], bR[i] ))
		self.OKAY 	= bool(sum([1 for m in self.models if not m.OKAY  ]) == 0)
def load_model_fits_txt_file(FILE):
	G 	= {} #sorted by unique ID 
	L 	= {} #by chrom and then a list of start and stops
	with open(FILE) as FH:
		for line in FH:
			if line[0]!= "#":
				S 	= segment(line,param_txt = True)
				if S.OKAY:
					G[S.ID] 	= S
					if S.chrom not in L:
						L[S.chrom] 	= list()
					L[S.chrom].append((S.start, S.stop, S))			
	for chrom in L:
		L[chrom].sort()
	return L,G

def load_model_fits_bed_file(FILE):
	G 	= {} #sorted by unique ID 
	L 	= {} #by chrom and then a list of start and stops
	with open(FILE) as FH:
		for line in FH:
			if line[0]!="#":
				S 	= segment(line, param_txt = False)

				if S.chrom not in L:
					L[S.chrom] 	= list()
				if 0.99 > S.models[0].wEM >0 and S.models[0].lam!=25:
					L[S.chrom].append((S.start, S.stop, S))
	for chrom in L:
		L[chrom].sort()
	return L,G

def load_Refseq(FILE, pad=1000):
	#just promoter
	G 	={}
	header 	= True
	with open(FILE) as FH:
		for line in FH:
			if not header:
				name, chrom, strand, start, stop 	= line.strip("\n").split("\t")[1:6]
				if strand 	== "-":
					start,stop 	= stop, start
				start, stop 	= int(start), int(stop)
				if chrom not in G:
					G[chrom] 	= list()
				G[chrom].append([start-pad, start+pad])
			else:
				header=False
	#now merge
	A={}
	for chrom in G:
		N = len(G[chrom])
		L = G[chrom]
		L.sort()
		A[chrom] = list()
		j=0
		while j < N:
			curr = L[j]
			while j < N and L[j][1]>curr[0] and L[j][0] < curr[1]:
				curr[0] = min(curr[0], L[j][0])
				curr[1] = max(curr[1], L[j][1])
				j+=1
			A[chrom].append(curr)
	return A

def load_ChIP_p53(FILE):
	G 	= {}
	with open(FILE) as FH:
		for line in FH:
			chrom,start, stop 	= line.strip("\n").split("\t")[:3]
			if chrom not in G:
				G[chrom] 		= list()
			start, stop 		= int(start),int(stop)
			G[chrom].append((start, stop))
	for chrom in G:
		G[chrom].sort()
	return G

def label(A,B, attr=""):
	G 	={}
	F 	= {}
	for chrom in A:
		if chrom in B:
			a,b 	= A[chrom], B[chrom]
			a.sort()
			b.sort()
			j,N 	= 0,len(b)
			for a_st, a_sp, a_M in a:
				while j < N and b[j][1] < a_st:
					j+=1
				if j < N and b[j][0] < a_sp:
					setattr(a_M, attr, True)
					if attr == "annotated":
						center 	= (b[j][0] + b[j][1])/2.					
						setattr(a_M, "dist", a_M.models[0].mu - center)
						pass
					string = chrom + ":" +str(b[j][0])+ "-"+ str(b[j][1])
					F[string] 	= 1
				elif j < N:
					string = chrom + ":" +str(b[j][0])+ "-"+ str(b[j][1])
					G[string] 	= 1















if __name__ == "__main__":
	pass	
