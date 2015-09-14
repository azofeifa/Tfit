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
		self.bF 	= bF
		self.bR 	= bR
		self.OKAY 	= True
		self.std 	= self.si + self.lam
		self.check()
	def check(self):
		parameters 	= ("mu","si","lam","pi","wEM","wF", "wR" ,"bF","bR" )
		for p in parameters:
			if math.isnan(getattr(self, p)):
				self.OKAY=False
		if self.lam == 25:
			self.OKAY=False
class segment:
	def __init__(self, line):
		line_array 			= line.strip("\n").split("\t")
		self.chrom,st_sp 	= line_array[0].split(":")
		self.start, self.stop = [int(x) for x in st_sp.split("-")]
		self.ID 			= line_array[1]
		self.ll, self.N 	= [float(x) for x in line_array[2].split(",")]
		self.models 		= list()
		self.OKAY 			= True
		self.add_models(line_array)
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
				S 	= segment(line)
				if S.OKAY:
					G[S.ID] 	= S
					if S.chrom not in L:
						L[S.chrom] 	= list()

					L[S.chrom].append((S.start, S.stop, S))			
	for chrom in L:
		L[chrom].sort()
	return L,G




if __name__ == "__main__":
	pass	
