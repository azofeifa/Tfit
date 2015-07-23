import numpy as np
from scipy.special import erf, erfc
import math as m
import time
import os 
import node
import sys
import multiprocessing as mp
def load_refseq(FILE):
	G 	= {}
	with open(FILE) as FH:
		header 	= True
		for line in FH:
			if not header:
				name, chrom, strand, start, stop 	= line.split("\t")[1:6]
				if chrom not in G:
					G[chrom] 	= list()
				G[chrom].append((int(start), int(stop), name+ "_" + strand))
			else:
				header=False
	for chrom in G:
		G[chrom] 	= node.tree(G[chrom])
	return G
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
	def __init__(self, k, ll , converged, diff):
		self.k 			= k
		self.ll 		= ll
		self.converged 	= converged
		self.diff 		= diff
		self.rvs 		= list()
	def insert_component(self, line):
		TYPE, rest 		= line.split(": ")
		if TYPE=="N":
			mu, si, l, w, pi 	= rest.strip("\n").split(",")
			mu, si, l, w, pi 	= float(mu), float(si), float(l), float(w), float(pi)
			component 			= EMG(mu, si, l, w, pi)
		elif TYPE=="U":
			a, b, w, pi 		= rest.strip("\n").split(",")
			a, b, w, pi 		= float(a), float(b), float(w), float(pi)
			component 			= UNI(a, b, w, pi)
		self.rvs.append(component)

	def pdf(self, x, st):
		return sum([rv.pdf(x, st) for rv in self.rvs])
	def __str__(self):
		return "~" + str(self.k) + "," + str(self.ll) + "," + str(self.converged) + "," + str(self.diff) + "\n"+ "\n".join([rv.__str__() for rv in self.rvs])+"\n"


class segment:
	def __init__(self, chrom, start, stop, N, annotation_N=0):
		self.annotations 		= list()
		self.chrom 				= chrom
		self.start 				= start
		self.stop 				= stop
		self.N 					= N
		self.models 			= {}
		self.k 					= None
		self.data_types 		= list()
		self.BINNED 			= False
		self.annotation_N 		= annotation_N
	def insert_model_info(self,line):
		k,ll,converged, diff 	= line[1:].strip("\n").split(",")
		k,ll,converged, diff 	= int(k), float(ll), bool(converged), float(diff)
		if k not in self.models:
			self.models[k] 		= list()
		self.k 					= k
		self.models[k].append(model(k,ll, converged, diff))
	def insert_component(self, line):
		self.models[self.k][-1].insert_component(line)
	def insert_data(self, x,y, data_type):
		if not hasattr(self, data_type):
			setattr(self, data_type, list())
			self.data_types.append(data_type)
		getattr(self, data_type).append((x,y))
	def bin_data_types(self, bins=500):
		self.BINNED=True
		for data_type in self.data_types:
			if data_type!= "dbSNP" and data_type != "ClinVar":
				LST 		= getattr(self, data_type)
				counts,edges= np.histogram([x for x,y in LST ],weights=[y for x,y in LST], bins=bins)
				edges 		= (edges[1:] + edges[:-1]) /2.
				edges 		-= self.start
				edges 		/= 100.
				setattr(self, data_type+"_binned", zip(edges, counts))
		
	def print_models(self):
		string 	= ""
		for k in self.models:
			for model in self.models[k]:
				string+=model.__str__()
		return string
	def print_segment_header(self):
		return "#" + self.chrom+":" + str(self.start) + "-" + str(self.stop) + "," + str(self.N) + "," + str(len(self.annotations)) + "\n"
	def print_data_types(self):
		if not self.BINNED:
			self.bin_data_types()
		string 	= ""
		for data_type in self.data_types:
			if data_type!= "dbSNP" and data_type != "ClinVar":
				string+= data_type+","+ str(getattr(self, data_type+ "_peak"))+ ","+  ":".join([str(x) + "," + str(y) for x,y in getattr(self, data_type+"_binned")]) + "\n"
			else:
				string+= data_type+","+ "1" + ","+  ":".join([str((x-self.start)/100.) + "," + y for x,y in getattr(self, data_type) ]) + "\n"
		return string
def load_file_segments(G, FILE, header=True, test=None): #this just adds the model fit parameters
	s 	= 0
	parameters 	= {}
	for line in FILE:
		if not header:
			if "#" == line[0]:
				chrom,info 		= line[1:].strip("\n").split(":")
				start_stop, N 	= info.split(",")
				start,stop 		= start_stop.split("-")
				if chrom not in G:
					G[chrom] 	= list()
				if test is not None and s > test:
					break
				G[chrom].append(segment(chrom,int(start),int(stop),float(N)))
				s+=1
			elif "~" == line[0]:
				G[chrom][-1].insert_model_info(line)
			else:
				G[chrom][-1].insert_component(line)
			
		else:
			data 		= line.strip("\n").split(",")
			parameters  = dict([d.split(": ") for d in data[2:]])
			header=False

def load_directory_segments(directory, test=None):
	G 	= {}
	for FILE in os.listdir(directory):
		load_file_segments(G, open(directory+FILE, "r"), test=test)
	return G

def load_data(FILE, G, data_type , test=None):
	if os.path.exists(FILE):	
		print "working on raw signal merger", data_type
		with open(FILE) as FH:
			prevChrom=""
			for line in FH:
				chrom,start, stop, coverage 	= line.strip("\n").split("\t")
				x 								= (float(stop) + float(start))  / 2.
				if chrom !=prevChrom:
					if chrom in G:
						j,N 	= 0, len(G[chrom])
					else:
						j,N 	= 0,0
				while j < N and G[chrom][j].stop < x:
					j+=1
				if j < N and G[chrom][j].start <= x:
					G[chrom][j].insert_data(x, float(coverage), data_type)	
					if test is not None:
						break

				prevChrom 	= chrom
	else:
		print FILE, "doesn't exist"
def classify_annotations(RF, G):
	for chrom in G:
		rf 	= RF[chrom]
		for I in G[chrom]:
			intervals 	= rf.searchInterval((I.start, I.stop))
			if intervals is not None:
				I.annotations 	= intervals
def load_peak_files(FILE, G, TYPE):
	if not FILE or not os.path.exists(FILE):
		for chrom in G:
			for I in G[chrom]:
				setattr(I, TYPE+"_peak", True)
	else:
		with open(FILE) as FH:
			A 		= {}
			for line in FH:
				chrom,start, stop 	= line.strip("\n").split("\t")[:3]
				start, stop 		= int(start), int(stop)
				if chrom not in A:
					A[chrom]=list()
				A[chrom].append((start, stop, ""))
		

		for chrom in G:
			if chrom in A:
				T  	= node.tree(A[chrom])
				for I in G[chrom]:
					F 	= T.searchInterval((I.start, I.stop))
					if F:
						setattr(I, TYPE+"_peak", True)
					else:
						setattr(I, TYPE+"_peak", False)	
			else:
				for I in G[chrom]:
					setattr(I, TYPE+"_peak", False)
	if not os.path.exists(FILE):
		print FILE, "peak file, doesn't exist"
						

def merge_files(G,FILES, TYPES, test=None):
	for FILE, TYPE in zip(FILES, TYPES):
		load_data(FILE, G, data_type=TYPE, test=test)
def merge_peak_files(G,FILES, TYPES, test=None):
	
	for FILE, TYPE in zip(FILES, TYPES):
		print "woking on peak merger", TYPE
		load_peak_files(FILE, G, TYPE)	

	
def insert_dbSNP(G, dbSNP_directory, data_type,test=None): #bed
	#need to make intervals into interval tree
	A  	= {}
	for chrom in G:
		A[chrom] = node.tree([(I.start, I. stop, I) for I in G[chrom]])
	for FILE in os.listdir(dbSNP_directory):
		if "bed" == FILE[-3:]:
			header 	= True
			with open(dbSNP_directory+FILE) as FH:
				for line in FH:
					if not header:
						chrom, start, stop, ID, zero, strand 	= line.strip("\n").split("\t")
						start, stop 	= int(start), int(stop)
						if chrom in A:
							finds 		= A[chrom].searchInterval((start, stop))
							if finds:
								for st, sp, I in finds:
									I.insert_data(start,ID, data_type)	
								if test is not None:
									break
					else:
						header=False
def insert_clinVarSNP(G, FILE, data_type, sig_level=4): #vcf
	A 	= {}
	N 	= 0.
	NN 	= 0.
	NNN = 0.
	with open(FILE) as FH:
		for chrom in G:
			A[chrom] = node.tree([(I.start, I. stop, I) for I in G[chrom]])
		for line in FH:
			if "#" != line[0]:
				chrom, pos, ID, REF, ALT,QUAL,FILTER, INFO 	= line.strip("\n").split("\t")

				chrom="chr"+chrom
				info_array 		= dict([d.split("=") for d in INFO.split(";") if len(d.split("="))==2 ])
				sig 	 		= int(info_array["CLNSIG"].split("|")[-1].split(",")[-1])
				if sig  >= sig_level and sig!=255:
					NN+=1
				if chrom in A:
					pos 		= int(pos )
					finds 	 	= A[chrom].searchPoint(pos)
					if finds:
						for st, sp, I in finds:
							I.insert_data(pos,ID, data_type)	
						NNN+=1
								
				N+=1


def write_out(G, out):
	FHW = open(out, "w")
	for chrom in G:
		for I in G[chrom]:
			if len(I.data_types) > 1:
				string 	= I.print_models()
				FHW.write(I.print_segment_header())
				FHW.write(I.print_models())
				FHW.write(I.print_data_types())
	FHW.close()









if __name__ == "__main__":
	TEST 					= None
	
	#========================================================================================================================
	#data file directories
	if len(sys.argv)==1:
		model_fits   		= "/Users/joeyazo/Desktop/Lab/gro_seq_files/HCT116/EMG_out_files/DMSO_ND_intervals_model_fits_2/"#on MAC
		refseq_file 		= "/Users/joeyazo/Desktop/Lab/genome_files/RefSeqHG19.txt"
		dbSNP_directory 	= "/Users/joeyazo/Desktop/Lab/dbSNP/"
		D 					= "/Users/joeyazo/Desktop/Lab/ENCODE/HCT116/"
		D2 					= "/Users/joeyazo/Desktop/Lab/gro_seq_files/HCT116/bed_graph_files/"
		out_dir 			= "/Users/joeyazo/Desktop/Lab/EMG_files/"
	else:
		model_fits   		= sys.argv[1] 
		refseq_file 		= sys.argv[2]
		dbSNP_directory 	= sys.argv[3]
		D 					= sys.argv[4]
		D2 					= sys.argv[5]
		out_dir 			= sys.argv[6]
		
	#========================================================================================================================
	#data files
	
	DNAse 				= D+"DNAse/bedgraph_files/wgEncodeUwDnaseHct116RawRep1.bedgraph"
	H3K27ac 			= D+"H3K27ac/bedgraph_files/ENCFF000VCU.bedgraph"
	H3K4me1 			= D+"H3K4me1/bedgraph_files/ENCFF000VCM.bedgraph"
	H3K4me3 			= D+"H3K4me3/bedgraph_files/ENCFF001FIM.bedgraph"
	gro_forward 		= D2+"DMSO2_3.pos.BedGraph"
	gro_reverse  		= D2+"DMSO2_3.neg.BedGraph"
	poll_II 			= D+"PolII/bedgraph_files/tPol_II_DMSO_150bp_genomeCovGraphPDMMR.BedGraph"
	CTCF 				= D+"CTCF/bedgraph_files/ENCFF000OYT.bedgraph"
	JunD 				= D+"JunD/bedgraph_files/ENCFF000PAA.bedgraph"
	Sin3A 				= D+"Sin3A/bedgraph_files/ENCFF000PBV.bedgraph"
	Sp1 				= D+"Sp1/bedgraph_files/ENCFF000PCF.bedgraph"
	ClinVar 			= dbSNP_directory+"clinvar.vcf"
	OUT_FILE 			= out_dir + "merged_data_file.txt"
	#========================================================================================================================
	#peak files 		
	DNAse_peak 			= D+"DNAse/peak_files/wgEncodeUwDnaseHct116PkRep1.narrowPeak"
	H3K27ac_peak 		= D+"H3K27ac/peak_files/HCT-116_H3K27Ac_narrowPeak.bed"
	H3K4me1_peak 		= D+"H3K4me1/peak_files/HCT-116_H3K4me1_narrowPeak.bed"
	H3K4me3_peak 		= D+"H3K4me3/peak_files/HCT116-DS16055.peaks.fdr0.01.hg19.bed"
	CTCF_peak 			= D+"CTCF/peak_files/SL12242_Peaks.bed.broadPeak"
	JunD_peak 			= D+"JunD/peak_files/SL12233_Peaks.bed.broadPeak"
	Sin3A_peak 			= D+"Sin3A/peak_files/SL12244_Peaks.bed.broadPeak"
	Sp1_peak 			= D+"Sp1/peak_files/SL12239_Peaks.bed.broadPeak"
	

	
	RF 					= load_refseq(refseq_file)
	print "loaded refseq annotations"
	#okay so local is way faster....
	G 					= load_directory_segments(model_fits, test=TEST)
	print "loaded model fits"

	insert_clinVarSNP(G, ClinVar, "ClinVar" )
	print "inserted ClinVar aata"

	insert_dbSNP(G, dbSNP_directory, "dbSNP",test=TEST)
	print "inserted SNP data"

	classify_annotations(RF, G)
	print "classified model fits by annotations"
	
	

	TYPES 				= ("DNAse", "H3K27ac", "H3K4me1", "H3K4me3", "gro_f", "gro_r", "poll","CTCF", "JunD", "Sin3A", "Sp1" )
	FILES 				= (DNAse, H3K27ac, H3K4me1,H3K4me3, gro_forward, gro_reverse, poll_II, CTCF, JunD, Sin3A,  Sp1 )
	peak_FILES 			= (DNAse_peak, H3K27ac_peak, H3K4me1_peak,H3K4me3_peak, "", "", "", CTCF_peak, JunD_peak, Sin3A_peak,  Sp1_peak )

	merge_peak_files(G,peak_FILES, TYPES, test=None)

	merge_files(G, FILES, TYPES,test=TEST)
	
	#write out
	write_out(G, OUT_FILE)	
















