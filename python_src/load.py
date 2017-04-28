import node
import random
import numpy as np
import matplotlib.pyplot as plt
import simulate,re
import time
import model_across

class info:
	def __init__(self, st, sp, unit, chrom=""):
		self.tot_st, self.tot_sp	= st, sp
		self.units 					= [unit]
		self.forward, self.reverse 	= list(), list()
		self.X  					= None
		self.okay 					= False
		self.chrom 					= chrom
		self.N 						= 0
	def update_bounds(self, st, sp, unit):				
		if st!=self.tot_st and sp != self.tot_sp:
			self.tot_st, self.tot_sp	= min(self.tot_st, st), max(self.tot_sp, sp)
			self.units.append(unit)
	def check(self):
		if len(self.forward) and len(self.reverse):
			self.okay 	= True 
		return self.okay
	def insert(self, x,y, strand):
		if strand == 0:
			self.forward.append((x,y))
		else:
			self.reverse.append((x,y))
	def bin(self, bins=300):
		self.check()
		if self.okay:
			self.X 		= simulate.BIN(self.forward, self.reverse, bins)
			self.N 		= np.sum(self.X[:,1:])
		self.forward, self.reverse 	= None, None

	def show(self):
		if self.okay:
			F 	= plt.figure(figsize=(15,10))
			ax 	= F.add_subplot(1,1,1)
			ax.set_title(str(self.tot_st) + "-" + str(self.tot_sp))
			ax.bar(self.X[:,0], self.X[:,1])
			ax.bar(self.X[:,0], -self.X[:,2])
			plt.show()	
	def print_info(self):
		txt 	= self.chrom + ":" + str(self.tot_st) + "-" + str(self.tot_sp) + ":" + str(self.N) + "\n"
		return txt
	def __str__(self):
		txt 	= self.chrom + ":" + str(self.tot_st) + "-" + str(self.tot_sp) + ":" + str(self.N) 
		return txt

def grab_specific_region(chrom_spec,start_spec, stop_spec, 
	pos_file="/Users/joazofeifa//Lab/gro_seq_files/HCT116/bed_graph_files/DMSO2_3.pos.BedGraph", 
	neg_file="/Users/joazofeifa//Lab/gro_seq_files/HCT116/bed_graph_files/DMSO2_3.neg.BedGraph",
	SHOW 	=False, bins=None):
	D 		= (list(), list())
	files 	= (pos_file, neg_file)
	for i,f in enumerate(files):
		FOUND 	=False
		with open(f) as FH:
			for line in FH:
				chrom,start, stop, coverage 	= re.split("\s+", line.strip("\n"))
				if chrom == chrom_spec:
					FOUND=True
				start, stop, coverage 			= float(start), float(stop), abs(float(coverage))
				if FOUND and (chrom!=chrom_spec or start>stop_spec ):
					break
				if FOUND and start < stop_spec and stop > start_spec:
					for j in range(int(start), int(stop)):
						D[i].append((j, coverage))
	if bins is None:
		return D
	if SHOW:
		X 	= simulate.BIN(D[0], D[1], 200)
		plt.title(chrom + ":" + str(start_spec) + "-" + str(stop_spec))
		plt.bar(X[:,0], X[:,1])
		plt.bar(X[:,0], -X[:,2])
		plt.show()
	X 		= simulate.BIN(D[0], D[1], bins)
	return X

def merge_intervals(G):
	IS 				= {}
	for chrom in G:
		G[chrom].sort()
		IS[chrom] 	= list()
		st, sp  	= G[chrom][0][0], G[chrom][0][1]
		I 			= info(st, sp, G[chrom][0])
		i 			= 1
		N 			= len(G[chrom])
		while i < N:
			while i < N and G[chrom][i][0] < I.tot_sp and G[chrom][i][1] > I.tot_st:
				I.update_bounds(G[chrom][i][0], G[chrom][i][1], G[chrom][i])
				i+=1
			IS[chrom].append(I)
			if i < N:
				I 	= info(G[chrom][i][0], G[chrom][i][1], G[chrom][i])
			i+=1
	return IS


def gene_annotations(FILE, SI=True, pad=0):
	G  	= {}
	with open(FILE) as FH:
		header 	= True
		for line in FH:
			if not header:
				bin, name, chrom, strand, start, stop 	= line.split("\t")[:6]
				strand 	= int(strand=="-")
				if chrom not in G:
					G[chrom]=list()
				G[chrom].append((int(start)-pad, int(stop)+pad, strand, name))
			else:
				header=False
	#merge?
	IS 		= merge_intervals(G)
	if SI: #filter for only single isoform genes
		single 	= {}
		N 		= sum([len(ISS) for ISS in IS.values()])
		t 		= 0
		for chrom in IS:
			for I in IS[chrom]:
				if len(I.units)==1:
					if chrom not in single:
						single[chrom]=list()
					single[chrom].append(I)			
					t+=1
	
	return IS
def FStitch_annotations(forward, reverse, merge=True, pad=0):
	G 	= {}
	for i,f in enumerate((forward, reverse)):
		header=True
		with open(f) as FH:
			for line in FH:
				if not header:
					chrom,start, stop, state 	= line.split("\t")[:4]
					if "ON" in state:
						if chrom not in G:
							G[chrom] 			= list()
						G[chrom].append((int(start)-pad, int(stop)+pad, i))
				else:
					header=False
	if not merge:
		return G
	IS 	= merge_intervals(G)
	return IS

def filter_single_overlaps(FS, RF):
	filtered	= {}
	for chrom in FS:
		if chrom in RF:
			fs,rf 	= FS[chrom], RF[chrom]
			i,N 	= 0, len(RF[chrom])
			for I in fs:
				while i < N and rf[i].tot_sp<I.tot_st:
					i+=1
				ot 	= 0
				j 	= i
				while j < N and rf[j].tot_st<I.tot_sp:
					o_st, o_sp 	= max(rf[j].tot_st, I.tot_st), min(rf[j].tot_sp, I.tot_sp)
					ot+=1
					j+=1
				if ot ==1:
					if chrom not in filtered:
						filtered[chrom] 	= list()
					if float(o_sp - o_st) / float(rf[j-1].tot_sp - rf[j-1].tot_st) > 0.5:
						filtered[chrom].append((I, j))
	final 			= {}
	for chrom in filtered:
		final[chrom] 	= list()
		for j,(I,i) in enumerate(filtered[chrom]):
			if j == 0 or j == len(filtered[chrom])-1:
				final[chrom].append(I)
			elif i!=filtered[chrom][j-1][1] and i!=filtered[chrom][j+1]:
				final[chrom].append(I)
	return final
def bedGraphFile(forward_file, reverse_file,intervals, write_out=True, test=True,SHOW=False):
	for s,FILE in enumerate((forward_file, reverse_file)):
		prev, i 	= "", 0
		print "working on", FILE
		with open(FILE) as FH:
			for line in FH:
				chrom,start, stop, coverage 	= re.split("\s+", line.strip("\n"))
				start, stop, coverage 			= int(start), int(stop), abs(float(coverage))
				if test and i>10:
					break
					
				if chrom != prev:
					i 	= 0
					if chrom in intervals:
						IS 	= intervals[chrom]
						N 	= len(IS)
					else:
						N       = 0
				while i < N and IS[i].tot_sp < start:
					i+=1
				if i < N and IS[i].tot_st < stop:
					o_st,o_sp 	= max(start, IS[i].tot_st), min(stop, IS[i].tot_sp)
					for o in range(o_st, o_sp):
						IS[i].insert(o,coverage, s)
				prev 	= chrom
	print "finished reading in FILES"
	if write_out is not None:
		print "writing out intervals/data for future parsing"
		print write_out
		out_file 	= write_out+"/out_format_file_"

		out_file 	= model_across.checkFileExists(out_file, 1)
		print out_file
		FHW			= open(out_file, "w")
	for chrom in intervals:
		for I in intervals[chrom]:
			if SHOW:
				I.show()
			if write_out is not None and I.check():
				FHW.write("#" + chrom + "," + str(I.tot_st) + "," + str(I.tot_sp) + "\n")
				FHW.write("~forward\n")
				for x,y in I.forward:
					FHW.write(str(x) + "," + str(y)  + "\n")
				FHW.write("~reverse\n")
				for x,y in I.reverse:
					FHW.write(str(x) + "," + str(y)  + "\n")
					
	if write_out is not None:
		FHW.close()
def formatted_file(FILE, bins, sc):
	I  	= None	
	D 	= list()
	with open(FILE) as FH:
		for i,line in enumerate(FH):

			if line[0]=="#":
				if I is not None:
					I.bin(bins=bins)
					D.append(I)
				
				chrom, st, sp 	= line[1:].strip("\n").split(",")

				st, sp 	= int(st), int(sp)

				if chrom == sc or sc == "all":
					I 		= info(st, sp, None, chrom=chrom)
				else:
					I 		= None
				forward, reverse 	= False, False
			elif "forward" in line:
				forward = True
				reverse = False
			elif "reverse" in line:
				forward = False
				reverse = True
			elif forward and I is not None:
				x,y 	= line.strip("\n").split(",")
				x,y 	= float(x), float(y)
				I.insert(x, y, 0)
			elif reverse and I is not None:
				x,y 	= line.strip("\n").split(",")
				x,y 	= float(x), float(y)
				I.insert(x, y, 1)
	return D
























