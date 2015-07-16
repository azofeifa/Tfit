import load,time
def write_out(IS, OUT):
	FHW 	= open(OUT, "w")
	for chrom in IS:
		for I in IS[chrom]:
			try:
				FHW.write(chrom + "\t" + str(I.tot_st) + "\t" + str(I.tot_sp) + "\n")
			except:
				FHW.write(chrom + "\t" + str(I.start) + "\t" + str(I.stop) + "\n")
				
	FHW.close()

class info:
	def __init__(self, chrom, start, stop):
		self.chrom 	= chrom
		self.start 	= start
		self.stop 	= stop
	def update_bounds(self, st, sp):				
		if st!=self.start and sp != self.stop:
			self.start, self.stop	= min(self.start, st), max(self.stop, sp)
	
def load(FILE):
	with open(FILE) as FH:
		header 	= True
		G 		= {}
		for line in FH:
			if not header:
				chrom, start, stop, state 	= line.strip("\n").split("\t")[:4]
				if "ON" in state:
					I 	= info(chrom, int(start), int(stop))
					if chrom not in G:
						G[chrom] 	= list()
					G[chrom].append(I)
			else:
				header=False
	return G

def merge_all(FILE_1,FILE_2,FILE_3,FILE_4):
	A,B,C,D 	= load(FILE_1),load(FILE_2),load(FILE_3),load(FILE_4)
	G 			= {}
	for chrom in A:
		if chrom in C and chrom in B and chrom in D:
			a,b,c,d 	= A[chrom], B[chrom], C[chrom], D[chrom]
			G[chrom] 	= list()
			for l in (a,b,c,d):
				for I in l:
					G[chrom].append((I.startN, I))
			G[chrom].sort()
	for chrom in G:
		G[chrom] 	= [I for i, I in G[chrom]]
	IS 	= {}
	for chrom in G:
		I 	= info(chrom, G[chrom][0].start, G[chrom][0].stop)
		i 	= 1
		A 	= G[chrom]
		N 	= len(A)
		IS[chrom] 	= list()
		while i < N:
			while i < N and A[i].start < I.stop and A[i].stop > I.start:
				I.update_bounds(A[i].start, A[i].stop)
				i+=1
			IS[chrom].append(I)
			if i < N:
				I 	= info(chrom, A[i].start, A[i].stop)
			i+=1
	return IS





	pass


if __name__=="__main__":

	#==================================================================
	# need to figure out WHAT kind of 
	# intervals we want to run the model on
	# Merged FStitch FILES is really the main thing
	# I think we can do all the filtering after the fact...

	DIR 			 	= "/Users/joeyazo/Desktop/Lab/gro_seq_files/HCT116/FStitch/"
	DMSO,NUTLIN,DMSO_NUTLIN 		= False, False,True
	if NUTLIN:
		FS_forward 			= DIR+"Nutlin2_3.sorted.fiveprime.pos_segs_IGV.bed"
		FS_reverse 			= DIR+"Nutlin2_3.sorted.fiveprime.neg_segs_IGV.bed"
		OUT 				= DIR+"Nutlin2_3_merged_FS.bed"
		FS 					= load.FStitch_annotations(FS_forward, FS_reverse, merge=True)
		write_out(FS, OUT)
	if DMSO:
		FS_forward 			= DIR+"DMSO2_3.sorted.fiveprime.pos_segs_IGV.bed"
		FS_reverse 			= DIR+"DMSO2_3.sorted.fiveprime.neg_segs_IGV.bed"
		OUT 				= DIR+"DMSO2_3_merged_FS.bed"
		FS 					= load.FStitch_annotations(FS_forward, FS_reverse, merge=True)
		write_out(FS, OUT)
	if DMSO_NUTLIN:#merge all!
		FILE_1 			= DIR+"Nutlin2_3.sorted.fiveprime.pos_segs_IGV.bed"
		FILE_2 			= DIR+"Nutlin2_3.sorted.fiveprime.neg_segs_IGV.bed"
		FILE_3 			= DIR+"DMSO2_3.sorted.fiveprime.pos_segs_IGV.bed"
		FILE_4 			= DIR+"DMSO2_3.sorted.fiveprime.neg_segs_IGV.bed"
		OUT 			= DIR+"DMSO2_3_Nutlin2_3_merged_FS.bed"
		
		IS 				= merge_all(FILE_1, FILE_2, FILE_3, FILE_4)
		write_out(IS, OUT)



