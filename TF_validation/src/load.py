import time
class bidir:
	def __init__(self, line,chrom,start, stop, params, TFS, raw):
		self.line 		= line
		self.chrom 		= chrom
		self.start 		= int(start)
		self.stop  		= int(stop)
		self.partner 	= None
		if raw is not None and raw:
			self.raw 		= dict([(tf.split(":")[0], float(tf.split(":")[1] )) for tf in raw.split(",") ])
		else:
			self.raw 	= {}
		self.TFS 		= {}
		self.CHIP 		= {}
		self.rawCHIP 	= {} 
		self._construct_TFS(TFS)
		self._construct_params(params)
	def add_coverage(self, TF, y ):
		if TF not in self.rawCHIP:
			self.rawCHIP[TF]=0
		self.rawCHIP[TF]+=y
	def get_density(self):
		return self.N / (self.stop - self.start)
	def _construct_TFS(self,tfs):
		if tfs:
			for TF in tfs.split(","):
				tf 	= TF.split("=")[0]
				self.TFS[tf]=list()
				for instance in TF.split("=")[1].split(";"):
					dist, pval, qval,strand 	= instance.split("_")
					dist, pval, qval 			= float(dist), float(pval), float(qval)
					self.TFS[tf].append((dist, pval, qval, strand))




	def _construct_params(self,params):
		if params:
			#100029024.056121_32.578913_20.000000_0.340978_0.888104_34927_350_-522981.937500
			mu,si,l,w,pi,N,fp,nll 	= [float(x) for x in params.split("_")]
			self.mu,self.si,self.l,self.w,self.pi,self.N,self.fp,self.nll 	= mu,si,l,w,pi,N,fp,nll
def read_in_motifs(FILE, window=0, RAW=False,test=False):
	G,M={},{}
	with open(FILE) as FH:
		header=True
		t=0
		for line in FH:
			if test and t > 1000:
				break
			t+=1
			if not header:
				start, stop 	= None, None
				if not RAW:
					chrom,start, stop 	= line.strip("\n").split("\t")[:3]
					raw 				= None
				elif RAW and len(line.split("\t")) > 3 :
					chrom,start, stop, raw 	= line.strip("\n").split("\t")
				elif RAW:
					TF, info 	= line.strip("\n").split("\t")
					M[TF] 		= float(info.split(",")[0]),float(info.split(",")[1])
				if start is not None:
					if chrom not in G:
						G[chrom]=list()
					x 	= (float(stop) + float(start)) /2.
					start 	= x  - (window / 2.)
					stop 	= x  + (window/ 2.)
					G[chrom].append((int(start), (bidir(line.strip("\n"), chrom, start, stop, "", "", raw ))) )
			else:
				header=False
	for chrom in G:
		G[chrom].sort()
		G[chrom]=[y for x,y in G[chrom]]
	if RAW:
		return G,M
	return G
def load_bidir(FILE, RAW=True, window=0,test=False):
	G,M 	= {},{}
	collect=False
	with open(FILE) as FH:
		t=0
		for line in FH:
			start, stop 	= None, None
			if test and t > 1000:
				break
			t+=1
			if not RAW and len(line.split("\t")) > 4:
				chrom,start, stop, params, TFS 	= line.strip("\n").split("\t")
				raw = None
			elif RAW and len(line.split("\t")) > 3 :
				chrom,start, stop, params, TFS, raw 	= line.strip("\n").split("\t")[:6]
			elif RAW:
				TF, info 	= line.strip("\n").split("\t")
				M[TF] 		= float(info.split(",")[0]),float(info.split(",")[1])

			if start is not None:
				x 	= (float(stop) + float(start)) /2.
				if window:
					start 	= x - (window / 2.)
					stop 	= x+ (window/ 2.)
				if chrom not in G:
					G[chrom]=list()
				G[chrom].append((int(start), (bidir(line.strip("\n"), chrom, start, stop, params, TFS, raw ))))
		
	for chrom in G:
		G[chrom].sort()
		G[chrom]=[ y for x,y in G[chrom]]
	if RAW:
		return G,M
	return G
def load_peaks(FILE):
	G	= {}
	with open(FILE) as FH:
		for line in FH:
			chrom,start, stop, 	= line.split("\t")[:3]
			if chrom not in G:
				G[chrom]=list()
			G[chrom].append((int(start)-1000, int(stop)+1000))
	for chrom in G:
		G[chrom].sort()
	return G
def overlap(BIDIR, CHIP, TF=""):
	chromosomes 	= BIDIR.keys()
	for chrom in chromosomes:
		if chrom in CHIP:
			a,b 	= BIDIR[chrom], CHIP[chrom]
			j,N 	= 0,len(b)
			for bid in a:
				while j < N and b[j][1] < bid.start:
					j+=1
				while j < N and b[j][0] < bid.stop:
					if TF not in bid.CHIP:
						bid.CHIP[TF]=0
					bid.CHIP[TF]+=1
					j+=1
		else:
			del BIDIR[chrom]
def add_ChIP_signal(B, FILES, OUT=""):
	FHW 	= open(OUT, "w")
	GLOBAL 	= {}
	for lbl,FILE in FILES:
		t=0
		L 		= 0
		GLOBAL[lbl]=[]
		prevstart=None
		TOTAL 	= 0
		print lbl
		with open(FILE) as FH:
			prevchrom=""
			prevstop=None
			for line in FH:
				chrom,start, stop, cov 	= line.strip("\n").split("\t")
				x,cov 					= (float(stop) + float(start) ) / 2., float(cov)
				cov*=(float(stop)-float(start))
				TOTAL+=cov
				if prevchrom!= chrom:
					print chrom,
					if prevstart is not None:

						L+=(prevstop-prevstart)
					prevstart=x
					if chrom in B:
						j,N 			 	= 0, len(B[chrom])
					else:
						j,N 				= 0,0
				while j < N and B[chrom][j].stop < x:
					j+=1
				k 	= j
				while k < N and B[chrom][k].start < x:
					B[chrom][k].add_coverage(lbl, cov)
					k+=1
				t+=1
				prevchrom 				= chrom
				prevstop 				= x
		GLOBAL[lbl] 	= (L,TOTAL)
		print 
	for lbl in GLOBAL:
		FHW.write(lbl+"\t" + str(GLOBAL[lbl][0]) + ","+ str(GLOBAL[lbl][1]) + "\n"  )
	for chrom in B:
		for b in B[chrom]:
			line 	= b.line
			line= line + "\t" +  ",".join([tf +":"+str(b.rawCHIP[tf])   for tf in b.rawCHIP ])
			FHW.write(line+"\n")
	FHW.close()

if __name__ == "__main__":
	make_raw 			= False
	DMSO2_3_FILE 		= "/Users/joazofeifa/Lab/TF_predictions/distances_crude/Allen2014_motif_distances.tsv"
	
	MAX_FIMO_FILE 		= "/Users/joazofeifa/Lab/TF_predictions/FIMO_OUT/MAX/MAX_f1_fimo_out.bed"
	SP1_FIMO_FILE 		= "/Users/joazofeifa/Lab/TF_predictions/FIMO_OUT/Sp1/SP1_f1_fimo_out.bed"
	
	DIR 				= "/Users/joazofeifa/Lab/TF_predictions/raw_TF_files/"
	
	rawMAX, rawSP1,rawZBTB33 	= DIR+"MAX_rep1.bedgraph", DIR+"Sp1_rep1.bedgraph", DIR+"ZBTB33_rep1.bedgraph"

	#DMSO2_3 		= load_bidir(DMSO2_3_FILE, RAW=False, window=1500)
	#MAX 			= read_in_motifs(MAX_FIMO_FILE, window=1500)
	#SP1 			= read_in_motifs(SP1_FIMO_FILE, window=1500)
	
	print "loaded files"
	OUT1 			= "/Users/joazofeifa/Lab/TF_predictions/distances_crude/Allen2014_rawCHIP_motif_distances.tsv"
		
	OUT2 			= "/Users/joazofeifa/Lab/TF_predictions/FIMO_OUT/MAX/fimo_rawCHIP.txt"
	OUT3 			= "/Users/joazofeifa/Lab/TF_predictions/FIMO_OUT/Sp1/fimo_rawCHIP.txt"
	
	if make_raw:
		OUT4 			= "/Users/joazofeifa/Lab/TF_predictions/FIMO_OUT/ZBTB33/fimo_rawCHIP.txt"
			
		#add_ChIP_signal(DMSO2_3, (("MAX", rawMAX), ("SP1", rawSP1 ), ("ZBTB33",rawZBTB33 )) , OUT=OUT1)
		DMSO2_3 			= None
		print "done :)"
	
		add_ChIP_signal(MAX, (("MAX", rawMAX), ("SP1", rawSP1 ) ) , OUT=OUT2)
		MAX 			= None
		print "done :)"
		add_ChIP_signal(SP1, (("MAX", rawMAX), ("SP1", rawSP1 ) ) , OUT=OUT3)
		SP1 			= None
		print "done :)"
	DMSO2_3 	= load_bidir(OUT1, RAW=True)
	MAX 		= read_in_motifs(OUT2, RAW=True)

		
		
	


