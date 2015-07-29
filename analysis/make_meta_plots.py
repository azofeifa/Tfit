import numpy as np
import matplotlib.pyplot as plt
import node
import merge_data_types as mdt
import time
def load(FILE, w_thresh=0.01, s_thresh=2, l_thresh=2., pi_thresh=0.05, PEAK=False):
	G 	= {}
	N 	= 0
	with open(FILE) as FH:
		collect=False
		for line in FH:
			if "#" == line[0]:
				collect = False
			elif "N: " == line[:3]:
				mu, si, l, w, pi 	= line[3:].strip("\n").split(",")
				mu, si, l, w, pi 	= float(mu), float(si), float(l), float(w), float(pi)
				if si < s_thresh and l < l_thresh and pi_thresh<= pi<=(1-pi_thresh):
					collect=True
			elif collect:
				lineArray = line.strip("\n").split(",")
				data_type, peak, data 	= lineArray[0], lineArray[1], lineArray[2:]
				if data_type not in G:
					G[data_type] 	= np.zeros((len(data), ))
				NN 	= sum(np.array([float(d.split("-")[1]) for d in data ]))
				if NN:
					G[data_type]+=((np.array([float(d.split("-")[1]) for d in data ]))/NN)
				N+=1
	for data_type in G:
		F 	= plt.figure(figsize=(15,10))
		plt.title(data_type)
		plt.plot(np.linspace(-4,4, len(G[data_type])-2), G[data_type][2:] , label="Meta Signal", linewidth=3. )
		plt.fill_between(np.linspace(-4,4, len(G[data_type])-2), np.zeros( (len(G[data_type])-2, ) ), G[data_type][2:],alpha=0.5, color="grey" )
		plt.scatter([0], [0], marker=(5,2),s=1000, label="Predicted Position of Pol Loading", color="blue")
		plt.xlabel("Standard Deviations (Genomic Coordinate)")
		plt.ylabel("Averaged Signal")
		plt.legend()
		plt.grid()
		plt.show()
class segment:
	def __init__(self, chrom ,start ,stop):
		self.chrom 	= chrom
		self.start 	= start
		self.stop 	= stop
		self.SNPS 	= list()
		self.bidirs = list()
		self.metas 	= list()
		self.normed = False
	def insert_model(self, mu, si, l):
		self.bidirs.append(( mu, si + (1. / l ) ))
	def bin(self, bins):
		if not self.normed:
			self.SNPS 	= [(x-self.start)/100. for x in self.SNPS]
			self.normed = True
		for mu, std in self.bidirs:

			X 		= np.zeros((bins, 2))
			X[:,0] 	= np.linspace(-4, 4, bins)*std + mu
			j 		= 0

			for pos in self.SNPS:
				while j < bins and X[j,0] < pos:
					j+=1
				if 1< j < bins:
					X[j-1,1]+=1
			self.metas.append(X)

			


def load2(FILE,RF, w_thresh=0.01, s_thresh=2, l_thresh=2., pi_thresh=0.05, PEAK=False):
	with open(FILE) as FH:
		G 		= {}
		t 		= 0
		for line in FH:
			if "#" == line[0]:
				chrom,stsp 	= line[1:].strip("\n").split(":")
				start, stop = [int(x) for x in stsp.split("-")]
				S 		= segment(chrom, start, stop)
				if chrom not in G:
					G[chrom] 	= list()
				collect = True
				G[chrom].append(S)
			

			elif "N: " == line[:3] and collect:
				mu, si, l, w, pi 	= line[3:].strip("\n").split(",")
				mu, si, l, w, pi 	= float(mu), float(si), float(l), float(w), float(pi)
				if si < s_thresh and l < l_thresh and pi_thresh<= pi<=(1-pi_thresh):
					sst, ssp 	= mu*100 + start - 100*(si + 1.0/ l)*1,mu*100 + start + 100*(si + 1.0/ l)*1

					if  RF[chrom].searchInterval((sst, ssp)):
						G[chrom][-1].insert_model(mu, si, l)

	return G



def insert_clinvar(FILE, G):
	A 	= {}
	with open(FILE) as FH:
		for line in FH:
			if "#"!=line[0]:
				chrom, pos 	= line.split("\t")[:2]
				chrom 		= "chr" + chrom
				pos 		= int(pos)
				if chrom not in A:
					A[chrom] 	= list()
				A[chrom].append(pos)
	for chrom in A:
		if chrom in G:

			a,g 	= A[chrom], G[chrom]
			a.sort()
			N,j 	= len(g),0
			for pos in a:
				while j < N and g[j].stop < pos:
					j+=1
				if j < N and  g[j].start<= pos <= g[j].stop:
					g[j].SNPS.append(pos)
	bins 	= 75
	counts 	= np.zeros((bins,))
	for chrom in G:
		for S in G[chrom]:
			S.bin(bins=bins)
			for i,X in enumerate(S.metas):
				if sum(X[:,1]):
					counts+=X[:,1]

	##plot
	plt.figure(figsize=(15,10))
	plt.title("SNP Density")
	plt.plot(np.linspace(-4,4, bins-5), counts[5:])
	plt.fill_between(np.linspace(-4,4, len(counts)-5), np.zeros( (len(counts)-5, ) ), counts[5:],alpha=0.5, color="grey" )
		
	plt.xlabel("Standard Deviations (Genomic Coordinate)")
	plt.ylabel("SNP Density")
	plt.scatter([0], [0], marker=(5,2),s=1000, label="Predicted Position of Pol Loading", color="blue")
		
	plt.legend()
	
	plt.grid()
	plt.show()




if __name__ == "__main__":
	READ 	= True
	READ2 	= False

	if READ2:
		FILE 	= "/Users/joeyazo/Desktop/Lab/EMG_files/counted_data_100_5_4_2_0.01_0.1"
		FILE2 	= "/Users/joeyazo/Desktop/Lab/dbSNP/clinvar.vcf"
		FILE3 	= "/Users/joeyazo/Desktop/Lab/genome_files/RefSeqHG19.txt"
		RF 		= mdt.load_refseq(FILE3)
		G 		= load2(FILE, RF)
		
		insert_clinvar(FILE2, G)

	
	if READ:
		FILE 	= "/Users/joeyazo/Desktop/Lab/EMG_files/counted_data_100_5_4_2_0.01_0.1"
		load(FILE)
