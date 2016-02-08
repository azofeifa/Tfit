import math
import matplotlib.pyplot as plt
import time
def load(FILE):
	G 	= {}
	with open(FILE) as FH:
		for line in FH:
			if line[0]!="#":
				chrom,start, stop, pars 	= line.strip("\n").split("\t")
				if chrom not in G:
					G[chrom]  = list()
				pvs 			= [float(x) for x in pars.split("|")[1].split(",")]
				if pvs[1] < 300 and 500 > pvs[2] > 50:
					G[chrom].append((int(start), int(stop), [float(x) for x in pars.split("|")[1].split(",")]))
	for chrom in G:
		G[chrom].sort()
	return G
def match_up(A,B):
	O 	= {}
	for chrom in A:
		if chrom in B:
			a,b 	= A[chrom], B[chrom]
			j,N 	= 0,len(b)
			O[chrom]=list()
			for start, stop, ps in a:
				while j < N and b[j][1] < start:
					j+=1
				if j < N and b[j][0] < stop:
					O[chrom].append((ps, b[j][2]))
	return O


def scatter(O):
#	plt.scatter([math.log(x[8],10) for chrom in O for x,y in O[chrom]],[math.log(y[8],10) for chrom in O for x,y in O[chrom]],s=5,alpha=0.5)
	print len([y[3] for chrom in O for x,y in O[chrom]])
	plt.scatter([x[5] for chrom in O for x,y in O[chrom]],[y[5] for chrom in O for x,y in O[chrom]],s=5,alpha=0.5)
	plt.show()



if __name__ == "__main__":
	FILE1	= "/Users/joazofeifa/Lab/gro_seq_files/Rubin2016/Rubin2015_DMSO-1_divergent_classifications.bed"
	FILE2 	= "/Users/joazofeifa/Lab/gro_seq_files/Rubin2016/Rubin2015_DMSO-2_divergent_classifications.bed"
	A 		= load(FILE1)
	B 		= load(FILE2)
	O 		= match_up(A,B)
	scatter(O)