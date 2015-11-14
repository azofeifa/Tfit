import sys, numpy as np 
import node, os
def write_filter_etc(CM, MO,FHW, i):
	G 	= {}
	T 	= 0
	with open(CM) as FH :
		for line in FH:
			line_array 	= line.strip("\n").split("\t")
			chrom,start, stop 	= line_array[1],line_array[2],line_array[3]
			start, stop 		= int(start),int(stop)
			FHW.write(line_array[1]+"\t" +line_array[2]+"\t"+line_array[3]+"\tCM_" + str(i)+"\n" )
			if chrom not in G:
				G[chrom]=list()
			G[chrom].append((start-1000, stop+1000))
			T+=1
			i+=1
	for chrom in G:
		G[chrom].sort()
		G[chrom]=node.tree(G[chrom])
	t=0
	with open(MO) as FH:
		for line in FH:
			line_array 	= line.strip("\n").split("\t")
			chrom,start, stop 	= line_array[1],line_array[2],line_array[3]
			start, stop 		= int(start),int(stop)
			U 					= np.random.uniform(0,1)
			if U < 0.1 and not G[chrom].searchInterval((start, stop)) and t <T :
				FHW.write(line_array[1]+"\t" +line_array[2]+"\t"+line_array[3]+"\tMO_" + str(i)+"\n")
				t+=1
				i+=1
	return i


def iterate(root, out):
	for TF_DIR in os.listdir(root):
		PATH 	= TF_DIR+"/peak_files/outfiles/MEME/"
		FHW 	= open(out+TF_DIR, "w")
		i 		= 1
		if os.path.exists(PATH):
			for fimo_dir in os.listdir(PATH):
				print fimo_dir
				if fimo_dir[:8]=="fimo_out" and os.path.exists(fimo_dir.split("_")[-1]+ "_fimo_out"): 
					ChIP_MOTIF 	= PATH+"/" +fimo_dir + "/fimo.txt"
					MOITF_ONLY 	= PATH+"/" +fimo_dir.split("_")[-1]+ "_fimo_out/fimo.txt"
					i 			= write_filter_etc(ChIP_MOTIF, MOITF_ONLY , FHW, i)
		FHW.close()



if __name__ == "__main__":
	root 	= sys.argv[1]
	out 	= sys.argv[2]
	iterate(root, out)