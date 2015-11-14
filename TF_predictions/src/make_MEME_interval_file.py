import sys, numpy as np 
import node, os
def collect_all_ChIP_motif_hits(FILES, FHW,i):
	for CM in FILES:
		G 	= {}
		T 	= 0
		header 	= True
		with open(CM) as FH :
			for line in FH:
				if not header:
					line_array 	= line.strip("\n").split("\t")
					chrom,start, stop 	= line_array[1],line_array[2],line_array[3]
					start, stop 		= int(start),int(stop)
					FHW.write(line_array[1]+"\t" +line_array[2]+"\t"+line_array[3]+"\tCM_"+ CM +","+ str(i)+"\n" )
					if chrom not in G:
						G[chrom]=list()
					G[chrom].append((start-1000, stop+1000))

					T+=1
					i+=1
				else:
					header=False
		for chrom in G:
			G[chrom].sort()
			G[chrom]=node.tree(G[chrom])
	return G,i



def write_filter_etc(MO,FHW, i, MODEL, G):
	t=0
	header=True
	T 		= 5000
	with open(MO) as FH:
		for line in FH:
			if not header:
				line_array 	= line.strip("\n").split("\t")
				chrom,start, stop 	= line_array[1],line_array[2],line_array[3]
				start, stop 		= int(start),int(stop)
				U 					= np.random.uniform(0,1)
				if U < 0.1 and chrom in G and not G[chrom].searchInterval((start, stop)) and t <T :
					FHW.write(line_array[1]+"\t" +line_array[2]+"\t"+line_array[3]+"\tMO_" + MODEL+","+str(i)+"\n")
					t+=1
					i+=1
				if t > T:
					break
			else:
				header=False
	return i


def iterate(root, out):
	for TF_DIR in os.listdir(root):
		PATH 	= root+"/"+TF_DIR+"/peak_files/outfiles/MEME/"
		FHW 	= open(out+TF_DIR, "w")
		i 		= 1
		if os.path.exists(PATH):
			CM_FILES 	= [PATH+"/" +fimo_dir + "/fimo.txt" for fimo_dir in os.listdir(PATH) if fimo_dir[:8]=="fimo_out" ]
			G,i 			= collect_all_ChIP_motif_hits(CM_FILES, FHW,i)
			print TF_DIR
			for fimo_dir in os.listdir(PATH):
				if  os.path.exists(PATH+fimo_dir.split("_")[-1]+ "_fimo_out"): 
					MOITF_ONLY 	= PATH+"/" +fimo_dir.split("_")[-1]+ "_fimo_out/fimo.txt"
					MODEL 		= fimo_dir.split("_")[-1]
					i 			= write_filter_etc(MOITF_ONLY , FHW, i, MODEL, G)
		FHW.close()



if __name__ == "__main__":
	root 	= sys.argv[1]
	out 	= sys.argv[2]
	iterate(root, out)