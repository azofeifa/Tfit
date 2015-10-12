import os, node
import numpy as np 
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt
import time
import sys,os
import matplotlib.patches as patches
sys.path.append("/Users/joazofeifa/Lab/")

from interval_package import intervals
import numpy as np
def write_out(G,OUT):
	FHW 	= open(OUT, "w")
	for chrom in G:
		for start, stop, INFO in G[chrom]:
			FHW.write(chrom+"\t" + str(start) + "\t" + str(stop) + "\t" + INFO + "\n")
	FHW.close()
def bed_file(FILE, FILTER=None):
	G 	= {}
	ct 	= 0.
	N 	= 0.
	with open(FILE) as FH:
		for line in FH:
			if "#" not in line[0]:
				chrom,start, stop, INFO 	= line.strip("\n").split("\t")
				if chrom not in G:
					G[chrom] 	= list()

				if FILTER is None or (chrom in FILTER  and not FILTER[chrom].searchInterval((int(start), int(stop)) )):
					G[chrom].append((int(start), int(stop), INFO))
					ct+=1
				N+=1
	return G
def refseq_table(FILE, TSS=True):
	G 	= {}
	with open(FILE) as FH:
		header 	= True
		for line in FH:
			if not header:
				N,chrom,strand, start, stop  	= line.strip("\n").split("\t")[1:6]
				if chrom not in G:
					G[chrom] 	= list()
				if TSS:
					if strand == "+":
						G[chrom].append((int(start)-500, int(start) + 500, N))
					else:
						G[chrom].append((int(stop)-500, int(stop) + 500, N))
				else:
					G[chrom].append((int(start) , int(stop)  , N))
			else:
				header=False
	FHW 	= open("/Users/joazofeifa/Lab/genome_files/TSS.bed","w")
	for chrom in G:
		for st, sp, N in G[chrom]:
			FHW.write(chrom + "\t" + str(st) + "\t" + str(sp) + "\t" + N + "\n" )
		G[chrom] 	= node.tree(G[chrom])
	FHW.close()
	return G
def clear_eRNA_files(ROOT):
	for DIR in os.listdir(ROOT):
		for FILE in os.listdir(ROOT+DIR+"/EMG_out_files/"):
			if "eRNAs" in FILE:
				os.remove(ROOT+DIR+"/EMG_out_files/" + FILE)
					
def make_eRNA_files(ROOT,R=None):
	for DIR in os.listdir(ROOT):
		if DIR != "HCT116":
			for FILE in os.listdir(ROOT+DIR+"/EMG_out_files/"):
				if "bidirectional" in FILE and "eRNAs_" not in FILE:
					G 	= bed_file(ROOT+DIR+"/EMG_out_files/" + FILE, FILTER=R)
					OUT = ROOT+DIR+"/EMG_out_files/" + "eRNAs_" + FILE
					write_out(G,OUT)
def load_eRNA_files(ROOT, R=None):
	A 	= {}
	for DIR in os.listdir(ROOT):
		A[DIR] 	= list()
		if DIR != "HCT116":
			for FILE in os.listdir(ROOT+DIR+"/EMG_out_files/"):
				if  "eRNAs_"  in FILE:
					G 	= bed_file(ROOT+DIR+"/EMG_out_files/" + FILE, FILTER=R)
					A[DIR].append((FILE, G))
	return A

def get_counts(A,B):
	O 	= 0
	for chrom in A:
		if chrom in B and len(B[chrom]) > 0 and len(A[chrom])> 0:
			a,b 	= A[chrom],B[chrom]
			ST 	= intervals.comparison((a,b ), verbose=False)
			O+=len(ST.find_overlaps(0,1))
	return O

def make_counts(A, OUT):
	intervals 	= [(D + "-" + str(j+1),rep,G) for D in A for j,(rep, G) in enumerate(A[D])]
	intervals 	= [(i,D, rep, G) for i,(D,rep, G) in enumerate(intervals)]
	header 		= "#" + ",".join([str(i) + "_" + D for i, D, rep, G in intervals]) + "\n"
	X 			= np.zeros((len(intervals), len(intervals)))
	for i in range(len(intervals)):
		print 100*float(i) / len(intervals), "% done"
		X[i,i] 	= sum([len(intervals[i][3][chrom]) for chrom in intervals[i][3]  ])
		for j in range(i+1, len(intervals)):
			X[i,j] 	= get_counts(intervals[i][3], intervals[j][3])
			X[j,i] 	= X[i,j]
	FHW 	= open(OUT, "w")
	FHW.write(header)
	for i in range(len(intervals)):
		for j in range(len(intervals)):
			if j + 1 < len(intervals):
				FHW.write(str(X[i,j]) + ",")
			else:
				FHW.write(str(X[i,j]) + "\n")
				

def load_overlap_counts(FILE):
	D=list()
	header = True
	i 		= 0
	with open(FILE) as FH:
		for line in FH:
			if not header:
				d 	= [float(x) for x in line.strip("\n").split(",")]
				D.append(d)
				i+=1
			else:
				IDS 	= dict([(int(x.split("_")[0]),x.split("_")[1]) for x in line[1:].strip("\n").split(",")])
				
				header=False
	D 	= np.array(D)
	return D,IDS
def transform_counts_distance(D):
	for i in range(D.shape[0]):
		D[i,:]/=D[i,i]
		D[i,:]=1-D[i,:]
	return D
def add_cell_type(ID):
	cell_types 	= {}
	cell_types["Allen2014"] 	= "HCT116"
	cell_types["Andersson2014"] = "HeLa"
	cell_types["Core2014"] 		= "K562"
	cell_types["Danko2015"] 	= "T-Cells"
	cell_types["Danko2014"] 	= "K562"
	cell_types["Duttke2015"] 	= "HeLa"
	cell_types["Fong2014"] 		= "HEK293T"
	cell_types["Gardini2014"] 	= "HeLa"
	cell_types["Hah2013"] 		= "MCF7"
	cell_types["Jin2013"] 		= "IMR90"
	cell_types["Jin2014"] 		= "LNCaP"
	cell_types["Le2013"] 		= "MCF7"
	cell_types["Li2013"] 		= "MCF7"
	cell_types["Liu2013"] 		= "HEK293T"
	cell_types["Puc2015"] 		= "LNCaP"
	cell_types["Margazzi2012"] 	= "A549"
	cell_types["Meng2014a"] 	= "B-Cells"
	cell_types["Saponaro2014"] 	= "HEK293T" 
	cell_types["Sigova2013"] 	= "H1-hESC"
	cell_types["Wang2013"] 		= "B-cells"
	cell_types["Yang2014"] 		= "LNCaP"
	cell_types["Luo2014"] 		= "AC16"
	cell_types["Chen2014"] 		= "H1-hESC"
	for key in cell_types:

		if key in ID:
			return cell_types[key]+"_" +ID
	return "Not Found?" + ID
def get_colors(ID):
	cell_types 	= {}
	cell_types["Allen2014"] 	= "HCT116"
	cell_types["Andersson2014"] = "HeLa"
	cell_types["Core2014"] 		= "K562"
	cell_types["Danko2015"] 	= "T-Cells"
	cell_types["Danko2014"] 	= "K562"
	cell_types["Duttke2015"] 	= "HeLa"
	cell_types["Fong2014"] 		= "HEK293T"
	cell_types["Gardini2014"] 	= "HeLa"
	cell_types["Hah2013"] 		= "MCF7"
	cell_types["Jin2013"] 		= "IMR90"
	cell_types["Jin2014"] 		= "LNCaP"
	cell_types["Le2013"] 		= "MCF7"
	cell_types["Li2013"] 		= "MCF7"
	cell_types["Liu2013"] 		= "HEK293T"
	cell_types["Puc2015"] 		= "LNCaP"
	cell_types["Margazzi2012"] 	= "A549"
	cell_types["Meng2014a"] 	= "B-Cells"
	cell_types["Saponaro2014"] 	= "HEK293T" 
	cell_types["Sigova2013"] 	= "H1-hESC"
	cell_types["Wang2013"] 		= "B-Cells"
	cell_types["Yang2014"] 		= "LNCaP"
	cell_types["Luo2014"] 		= "AC16"
	cell_types["Core2008"] 		= "IMR90"
	
	colors 						= {}
	colors["HeLa"] 				= "red"
	colors["T-Cells"] 			= "blue"
	colors["B-Cells"] 			= "blue"
	colors["LNCaP"] 			= "green"
	colors["MCF7"] 				= "purple"
	colors["HEK293T"] 			= "teal"



	for key in cell_types:

		if key in ID:
			if cell_types[key] in colors:
				return colors[cell_types[key]]
			return "black"
	return "black"


def add_legend(ax):
	most_x, most_y 	= (300,2.15)
	corner 	= (most_x,most_y)
	width 	= 30
	height 	= 0.15
	fs 		= 14
	patch 	=   patches.Rectangle(
        (corner[0], corner[1]),   # (x,y)
        width,          # width
        height,          # height
        color="red"
    )
	rx, ry = patch.get_xy()
	cx = rx + patch.get_width()/2.0
	cy = ry + patch.get_height()/2.0
	ax.annotate("HeLa", (cx, cy), color='w', 
            fontsize=fs, ha='center', va='center')

	ax.add_patch(patch)

	corner 	= (most_x,most_y-height*1.1)
	patch 	=   patches.Rectangle(
        (corner[0], corner[1]),   # (x,y)
        width,          # width
        height,          # height
        color="blue"
    )
	rx, ry = patch.get_xy()
	cx = rx + patch.get_width()/2.0
	cy = ry + patch.get_height()/2.0
	ax.annotate("T/B cell", (cx, cy), color='w', 
            fontsize=fs, ha='center', va='center')
	
	ax.add_patch(patch)
	corner 	= (most_x+width*1.025,most_y )
	patch 	=   patches.Rectangle(
        (corner[0], corner[1]),   # (x,y)
        width,          # width
        height,          # height
        color="teal"
    )
	rx, ry = patch.get_xy()
	cx = rx + patch.get_width()/2.0
	cy = ry + patch.get_height()/2.0
	ax.annotate("HEK293T", (cx, cy), color='w', 
            fontsize=fs, ha='center', va='center')
	
	ax.add_patch(patch)

	ax.add_patch(patch)
	corner 	= (most_x+width*1.025,most_y-height*1.1)
	patch 	=   patches.Rectangle(
        (corner[0], corner[1]),   # (x,y)
        width,          # width
        height,          # height
        color="green"
    )
	rx, ry = patch.get_xy()
	cx = rx + patch.get_width()/2.0
	cy = ry + patch.get_height()/2.0
	ax.annotate("LNCaP", (cx, cy), color='w', 
            fontsize=fs, ha='center', va='center')
	
	ax.add_patch(patch)

	ax.add_patch(patch)
	corner 	= (most_x,most_y-height*2.2)
	patch 	=   patches.Rectangle(
        (corner[0], corner[1]),   # (x,y)
        width,          # width
        height,          # height
        color="purple"
    )
	rx, ry = patch.get_xy()
	cx = rx + patch.get_width()/2.0
	cy = ry + patch.get_height()/2.0
	ax.annotate("MCF7", (cx, cy), color='w', 
            fontsize=fs, ha='center', va='center')
	
	ax.add_patch(patch)

	ax.add_patch(patch)
	corner 	= (most_x+width*1.025,most_y-height*2.2)
	patch 	=   patches.Rectangle(
        (corner[0], corner[1]),   # (x,y)
        width,          # width
        height,          # height
        color="black"
    )
	rx, ry = patch.get_xy()
	cx = rx + patch.get_width()/2.0
	cy = ry + patch.get_height()/2.0
	ax.annotate("Other", (cx, cy), color='w', 
            fontsize=fs, ha='center', va='center')
	
	ax.add_patch(patch)

def annotate_branches(ddata, ax):
	for k,(i, d) in enumerate(zip(ddata['icoord'], ddata['dcoord'])):
		x = 0.5 * sum(i[1:3])
		y = d[1]
		name = ""		
		xytext=None
		if k == 9:
			name 	= "MCF7/\nHEK293T"
			xytext 	= (-90,90)
		if k == 30:
			name 	= "T/B cell"
			xytext 	= (40,90)
		if k == 25:
			name 	= "HeLa"
			xytext 	= (-80,30)
		if k == 14:
			name 	= "LNCaP"
			xytext 	= (90,70)
		
			

		if name:

			ax.plot(x, y, 'ro' , alpha=0.5)
			ax.annotate(name , (x, y), xytext=xytext,
				textcoords='offset points',
				va='top', ha='center', arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"), fontsize=18)	



def cluster(D,IDS):
	F 	= plt.figure(figsize=(15,10))
	ax 	= F.add_subplot(1,1,1)
	OIDS 	= [IDS[i] for i in range(len(IDS))]
	IDS 	= [add_cell_type(IDS[i]) for i in range(len(IDS))]
	colors  = [get_colors(IDS[i]) for i in range(len(IDS))]
	D 	= transform_counts_distance(D)
	linkage_matrix = linkage(D, "complete")
			
	ddata = dendrogram(linkage_matrix,
 	                  labels=OIDS , color_threshold=0,leaf_rotation= 85, leaf_font_size=14 )
	annotate_branches(ddata, ax)
	
	ax.set_ylim(0, 2.345)
	ylbls = ax.get_xmajorticklabels()
	ticks = ax.get_xticks()
	for i,lbl in enumerate(ylbls):
		c 	= get_colors(add_cell_type(lbl.get_text() ))
		lbl.set_color(c)
		if c!= "black":
			lbl.set_weight("bold")
		lbl.set_text(OIDS[i])
	add_legend(ax)
	ax.set_yticklabels([])
	plt.tight_layout()
	plt.show()

if __name__ == "__main__":
	ROOT 			= "/Users/joazofeifa/Lab/gro_seq_files/"
	R 				= refseq_table("/Users/joazofeifa/Lab/genome_files/RefSeqHG19.txt")
	DISTANCE_FILE 	= os.getcwd() + "/Distance.txt"
	MAKE_ERNAS 		= False
	MAKE_COUNTS  	= False
	CLUSTER 		= True

	if MAKE_ERNAS:
		clear_eRNA_files(ROOT)
		make_eRNA_files(ROOT,R=R)
	if MAKE_COUNTS:
		A 			= load_eRNA_files(ROOT)
		make_counts(A,DISTANCE_FILE)
	if CLUSTER:	
		D,IDS 			= load_overlap_counts(DISTANCE_FILE)
		cluster(D,IDS)


















						
