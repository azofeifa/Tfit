import os
import re, numpy as np, math, matplotlib.pyplot as plt
def convert(FILE, results_dir, out_dir):
	header 	= "track name=Strand_-description=\"FStitch\" visibility=2 useScore=2 cgGrades=50 cgColour1=white cgColour2=yellow cgColour3=red height=30\n"
	FHW 	= open(out_dir+FILE+".bed", "w")
	FHW.write(header)
	with open(results_dir+FILE) as FH:
		prevChrom, prevState 	= "",  ""
		st, sp 					= None, None
		for line in FH:
			chrom, start, stop, state, cov1,cov2 	= re.split("\s+", line.strip("\n"))
			if chrom!=prevChrom or state != prevState:
				if st is not None and sp is not None:
					if prevState=="0":
						RGB, INFO 	= "255,0,0", "100"
					elif prevState=="1":
						RGB, INFO 	= "0,255,0", "200"
					elif prevState=="2":
						RGB, INFO 	= "0,0,255", "300"

						

					FHW.write(prevChrom+"\t" + st + "\t" + sp + "\t" + prevState + "\t"  + INFO + "\t-\t" + st + "\t"
					+ sp + "\t" + RGB + "\t.\n" )

				st 			= start
				prevState 	= state
			sp 	= stop
			prevChrom=chrom
	FHW.close()
def iterate_files(results_dir, out_dir):
	for FILE in os.listdir(results_dir):
		if os.path.isfile(results_dir+FILE):
			print FILE
			convert(FILE, results_dir, out_dir)
def graph_var():

	scalar 	= 0.1
	f 	= lambda x : 0.933*math.exp(-x*0.633)+ 0.403
	for scalar in np.linspace(0.09, 0.1, 10):
		xs 	= np.linspace(0,1000, 1000)
		xs+=1.0
		xs 	= np.log(xs)
		ys 	= map(f, xs)
		plt.plot(xs, ys, label=str(scalar))
	plt.legend()
	plt.show()






if __name__ == "__main__":
	print os.getcwd()
	results_dir 	= "/Users/joeyazo/Desktop/Lab/SMART-Project/results/"
	out_dir 		= "/Users/joeyazo/Desktop/Lab/SMART-Project/results/igv_results/"
	#iterate_files(results_dir, out_dir)
	graph_var()