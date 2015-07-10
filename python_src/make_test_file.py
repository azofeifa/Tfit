import load, simulate
def run(FILE):
	FHW 	= open(FILE, "w")
	X 		= load.grab_specific_region("chr1",6229860,6303055, SHOW=False, bins=1000 )
	print min(X[:,0]), max(X[:,0])
	FHW.write("#chr1,6229860,6303055\n")
	FHW.write("~forward\n")
	for i in range(X.shape[0]):
		FHW.write(str(X[i,0]) + "," + str(X[i,1]) + "\n")
	FHW.write("~reverse\n")
	for i in range(X.shape[0]):
		FHW.write(str(X[i,0]) + "," + str(X[i,2]) + "\n")
	X 	= simulate.runOne(mu=0, s=0.1, l=3, lr=100, ll=-50, we=0.5,wl=0.25, wr=0.25, pie=0.5, pil=0.1, pir=0.9, 
		N=1000, SHOW=False, bins=1000, noise=True )
	X[:,0]+=6303055
	X[:,0]*=100.
	st, sp 	= 	X[0,0], X[-1,0]
	print st, sp
	FHW.write("#chrN,"+str(st) + "," + str(sp) + "\n")
	FHW.write("~forward\n")
	for i in range(X.shape[0]):
		FHW.write(str(X[i,0]) + "," + str(X[i,1]) + "\n")
	FHW.write("~reverse\n")
	for i in range(X.shape[0]):
		FHW.write(str(X[i,0]) + "," + str(X[i,2]) + "\n")
	FHW.close()
	
	
	
	
	



if __name__=="__main__":
	FILE 	= "/Users/joeyazo/Desktop/Lab/gro_seq_files/HCT116/EMG_out_files/test_file_2.tsv"
	run(FILE)