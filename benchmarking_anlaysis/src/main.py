import load, draw
def main():
	BM 			= False
	RI 			= True
	if BM:


		pando 		= "/Users/joeyazo/Desktop/Lab/benchmarking_files/pando/"
		vieques 	= "/Users/joeyazo/Desktop/Lab/benchmarking_files/vieques/"
		
		A 			= load.all_files(vieques)
		G 			= load.all_files(pando)
		draw.wall_vs_np(G,A)
	if RI:
		directory 	= "/Users/joeyazo/Desktop/Lab/EMG_files/RI_files/"
		G 			= load.RI_directory(directory)
		draw.delta_ll_vs_RI(G)
if __name__ == "__main__":
	main()