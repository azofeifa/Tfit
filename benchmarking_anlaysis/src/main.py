import load, draw
def main():
	pando 		= "/Users/joeyazo/Desktop/Lab/benchmarking_files/pando/"
	vieques 	= "/Users/joeyazo/Desktop/Lab/benchmarking_files/vieques/"
	
	A 			= load.all_files(vieques)
	G 			= load.all_files(pando)
	draw.wall_vs_np(G,A)
	pass


if __name__ == "__main__":
	main()