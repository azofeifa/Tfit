import sys, BIC, load
def check_bidir_component(rv):
	si, l, w 	= rv.si, rv.l, rv.w
	pass

def run(G, penality, diff_threshold):
	for I in G:
		model 	= BIC.get_best_model(I, penality , diff_threshold)
		bidirs 	= [for rv in rvs if check_bidir_component(rv)]







if __name__ == "__main__":
	if len(sys.argv)==1:
		merge_data_out 		= "/Users/joeyazo/Desktop/Lab/gro_seq_files/HCT116/merged_data_file_20.txt"
		G 					= load.merge_data_out(merge_data_out)
		penality 			= 100
		diff_threshold 		= 5
		run(G, penality, diff_threshold)
