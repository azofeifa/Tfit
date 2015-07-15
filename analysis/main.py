import load
import look_at_parameters as lap
import numpy as np
import BIC
import display_model_fits as dmf
def run(root):

	display_fits 	= True
	parameters 		= False
	if display_fits:
		out_dir 	= "/Users/joeyazo/Desktop/Lab/gro_seq_files/HCT116/EMG_out_files/"
		model_file 	= out_dir+"model_fits_out_all_4"
		data_file 	= out_dir+"test_file_2.tsv"

		intervals 	= load.EMG_out(model_file)
		load.insert_data(data_file, intervals)
		dmf.display(intervals,bins=300)

	if parameters:


		EMG_out_FILE 	= root + "gro_seq_files/HCT116/EMG_out_files/EMG_model_fits_all_0"
		parameters 		= False
		BIC_analysis 	= True 
		#only supports loading one at a time
		fits 			= load.EMG_out(EMG_out_FILE)
		if parameters:
			lap.run(fits, spec=None, 
				weight_thresh=0.1,retry_tresh=0,
				converged=True)
		if BIC_analysis:
			BIC.run(fits)

if __name__=="__main__":
	root 	= "/Users/joeyazo/Desktop/Lab/" 
	run(root)
