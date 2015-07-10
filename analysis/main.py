import load
import look_at_parameters as lap
import numpy as np
import BIC
def run(root):
	EMG_out_FILE 	= root + "gro_seq_files/HCT116/EMG_out_files/EMG_model_fits_all_0"
	parameters 		= True
	BIC_analysis 	= False
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
