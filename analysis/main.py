import load
import look_at_parameters as lap
import numpy as np
import BIC
import display_model_fits as dmf
import correlations
def run(root):

	display_fits 	= False
	parameters 		= False
	correlation 	= True
	if correlation:
		DIR 			="/Users/joazofeifa/Lab/gro_seq_files/HCT116/EMG_out_files/"
		DMSO1hr101911 	="DMSO1hr101911_model_fits/model_fits.txt"
		DMSO1027 		="DMSO1027_1212_model_fits/model_fits.txt"
		Ma6_NoIndex 	="Ma6_NoIndex_L008_R1_001/model_fits.txt"
		DMSO2_3 		="DMSO2_3_model_fits/model_fits.txt"
		Nutlin2_3 		= "Nutlin2_3_model_fits/model_fits.txt"
		DMSO1hr101911_L,DMSO1hr101911_G = load.load_model_fits_txt_file(DIR+DMSO1hr101911)
		DMSO1027_L,DMSO1027_G 			= load.load_model_fits_txt_file(DIR+DMSO1027)
		Ma6_NoIndex_L,Ma6_NoIndex_G 	= load.load_model_fits_txt_file(DIR+Ma6_NoIndex)
		DMSO2_3_L,DMSO2_3_G 			= load.load_model_fits_txt_file(DIR+DMSO2_3)
		Nutlin2_3_L,Nutlin2_3_G 		= load.load_model_fits_txt_file(DIR+Nutlin2_3)
		
		correlations.run(DMSO2_3_L,DMSO2_3_L,DMSO1027_L,Nutlin2_3_G )

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
