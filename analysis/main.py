import load
import look_at_parameters as lap
import numpy as np
import BIC
import display_model_fits as dmf
import correlations
import density_plots
def run(root):

	display_fits 	= False
	parameters 		= False
	correlation 	= False
	correlation_BO 	= True
	if correlation_BO:
		DIR 			="/Users/joazofeifa/Lab/gro_seq_files/HCT116/EMG_out_files/"
		DMSO1hr101911 	="DMSO1hr101911_model_fits/EMG-4_bidirectional_hits_intervals.bed"
		DMSO1027 		="DMSO1027_1212_model_fits/EMG-3_bidirectional_hits_intervals.bed"
		Ma6_NoIndex 	="Ma6_NoIndex_L008_R1_001/EMG-6_bidirectional_hits_intervals.bed"
		DMSO2_3 		="DMSO2_3_model_fits/EMG-1_bidirectional_hits_intervals.bed"
		Nutlin2_3 		= "Nutlin2_3_model_fits/EMG-2_bidirectional_hits_intervals.bed"
		
		RefSeq 			= "/Users/joazofeifa/Lab/genome_files/RefSeqHG19.txt"
		ChIP_p53 		= "/Users/joazofeifa/Lab/ACM_IEEE_Paper_analysis/files/bedFiles/Atleast7of7.bedbothstrands.bed_norefgene.bed"
		ChIP_p53 		= "/Users/joazofeifa/Lab/nutlin_bidirectional_hits_intervals_091715.bed.count.bed.h.bed.namescoreDMSObi.resSig.txt.bed.txt"
		DMSO_forward 	= "/Users/joazofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/DMSO2_3.pos.BedGraph"
		DMSO_reverse 	= "/Users/joazofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/DMSO2_3.neg.BedGraph"
		Nutlin_forward 	= "/Users/joazofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/Nutlin2_3.sorted.pos.BedGraph"
		Nutlin_reverse 	= "/Users/joazofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/Nutlin2_3.sorted.neg.BedGraph"

		DMSO1hr101911_L,DMSO1hr101911_G = load.load_model_fits_bed_file(DIR+DMSO1hr101911)
		DMSO1027_L,DMSO1027_G 			= load.load_model_fits_bed_file(DIR+DMSO1027)
		Ma6_NoIndex_L,Ma6_NoIndex_G 	= load.load_model_fits_bed_file(DIR+Ma6_NoIndex)
		DMSO2_3_L,DMSO2_3_G 			= load.load_model_fits_bed_file(DIR+DMSO2_3)
		Nutlin2_3_L,Nutlin2_3_G 		= load.load_model_fits_bed_file(DIR+Nutlin2_3)
		R 								= load.load_Refseq(RefSeq)
		ChIP 							= load.load_ChIP_p53(ChIP_p53)
		load.label(DMSO2_3_L, R, "annotated")
		load.label(DMSO2_3_L, ChIP, "p53_site")
		load.label(Nutlin2_3_L, ChIP, "p53_site")
		load.label(Nutlin2_3_L, R, "annotated")
		# density_plots.insert_bedgraph(DMSO2_3_L,(DMSO_forward,DMSO_reverse ))
		# density_plots.insert_bedgraph(Nutlin2_3_L,(Nutlin_forward,Nutlin_reverse ))




		overlaps 						= correlations.match_UP(Ma6_NoIndex_L, DMSO1hr101911_L)
#		density_plots.plot_density(overlaps)
#		correlations.p53_binding(Nutlin2_3_L, DMSO2_3_L, overlaps)
#		correlations.label_p53(overlaps, attr="lam", LOG=True )
#		correlations.promoter_differences_test((DMSO2_3_L,Nutlin2_3_L))		
#		correlations.p53_differences_test((Nutlin2_3_L,))		

#		correlations.si_lam(overlaps)
		correlations.run(overlaps, attr="si", LOG=False	 )
	if correlation:
		DIR 			="/Users/joazofeifa/Lab/gro_seq_files/HCT116/EMG_out_files/"
		DMSO1hr101911 	="DMSO1hr101911_model_fits/model_fits.txt"
		DMSO1027 		="DMSO1027_1212_model_fits/model_fits.txt"
		Ma6_NoIndex 	="Ma6_NoIndex_L008_R1_001/model_fits.txt"
		DMSO2_3 		="DMSO2_3_model_fits/model_fits.txt"
		Nutlin2_3 		= "Nutlin2_3_model_fits/model_fits.txt"
		DMSO1hr101911_L,DMSO1hr101911_G = load.load_model_fits_bed_file(DIR+DMSO1hr101911)
		DMSO1027_L,DMSO1027_G 			= load.load_model_fits_bed_file(DIR+DMSO1027)
		Ma6_NoIndex_L,Ma6_NoIndex_L 	= load.load_model_fits_bed_file(DIR+Ma6_NoIndex)
		DMSO2_3_L,DMSO2_3_G 			= load.load_model_fits_bed_file(DIR+DMSO2_3)
		Nutlin2_3_L,Nutlin2_3_G 		= load.load_model_fits_bed_file(DIR+Nutlin2_3)
		
		correlations.run(DMSO2_3_L,DMSO2_3_L,Ma6_NoIndex_L,Ma6_NoIndex_L )

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
