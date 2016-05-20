# Tfit
Transcription fit (Tfit) implements a mixture model to identify sites of bidirectional transcription in Global Run-On followed by sequencing data. Tfit is separate by two modules: (1) bidir and (2) model. These are invoked as below. 

$Tfit bidir \<list of parameters and flags\>

$Tfit model \<list of parametres and flags\>

The bidir module will compute local likelihood statistics given a template mixture model (trained either from promoter associated transcription) or specified explicitly by the user (the former is encoraged). This method is fast and will finish in about 10 minutes on a single node single CPU machine. The output will be a bed file corresponding to areas of possible bidirectional transcription. This output is dicussed heavily in later sections

The model module will compute full MLE estimates of the mixture model at user specified regions of the genome provided as a bed file. This bed file may be the output from the bidir module. It is recomended that users fit MLE estimates to the output of the bidir module as this will greatly decrease false positives. Two files will output from this module: (A) a \<job_name\>_K_models_MLE.tsv and a \<job_name\>_divergent_classifications.bed both containing information regarding the location, spreading, pausing probability and strand probability of bidirectional events. This output is dicussed heavily in later sections. Computation time of this module depends on model complexity bounds. If the user is profiling only for bidirectional transcripts, it is recommended that -minK and -maxK flags be set to 1. In this case, computation take will take around 2-3 hours on a single node, single CPU machine.  

Please note that for both the bidir and model module computation time will decrease linearly with the number of available CPUs. Computation time is discussed throughout the below sections.  

##System Requirements and Makefile

##File Formats

##Bidir Module and List of Parameters 

##Model Module and List of Parameters

##References



