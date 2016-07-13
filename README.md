# Tfit
Transcription fit (Tfit) implements a finite mixture model to identify sites of bidirectional or divergent transcription in nascent trancsription assays such as Global Run-On (GRO) and Precision Run-on (PRO) followed by sequencing data. Tfit is separated by two modules: (1) bidir and (2) model. Output from these modules can be exported to any genome browser, bed file examples are provided below.  

![Alt text](https://github.com/azofeifa/Tfit/blob/master/images/Example_Snapshot.png)

These are invoked as below.

$Tfit bidir \<list of parameters and flags\>

$Tfit model \<list of parameters and flags\>

The bidir module will compute local likelihood statistics given a template mixture model (estimated either from promoter associated transcription) or specified explicitly by the user (the former is encoraged). This method is fast and will finish in about 30 minutes on a single node single CPU machine. The output will be a bed file corresponding to areas of possible bidirectional transcription. This output is dicussed heavily in later sections

The model module will compute full MLE estimates of the mixture model at user specified regions of the genome provided as a bed file. This bed file may be the output from the bidir module. It is recomended that users fit MLE estimates to the output of the bidir module as this will greatly decrease false positives. Two files will output from this module: (A) a \<job_name\>_K_models_MLE.tsv and a \<job_name\>_divergent_classifications.bed both containing information regarding the location, spreading, pausing probability and strand probability of bidirectional events. This output is dicussed heavily in later sections. Computation time of this module depends on model complexity bounds. If the user is profiling only for bidirectional transcripts, it is recommended that -mink and -maxk flags be set to 1. In this case, computation take will take around 2-3 hours on a single node, single CPU machine.  

Please note that for both the bidir and model module computation time will decrease linearly with the number of available CPUs. Computation time is discussed throughout the below sections.  

##System Requirements and Makefile
Transcription Fit (Tfit) is written in the C/C++ program language that requires GNU compilers 4.7.3 or greater. Tfit uses the popular openMPI framework to perform massive parallelization via multithreading on multiple core, single node systems or multiple core, multiple node compute clusters. 

After cloning this repo, please change directory into /where/you/clone/this/repo/Tfit/src/ and run make.  

$cd  /where/you/clone/this/repo/Tfit/src/

$make clean

$make

If the program compilled sucessfuly you should see the below output.

-------------------

GCC version: 5.1.0

...


successfully compiled

-------------------


If your program, did not compile properly it is likely that you do not have the correct dependencies. The three most significant dependencies are listed below. 

1) c++11 (this ships with most recent versions of GCC, please visit https://gcc.gnu.org/install/)

2) openmp (this ships with most recent versions of GCC, please visit https://gcc.gnu.org/install/)

3) MPI (this needs to installed and configured and serves as a wrapper for GCC, please visit https://www.open-mpi.org/faq/)


















##File Formats

![Alt text](https://github.com/azofeifa/Tfit/blob/master/images/bedgraph_joint_example.png)

![Alt text](https://github.com/azofeifa/Tfit/blob/master/images/bed_file_example.png)

![Alt text](https://github.com/azofeifa/Tfit/blob/master/images/config_file_example.png)

##Bidir Module
Bidir module utilizes a penalized likelihood ratio (LLR) test to scan across the genome for areas resembling bidirectional transcription. Neighoring genomic coordinates where the LLR exceeds some user defined threshold (-bct flag) are joined and are returned as a bed file (chrom[tab]start[tab]stop[newline]). 

Perhaps the most important user input file is the "BedGraph" File corresponding to the forward and reverse strad. This file is simple (chrom[tab]start[tab]stop[tab]coverage[newline]). Importantly, start and stop should be integer valued. Although coveraged can be a float, we recommend this to be an integer as well, i.e. normalization by millions mapped for example made lead to numerical degeneracies in the optimization routine. An example bedgraph file is listed below. 

![Alt text](https://github.com/azofeifa/Tfit/blob/master/images/bedgraph_single_example.png)


The critical input parameters are listed below:

1. -i	= \</path/to/BedGraphFile_forward strand> “BedGraph File from above”
2. -j = \</path/to/BedGraphFile_forward strand> “BedGraph File from above”
3. -N = job_name, simple a string
4. -o = \</path/to/output/directory>
5. -log_out = \</path/to/logoutput/directory> "As the program runs this file will be updated with progress 

The non-critical input parameters are listed below, these all have default settings.
1. -tss = \</path/to/bedfile/of/promoter/locations/ (promoter locations are provided for hg19 and mm10 in the annotations/ directory of this repo, it is recommended to optimize your template density function by promoter or TSS associated regions 
2. -chr = a [string] where the bidir module will only run on specified chromosome (default is "all")
3. -bct = a [floating point value], this is the LLR threshold, the default and recommended is 1
^^^If TSS is not provided, the user can manually enter the template parameters of interest^^^ 
4. -lambda = a [floating point value], this is the entry length parameter for the EMG density function (default = 200 bp)   
5. -sigma  = a [floating point value], this is the variance parameter for the EMG density function (default = 10 bp)
6. -pi     = a [floating point value], this is the strand bias parameter for the EMG density function (default = 0.5)
7. -w      = a [floating point value], this is the pausing probability parameter for the EMG density function (default = 0.5)
8. -foot_print = a [floating point value], this is the foot_print parameter for the EMG density function (default = 50)





##Model Module and List of Parameters

##References




