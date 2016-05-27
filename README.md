# Tfit
Transcription fit (Tfit) implements a mixture model to identify sites of bidirectional transcription in nascent trancsription assays such as Global Run-On (GRO) and Precision Run-on (PRO) followed by sequencing data. Tfit is separated by two modules: (1) bidir and (2) model. Example output from these modules can be exported to any genome browser, bed file examples are provided below.  

![Alt text](https://github.com/azofeifa/Tfit/blob/master/images/Example_Snapshot.png)

These are invoked as below.

$Tfit bidir \<list of parameters and flags\>

$Tfit model \<list of parametres and flags\>

The bidir module will compute local likelihood statistics given a template mixture model (estimated either from promoter associated transcription) or specified explicitly by the user (the former is encoraged). This method is fast and will finish in about 30 minutes on a single node single CPU machine. The output will be a bed file corresponding to areas of possible bidirectional transcription. This output is dicussed heavily in later sections

The model module will compute full MLE estimates of the mixture model at user specified regions of the genome provided as a bed file. This bed file may be the output from the bidir module. It is recomended that users fit MLE estimates to the output of the bidir module as this will greatly decrease false positives. Two files will output from this module: (A) a \<job_name\>_K_models_MLE.tsv and a \<job_name\>_divergent_classifications.bed both containing information regarding the location, spreading, pausing probability and strand probability of bidirectional events. This output is dicussed heavily in later sections. Computation time of this module depends on model complexity bounds. If the user is profiling only for bidirectional transcripts, it is recommended that -minK and -maxK flags be set to 1. In this case, computation take will take around 2-3 hours on a single node, single CPU machine.  

Please note that for both the bidir and model module computation time will decrease linearly with the number of available CPUs. Computation time is discussed throughout the below sections.  

##System Requirements and Makefile
Transcription Fit (Tfit) is written in the C/C++ program language that requires GNU compilers 4.7.3 or greater. Tfit uses the popular openMPI framework to perform massive parallelization via multithreading on multiple core, single node systems or multiple core, multiple node compute clusters. 

After cloning this repo, please change directory into /where/you/clone/this/repo/Tfit/src/ and run make.  

$cd  /where/you/clone/this/repo/Tfit/src/

$make clean

$make

If the program compilled sucessfuly you should see the below output.

\=========================================

GCC version: 5.1.0

...

\=
successfully compiled


If your program, did not compile properly it is likely that you do not have the correct dependencies. The three most significant dependencies are listed below. 

1) c++11 (this ships with most recent versions of GCC, please visit https://gcc.gnu.org/install/)

2) openmp (this ships with most recent versions of GCC, please visit https://gcc.gnu.org/install/)

3) MPI (this needs to installed and configured and serves as a wrapper for GCC, please visit https://www.open-mpi.org/faq/)


















##File Formats

![Alt text](https://github.com/azofeifa/Tfit/blob/master/images/bedgraph_single_example.png)

![Alt text](https://github.com/azofeifa/Tfit/blob/master/images/bedgraph_joint_example.png)

![Alt text](https://github.com/azofeifa/Tfit/blob/master/images/bed_file_example.png)

![Alt text](https://github.com/azofeifa/Tfit/blob/master/images/config_file_example.png)

##Bidir Module and List of Parameters 

##Model Module and List of Parameters

##References




