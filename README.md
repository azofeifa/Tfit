# Tfit
Transcription fit (Tfit) implements a finite mixture model to identify sites of bidirectional or divergent transcription in nascent transcription assays such as Global Run-On (GRO) and Precision Run-on (PRO) followed by sequencing data. Tfit is separated by two modules: (1) bidir and (2) model. Output from these modules can be exported to any genome browser, bed file examples are provided below.  

![Alt text](https://github.com/azofeifa/Tfit/blob/master/images/Example_Snapshot.png)

The bidir module will compute local likelihood statistics given a fixed template mixture model (estimated either from promoter associated transcription) or specified explicitly by the user (the former is encouraged). This method is fast and will finish in about 30 minutes on a single node single CPU machine. The output will be a bed file corresponding to areas of possible bidirectional transcription. This output is discussed heavily in later sections

The model module will compute full MLE estimates of the mixture model at user specified regions of the genome provided as a bed file. This bed file may be the output from the bidir module. It is recommended that users fit MLE estimates to the output of the bidir module as this will decrease false positives. Two files will output from this module: (A) a \<job_name\>_K_models_MLE.tsv and a \<job_name\>_divergent_classifications.bed both containing information regarding the location, spreading, pausing probability and strand probability of bidirectional events. This output is discussed heavily in later sections. Computation time of this module depends on model complexity bounds. If the user is profiling only for bidirectional transcripts, it is recommended that -mink and -maxk flags be set to 1. In this case, computation take will take around 2-3 hours on a single node, single CPU machine.  

Please note that for both the bidir and model module computation time will decrease linearly with the number of available CPUs. Computation time is discussed throughout the below sections.  

##System Requirements and Makefile
Transcription Fit (Tfit) is written in the C/C++ program language that requires GNU compilers 4.7.3 or greater. Tfit uses the popular openMPI framework to perform massive parallelization via multithreading on multiple core, single node systems or multiple core, multiple node compute clusters. 

After cloning this repo, please change directory into /where/you/clone/this/repo/Tfit/src/ and run make.  

$cd  /where/you/clone/this/repo/Tfit/src/

$make clean

$make

If the program compiled successfully you should see the below output.

-------------------

GCC version: 5.1.0

...


successfully compiled

-------------------


If your program, did not compile properly it is likely that you do not have the correct dependencies. The three significant dependencies are listed below. 

1) c++11 (this ships with most recent versions of GCC, please visit https://gcc.gnu.org/install/)

2) openmp (this ships with most recent versions of GCC, please visit https://gcc.gnu.org/install/)

3) MPI (this needs to installed and configured and serves as a wrapper for GCC, please visit https://www.open-mpi.org/faq/)

##Running and Installing Through Docker Container
Running Tfit via a Docker container requires that Docker is installed and running. In Windows and Mac OSX you must run the Docker/Tfit script from the Docker Quickstart Terminal (https://www.docker.com/products/docker-toolbox).

To run Tfit via docker, simply run the Docker/Tfit script with the same arguments and options as you would use when running the standard Tfit program. This script will check if you have the biofrontiers/fstitch_tfit:latest image and download it if you do not. It will then pass the arguments you provide to the Tfit program inside the container and output the results and logs to the location you specify.

In order to ensure you are using the latest version of the container, you may run the following command before executing the Tfit script:

docker pull biofrontiers/fstitch_tfit:latest

Advanced users may want to build their own Docker images using the provided Dockerfile. To do so, run the following command from within the Docker directory:

docker build -t [image name] .

And alter the Docker/Tfit script, changing the REPOSITORY and TAG variables so that they match your built image.



##Bidir Module
The bidir module scans across the genome for areas resembling bidirectional transcription by comparing a fixed template mixture model (user provided parameters or parameters estimated from promoter regions) to a noise model (uniform distribution) by a Likelihood ratio score (LLR). In brief, the template mixture model is parameterized by -lambda (entry length or amount of skew), -sigma (variance in loading, error), -pi (strand bias, probability of forward strand data point) and -w (pausing probability, how much bidirectional signal to elongation/noise signal). Neighboring genomic coordinates where the LLR exceeds some user defined threshold (-bct flag) are joined and are returned as a bed file (chrom[tab]start[tab]stop[newline]). An example of a bed file is provided below:

![Alt text](https://github.com/azofeifa/Tfit/blob/master/images/bed_file_example.png)


Perhaps the most important user input file is the "BedGraph" File corresponding to the forward and reverse strand. This file is simple (chrom[tab]start[tab]stop[tab]coverage[newline]). Importantly, start and stop should be integer valued. Although coverage can be a float, we recommend this to be an integer as well, i.e. normalization by millions mapped for example may lead to numerical degeneracies in the optimization routine. An example bedgraph file is listed below. 

![Alt text](https://github.com/azofeifa/Tfit/blob/master/images/bedgraph_joint_example.png)


The critical input parameters are listed below:

1. -i	= \</path/to/BedGraphFile_forward strand> “BedGraph File from above”
2. -j = \</path/to/BedGraphFile_forward strand> “BedGraph File from above”
3. -ij= \</path/to/BedGraphFile_joint_strand> (if -i and -j are not specified) “BedGraph File from above but reverse strand reads are specified by negative coverage values and forward strand reads are specified by positive coverage values, NOTE: either -ij is specified or both -i and -j are specified. An example of this combined joint bedgraph file is below. 
4. -N = job_name, simple a string
5. -o = \</path/to/output/directory>
6. -log_out = \</path/to/logoutput/directory> "As the program runs this file will be updated with progress 

Example of joint forward and reverse strand bedgraph file:

![Alt text](https://github.com/azofeifa/Tfit/blob/master/images/bedgraph_single_example.png)


The non-critical input parameters are listed below, these all have default settings.

1. -tss = \</path/to/bedfile/of/promoter/locations/ (promoter locations are provided for hg19 and mm10 in the annotations/ directory of this repo, it is recommended to optimize your template density function by promoter or TSS associated regions 
2. -chr = a [string] where the bidir module will only run on specified chromosome (default is "all")
3. -bct = a [floating point value], this is the LLR threshold, the default and recommended is 1
^^^If TSS is not provided, the user can manually enter the template parameters of interest^^^ 
4. -lambda = a [floating point value], this is the entry length parameter for the EMG density function (default = 200 bp)   
5. -sigma  = a [floating point value], this is the variance parameter for the EMG density function (default = 10 bp)
6. -pi     = a [floating point value], this is the strand bias parameter for the EMG density function (default = 0.5)
7. -w      = a [floating point value], this is the pausing probability parameter for the EMG density function (default = 0.5)

After the bidir model finishes a bed file will appear in the user specified output directory called: [-N]_prelim_bidir_hits.bed. This file will contain genomic intervals of interest where divergent transcription is likely occurring. This file may be used for downstream analysis or taken at face value. Last, but not least, the bidir module can be invoked as below:

$Tfit bidir \<list of parameters and flags\>

Example output from the bidir module,i.e. [-N]_prelim_bidir_hits.bed, is provided below. 

![Alt text](https://github.com/azofeifa/Tfit/blob/master/images/Prelim_Bidir_Example.png)


##Model Module and List of Parameters
Unlike the "bidir" module which utilizes an average or template version of the mixture model to scan the entire genome quickly, the "model" module will attempt to find (by maximum likelihood estimation, MLE) the best set of parameters (sigma,lambda, pi, w) on a per region of interest basis. Such a routine is especially valuable if it is believed that pausing or strand bias is changing following experimental perturbation. In addition, running the model module on the prelim_bidir_hits.bed file will greatly decrease the number of false positives as the MLE estimates will more accurately reflect the underlying structure of the region of interest rather than a static template model. 

In short, MLE estimates are computed by the EM algorithm which is a convergent optimization method found commonly in gaussian mixture modeling and cluster analysis. Given this, the user may specify sets of parameters that are specific to the EM routine such as number of random seeds (-rounds), maximum number of iterations (-mi) and convergence threshold (-ct). A list of critical and non-critical parameters are given below.   

The critical input parameters are listed below:

1. -i	= \</path/to/BedGraphFile_forward strand> “BedGraph File from above”
2. -j = \</path/to/BedGraphFile_forward strand> “BedGraph File from above”
3. -ij= \</path/to/BedGraphFile_joint_strand> (if -i and -j are not specified) “BedGraph File from above but reverse strand reads are specified by negative coverage values and forward strand reads are specified by positive coverage values, NOTE: either -ij is specified or both -i and -j are specified. An example of this combined joint bedgraph file is below. 
4. -k = \</path/to/bed_file of regions of interest> "bed" file from above, this output may be from _prelim_bidir_hits.bed or other genomic intervals of interest, TF peaks from MACS, RefSeq annotations for genes etc. 
5. -N = job_name, simple a string
6. -o = \</path/to/output/directory>
7. -log_out = \</path/to/logoutput/directory> "As the program runs this file will be updated with progress 

The non-critical input parameters are listed below, these all have default settings.

1. -mink = a [integer value], minimum number of finite mixtures to consider (default = 1)
2. -maxk = a [integer value], maximum number of finite mixtures to consider (default = 1)
3. -rounds = a [integer value], number of random seeds to the EM (default = 5)
4. -ct = a [floating point value], convergence threshold where EM halts (default = 0.0001)
5. -mi = a [integer value],maximum number of EM iterations after which EM will halt (default = 2000)

After the model module has finished, Tfit will output two files in the user specified output directory: [-N]_K_models_MLE.tsv and [-N]_divergent_classifications.bed. 

The first file, [-N]_K_models_MLE.tsv, gives a detailed account for each interval of interest from -k input parameters, a list for each finite mixture model fits -mink to -maxk, the log-likelihood score, and MLE estimates for the center of the bidirectional transcript (mu), variance in mu (sigma), entry length (lambda), strand bias (pi), distance between bidirectional peaks (foot print) and all the associated mixture weights (basically w). In addition, stats on the number of reads over that interval etc. This can be used for manual Bayesian Information Criterion (BIC) calculations. An example output of this file where -mink = 1 and -maxk = 10 is given below:

![Alt text](https://github.com/azofeifa/Tfit/blob/master/images/K_models_example.png)

The second file,[-N]_divergent_classifications.bed, provides a new bed file, where the center of each bed interval corresponds to the center of the bidirectional peak and the width corresponds the estimated standard deviation around that estimate (essentially sigma + lambda) following BIC model comparison. This again can be used for downstream analysis and comprises the most accurate set of bidirectional prediction centers that Tfit can currently offer. An example output of this file is given below.

![Alt text](https://github.com/azofeifa/Tfit/blob/master/images/Divergent_Example.png)

##Chaining the bidir and model module
Discussed thus far, profiling for divergent or bidirectional transcription events may be achieved by first running the bidir module and then using the output (_prelim_bidir_hits.bed) as input to the model module. For convenience, these modules can be chained and all three files (_prelim_bidir_hits.bed, _divergent_classifications.bed, _K_models_MLE.tsv) will output from one single call. This is invoked like below:

$Tfit bidir -MLE 1 \<list of other parameters and flags\>

##The config file
At this point, we have discussed all the necessary parameters to run Tfit. However, specifying each of these on the command line is tedious. To this end, the user may specify a config file. This is invoked like below.

$Tfit bidir -config \</path/to/config/file.txt

An example of the config file is below. 

![Alt text](https://github.com/azofeifa/Tfit/blob/master/images/config_file_example.png)

The structure of the config file should remain this way (i.e. "-flag = value"). Statements following a # sign are appropriately ignored.    

Please keep in mind that any parameters specified after the -config flag will overwrite those parameters specified in the -config file. Similarly, any parameters and flags specified before the -config flag will overwrite the config file.

##Utilizing openMP and MPI
Tfit is written using openMP and MPI to perform massive parallelization. If your institution has a large compute cluster, than Tfit will operate well across multiple cores and nodes. To invoke 4 MPI processes run:

$mpirun -np 4 Tfit bidir -config \</path/to/config/file.txt

To invoke Tfit to compute across 3 specific cores or CPUS, invoke.

$Tfit bidir -config \</path/to/config/file.txt -cores 3

Please keep in mind that if you are running on a single node machine. mpirun will utilize CPUs and thus the user should take into account overhead of specifying -cores and -np. The default options for -np and -cores are both 1 respectively. 

If your institution has a compute cluster, please consult your IT staff about the specific job allocation software. If you are using Torque/Maui where computation resources can be specified by PBS directives, then the below job script is sufficient to run Tfit across multiple nodes. 

![Alt text](https://github.com/azofeifa/Tfit/blob/master/images/JobSubmissionExample.png)

Again, this highly system dependent. But openMPI is a well maintained library with lots of resources to help aid in getting Tfit up and running on your compute cluster, please consult https://www.open-mpi.org/faq/ for further reference. 


##Bleeding Edge Parameters 
Tfit is on going effort in development both in terms of software and additions to the mathematics behind the mixture model inference procedure. Currently, I have two additions to the model that have led to better performance. However unlike the methodology we have outlined above, these additions are not published (paper in preparation).  

The first addition is maximum a-posteriori estimation of the parameters lambda, sigma, pi and w. These are achieved in the normal fashion by conjugate priors. These are given below and are only utilized in the "model" module. 

1. -ALPHA_0 prior parameter 1 for sigma (recommended to be set at 1)
2. -BETA_0 prior parameter 2 for sigma (recommended to be set at 1)
3. -ALPHA_1 prior parameter 1 for lambda 
4. -BETA_1 prior parameter 2 for lambda 
5. -ALPHA_2 symmetric prior on mixing weights (the higher this value the more the algorithm will attempt to components of equal mixing weights, recommended 100)
6. -ALPHA_3 symmetric prior on the strand bias (the higher this value the more the algorithm will attempt to find bidirectional events with equal strand bias, very useful! recommended 100)

Please note that the priors for sigma and lambda relate to Gamma distribution parameterized shape and rate parameters where the expected value of the gamma is ALPHA_0 /BETA_0 for sigma and ALPHA_1 / BETA_1 for lambda.

The second addition is estimation for the distance between the two radiating peaks of the bidirectional, informally at this point called the foot print parameter. This is currently integrated into the EM algorithm but will be ignored if -foot_print is not invoked or set to zero in either the call to the bidir or model module. 

##Questions and Comments
Please email me at joseph[dot]azofeifa[at]colorado[dot]edu with questions and comments! Or additionally open an "issue" through the GitHub site. 



##References




