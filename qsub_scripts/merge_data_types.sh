#Name the job
#PBS -N merging_data_types
### Specify the number of nodes/cores
#PBS -l nodes=1:ppn=1

### Allocate the amount of memory needed
#PBS -l pmem=1gb

### Set your expected walltime
#PBS -l walltime=12:00:00

### Setting to mail when the job is complete
#PBS -e /Users/azofeifa/qsub_errors/                                                                                              
#PBS -o /Users/azofeifa/qsub_stdo/  

### Set your email address
#PBS -m ae
#PBS -M jgazofeifa@gmail.com



### Choose your shell 
#PBS -S /bin/sh
### Pass enviroment variables to the job
#PBS -V



### ===================
### what machine?
### ===================
vieques_pando=true ###unix compute clusters
mac=false ###macOS
if [ "$vieques_pando" = true ] ; then ###load modules 
    module load numpy_1.9.2
    module load scipy_0.14.0
    module load matplotlib_1.3.1
fi


src=/Users/azofeifa/Lab/EMG/analysis/load.py
root=/Users/azofeifa/
model_first_directory=${root}Lab/gro_seq_files/HCT116/EMG_out_files/DMSO_ND_intervals_model_fits_2/
refseq_file=${root}genome_files/RefSeqHG19.txt
dbSNP_directory=${root}dbSNP/
ENCODE_directory=${root}ENCODE/HCT116/
bedgraph_gro_directory=${root}Lab/gro_seq_files/HCT116/bed_graph_files/
out_directory=Lab/gro_seq_files/HCT116/
python $src $model_first_directory $refseq_file $dbSNP_directory $ENCODE_directory $bedgraph_gro_directory $out_directory


