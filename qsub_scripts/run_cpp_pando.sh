#Name the job
#PBS -N EMG_Model_Fit

### Specify the number of nodes/cores
#PBS -l nodes=1:ppn=64

### Allocate the amount of memory needed
#PBS -l mem=10gb

### Set your expected walltime
#PBS -l walltime=12:00:00

### Setting to mail when the job is complete
#PBS -e /Users/azofeifa/qsub_errors/EMG/                                                                                              
#PBS -o /Users/azofeifa/qsub_stdo/EMG/  

### Set your email address
#PBS -M jgazofeifa@gmail.com

###
#PBS -t 1-64


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
    module load gcc_4.9.2 ###load most recent gcc compiler
fi


src=/Users/azofeifa/Lab/EMG/CPP_src/EMGU
config_file=/Users/azofeifa/Lab/EMG/cpp_config_files/model_config.txt
formatted_in_file=/Users/azofeifa/Lab/gro_seq_files/HCT116/EMG_out_files/fs_only_merged_32.tsv
out_directory=/Users/azofeifa/Lab/gro_seq_files/HCT116/EMG_out_files/cpp_model_fits_03/
maxK=5
rounds=64
r_mu=3
mi=1000
$src $config_file  -i $formatted_in_file -o $out_directory -np $PBS_ARRAYID -rounds $PBS_ARRAYID -r_mu $r_mu -mi $mi


