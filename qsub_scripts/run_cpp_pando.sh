#Name the job
#PBS -N EMG_Model_Fit

### Specify the number of nodes/cores
#PBS -l nodes=1:ppn=64

### Allocate the amount of memory needed
#PBS -l pmem=1gb

### Set your expected walltime
#PBS -l walltime=48:00:00

### Setting to mail when the job is complete
#PBS -e /Users/azofeifa/qsub_errors/EMG/                                                                                              
#PBS -o /Users/azofeifa/qsub_stdo/EMG/  

### Set your email address
#PBS -m a
#PBS -M jgazofeifa@gmail.com

#PBS -t 1-22


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
formatted_in_file=/Users/azofeifa/Lab/gro_seq_files/HCT116/EMG_out_files/EMG_formmated_files/DMSO_ND_intervals_formatted.tsv
out_directory=/Users/azofeifa/Lab/gro_seq_files/HCT116/EMG_out_files/DMSO_ND_intervals_model_fits_2/
chr=chr$PBS_ARRAYID
$src $config_file  -i $formatted_in_file -o $out_directory -chr $chr


