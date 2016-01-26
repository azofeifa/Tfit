#PBS -S /bin/bash


#PBS -m ae
#PBS -M joseph.azofeifa@colorado.edu

#PBS -e /Users/azofeifa/qsub_errors/EMG/
#PBS -o /Users/azofeifa/qsub_stdo/EMG/

#PBS -l walltime=48:00:00
#PBS -l nodes=10:ppn=64
#PBS -l mem=200gb

hostlist=$( cat $PBS_NODEFILE | sort | uniq | tr '\n' ',' | sed -e 's/,$//' )
# -- OpenMP environment variables --
OMP_NUM_THREADS=64
export OMP_NUM_THREADS
module load gcc_4.9.2
module load mpich_3.1.4

#================================================================
#paths to config and src
src=/Users/azofeifa/Lab/EMG/CPP_src/EMGU
config_file=/Users/azofeifa/Lab/EMG/cpp_config_files/config_file.txt
tmp_log_directory=/Users/azofeifa/EMG_log_files/

#================================================================
#calling command
cmd="mpirun -np $PBS_NUM_NODES -hosts ${hostlist}"
$cmd $src bidir -config $config_file -MLE 1   -i $forward -j $reverse  -o $out_directory -log_out $tmp_log_directory -N $name -tss $genome
#================================================================