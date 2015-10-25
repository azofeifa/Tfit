#PBS -S /bin/bash

#PBS -N Allen2014_boot_DMSO2_3

#PBS -m ae
#PBS -M joseph.azofeifa@colorado.edu

#PBS -e /Users/azofeifa/qsub_errors/EMG/
#PBS -o /Users/azofeifa/qsub_stdo/EMG/

#PBS -l walltime=12:00:00
#PBS -l nodes=40:ppn=64
#PBS -l mem=800gb

hostlist=$( cat $PBS_NODEFILE | sort | uniq | tr '\n' ',' | sed -e 's/,$//' )
# -- OpenMP environment variables --
OMP_NUM_THREADS=64
export OMP_NUM_THREADS
module load gcc_4.9.2
module load mpich_3.1.4

#================================================================
#paths to config and src
src=/Users/azofeifa/Lab/EMG/CPP_src/EMGU
config_file=/Users/azofeifa/Lab/EMG/cpp_config_files/bootstrap_config.txt
tmp_log_directory=/Users/azofeifa/EMG_log_files/
#================================================================
#FINAL PARAMETERS

bedgraph_directory=/Users/azofeifa/Lab/gro_seq_files/Allen2014/bed_graph_files/
MLE=/Users/azofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/DMSO2_3-1_bidirectional_hits_intervals.bed
forward=DMSO2_3.pos.BedGraph
reverse=DMSO2_3.neg.BedGraph
out_directory=/Users/azofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/




#================================================================
#calling command
cmd="mpirun -np $PBS_NUM_NODES -hosts ${hostlist}"
$cmd $src $config_file -i $MLE  -j ${bedgraph_directory}$forward -k ${bedgraph_directory}$reverse  -o $out_directory -log_out $tmp_log_directory -N $PBS_JOBNAME 
#================================================================