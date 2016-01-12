#PBS -S /bin/bash

#PBS -N Pol_II_DMSO
#PBS -m ae
#PBS -M joseph.azofeifa@colorado.edu

#PBS -e /Users/azofeifa/qsub_errors/EMG/
#PBS -o /Users/azofeifa/qsub_stdo/EMG/

#PBS -l walltime=24:00:00
#PBS -l nodes=32:ppn=64
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
config_file=/Users/azofeifa/Lab/EMG/cpp_config_files/single_config.txt
tmp_log_directory=/Users/azofeifa/EMG_log_files/

#================================================================
#INPUT files
bedgraph_dir=/Users/azofeifa/Lab/ChIP/HCT116/Pol_II/
interval_dir=/Users/azofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/

bedgraph_file=Pol2_DMSO_HG19_150bp_genomeCovPDMMR.BedGraph
interval_file=DMSO2_3-1_bidirectional_hits_intervals.bed

#================================================================
#OUTPUT diretories
out_directory=/Users/azofeifa/Lab/ChIP/HCT116/EMG_out_files/
tmp_log_directory=/Users/azofeifa/EMG_log_files/


#================================================================
#calling program
cmd="mpirun -np $PBS_NUM_NODES -hosts ${hostlist}"
$cmd $src $config_file -i ${bedgraph_dir}$bedgraph_file -j ${interval_dir}$interval_file   -o $out_directory -log_out $tmp_log_directory -N $PBS_JOBNAME 