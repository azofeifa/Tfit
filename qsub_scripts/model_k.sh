#PBS -S /bin/bash

#PBS -N fp_si
#PBS -m ae
#PBS -M joseph.azofeifa@colorado.edu

#PBS -e /Users/azofeifa/qsub_errors/EMG/
#PBS -o /Users/azofeifa/qsub_stdo/EMG/

#PBS -l walltime=48:00:00
#PBS -l nodes=18:ppn=64
#PBS -l mem=2000gb

hostlist=$( cat $PBS_NODEFILE | sort | uniq | tr '\n' ',' | sed -e 's/,$//' )
# -- OpenMP environment variables --
OMP_NUM_THREADS=64
export OMP_NUM_THREADS
module load gcc_4.9.2
module load mpich_3.1.4
#MAX_Peaks.bed.broadPeak
#SP1_Peaks.bed.broadPeak
#TEAD4_Peaks.bed.broadPeak
#ZBTB33_Peaks.bed.broadPeak
#================================================================
#paths to config and src
src=/Users/azofeifa/Lab/EMG/CPP_src/EMGU
config_file=/Users/azofeifa/Lab/EMG/cpp_config_files/model_config.txt
tmp_log_directory=/Users/azofeifa/EMG_log_files/

BG_DIR=/Users/azofeifa/Lab/gro_seq_files/Allen2014/bed_graph_files/
SIM_BG_DIR=/Users/azofeifa/Lab/simulation_analysis/bedgraph_files/
TF_INT_DIR=/Users/azofeifa/Lab/ChIP/HCT116/TFS/peak_files/



GENE_INT_DIR=/Users/azofeifa/Lab/gro_seq_files/Allen2014/interval_files/
SIM_INT_DIR=/Users/azofeifa/Lab/simulation_analysis/interval_files/

forward=${SIM_BG_DIR}fp_by_si.pos.bedgraph
reverse=${SIM_BG_DIR}fp_by_si.neg.bedgraph

#intervals=${TF_DIR}ZBTB33_Peaks.bed.broadPeak
intervals=${SIM_INT_DIR}fp_by_si_intervals.bed
out=/Users/azofeifa/Lab/ChIP/HCT116/EMG_out_files/
out=/Users/azofeifa/Lab/simulation_analysis/EMG_out_files/
log_out=/Users/azofeifa/EMG_log_files/

#HYPER parameters


#================================================================
#calling command
cmd="mpirun -np $PBS_NUM_NODES -hosts ${hostlist}"
$cmd $src $config_file -i $forward -j $reverse -k $intervals -o $out -log_out $log_out -N $PBS_JOBNAME -pad 0 
#================================================================