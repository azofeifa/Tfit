#PBS -S /bin/bash

#PBS -N DMSO2_3

#PBS -m ae
#PBS -M joseph.azofeifa@colorado.edu

#PBS -e /Users/azofeifa/qsub_errors/EMG/
#PBS -o /Users/azofeifa/qsub_stdo/EMG/

#PBS -l walltime=05:00:00
#PBS -l nodes=16:ppn=64
#PBS -l mem=180gb
hostlist=$( cat $PBS_NODEFILE | sort | uniq | tr '\n' ',' | sed -e 's/,$//' )
# -- OpenMP environment variables --
OMP_NUM_THREADS=64
export OMP_NUM_THREADS
module load gcc_4.9.2
module load mpich_3.1.4

### ===================
### what machine?
### ===================


src=/Users/azofeifa/Lab/EMG/CPP_src/EMGU
config_file=/Users/azofeifa/Lab/EMG/cpp_config_files/bidir_config.txt
bedgraph_directory=/Users/azofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/
interval_directory=/Users/azofeifa/Lab/gro_seq_files/HCT116/FStitch/
out_directory=/Users/azofeifa/Lab/gro_seq_files/HCT116/EMG_out_files/DMSO2_3_model_fits/
tmp_log_directory=/Users/azofeifa/EMG_log_files/

forward_bedgraph=DMSO2_3.pos.BedGraph
reverse_bedgraph=DMSO2_3.neg.BedGraph
interval_file=DMSO2_3_Nutlin2_3_merged_FStitch.bed
# -- OpenMP environment variables --

# -- program invocation here --
chrom=all
cmd="mpirun -np $PBS_NUM_NODES -hosts ${hostlist}"


$cmd $src $config_file -i ${bedgraph_directory}$forward_bedgraph -j ${bedgraph_directory}$reverse_bedgraph -f ${interval_directory}$interval_file -chr $chrom   -o $out_directory -log_out $tmp_log_directory