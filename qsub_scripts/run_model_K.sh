#PBS -S /bin/bash

#PBS -N Allen2014_DMSO2_3_model_K

#PBS -m ae
#PBS -M joseph.azofeifa@colorado.edu

#PBS -e /Users/azofeifa/qsub_errors/EMG/
#PBS -o /Users/azofeifa/qsub_stdo/EMG/

#PBS -l walltime=24:00:00
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
config_file=/Users/azofeifa/Lab/EMG/cpp_config_files/model_config.txt
tmp_log_directory=/Users/azofeifa/EMG_log_files/
hg18=/Users/azofeifa/Lab/genome_files/RefSeqHG18.bed
hg19=/Users/azofeifa/Lab/genome_files/RefSeqHG19.bed
DMSO2_3_Nutlin_2_3_FStitch=/Users/azofeifa/Lab/gro_seq_files/Allen2014/interval_files/DMSO2_3_Nutlin2_3_merged_FS_labeled.bed
#================================================================
#bedgraph directories
Allen2014=/Users/azofeifa/Lab/gro_seq_files/Allen2014/bed_graph_files/
#================================================================
Allen2014_DMSO2_3_forward=DMSO2_3.pos.BedGraph
Allen2014_DMSO2_3_reverse=DMSO2_3.neg.BedGraph
Allen2014_DMSO1027_1212_forward=DMSO1027_1212.pos.BedGraph
Allen2014_DMSO1027_1212_reverse=DMSO1027_1212.pneg.BedGraph
Allen2014_DMSO1hr101911_forward=DMSO1hr101911_lane8.pos.BedGraph
Allen2014_DMSO1hr101911_reverse=DMSO1hr101911_lane8.pneg.BedGraph
Allen2014_DMSOp531027_1212_forward=DMSOp531027_1212.pos.BedGraph
Allen2014_DMSOp531027_1212_reverse=DMSOp531027_1212.pneg.BedGraph
Allen2014_Ma6_NoIndex_forward=Ma6_NoIndex_L008_R1_001.pos.BedGraph
Allen2014_Ma6_NoIndex_reverse=Ma6_NoIndex_L008_R1_001.pneg.BedGraph
Allen2014_Nutlin2_3_forward=Nutlin2_3.sorted.pos.BedGraph
Allen2014_Nutlin2_3_reverse=Nutlin2_3.sorted.neg.BedGraph
Allen2014_Nutlin30_30_forward=Nutlin30_30.pos.BedGraph
Allen2014_Nutlin30_30_reverse=Nutlin30_30.pneg.BedGraph
#================================================================
#out and temp directory
EMG_out=/Users/azofeifa/Lab/gro_seq_files/
Allen2014_out=${EMG_out}Allen2014/EMG_out_files/
#================================================================
#FINAL PARAMETERS

bedgraph_directory=$Allen2014
forward=$Allen2014_DMSO2_3_forward
reverse=$Allen2014_DMSO2_3_reverse
out_directory=$Allen2014_out
intervals=$DMSO2_3_Nutlin_2_3_FStitch

#================================================================
#calling command
cmd="mpirun -np $PBS_NUM_NODES -hosts ${hostlist}"
$cmd $src $config_file -i ${bedgraph_directory}$forward -j ${bedgraph_directory}$reverse -k $intervals -o $out_directory -log_out $tmp_log_directory -N $PBS_JOBNAME 
#================================================================j