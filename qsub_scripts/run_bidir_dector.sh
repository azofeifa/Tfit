#PBS -S /bin/bash

#PBS -N Li2013_MCF7_siSMC3_EtoH

#PBS -m ae
#PBS -M joseph.azofeifa@colorado.edu

#PBS -e /Users/azofeifa/qsub_errors/EMG/
#PBS -o /Users/azofeifa/qsub_stdo/EMG/

#PBS -l walltime=12:00:00
#PBS -l nodes=10:ppn=64
#PBS -l mem=700gb

hostlist=$( cat $PBS_NODEFILE | sort | uniq | tr '\n' ',' | sed -e 's/,$//' )
# -- OpenMP environment variables --
OMP_NUM_THREADS=64
export OMP_NUM_THREADS
module load gcc_4.9.2
module load mpich_3.1.4

#================================================================
#paths to config and src
src=/Users/azofeifa/Lab/EMG/CPP_src/EMGU
config_file=/Users/azofeifa/Lab/EMG/cpp_config_files/bidir_config.txt
#================================================================
#Interval Files not needed yet
interval_directory=/Users/azofeifa/Lab/gro_seq_files/HCT116/FStitch/
tmp_log_directory=/Users/azofeifa/EMG_log_files/
interval_file=DMSO2_3_Nutlin2_3_merged_FStitch.bed
#================================================================
#bedgraph directories
pubgro=/projects/dowellLab/groseq/pubgro/
bedgraph_directory=/Users/azofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/ #allen,HCT116
Core2014=${pubgro}Core2014/bedGraph_from_bw/
Danko2014=${pubgro}Danko2014/bedGraph_from_bw/
Li2013=${pubgro}Li2013/forFStitch/
Liu2013=${pubgro}Liu2013/
Meng2014b=${pubgro}Meng2014b/
Zhang2015=${pubgro}Zhang2015/
#================================================================
#bedgraph files

core_2014_K562_forward=GSM1480325_K562_GROseq_plus.bigWig.bedGraph.mod.bedGraph
core_2014_K562_reverse=GSM1480325_K562_GROseq_minus.bigWig.bedGraph.mod.bedGraph
core_2014_GM12878_forward=GSM1480326_GM12878_GROseq_plus.bigWig.bedGraph.mod.bedGraph
core_2014_GM12878_reverse=GSM1480326_GM12878_GROseq_minus.bigWig.bedGraph.mod.bedGraph

danko_2014_ac16_unt_forward=GSE66031_ac16.unt.all_plus.bw.bedGraph.mod.bedGraph
danko_2014_ac16_unt_reverse=GSE66031_ac16.unt.all_minus.bw.bedGraph.mod.bedGraph
danko_2014_MCF7_40m_forward=GSE66031_MCF7.40m.all_plus.bw.bedGraph.mod.bedGraph
danko_2014_MCF7_40m_reverse=GSE66031_MCF7.40m.all_minus.bw.bedGraph.mod.bedGraph
danko_2014_MCF7_unt_forward=GSE66031_MCF7.unt.all_plus.bw.bedGraph.mod.bedGraph
danko_2014_MCF7_unt_reverse=GSE66031_MCF7.unt.all_minus.bw.bedGraph.mod.bedGraph
danko_2014_cd4_forward=GSM1613181_cd4_plus.bw.bedGraph.mod.bedGraph
danko_2014_cd4_reverse=GSM1613181_cd4_minus.bw.bedGraph.mod.bedGraph
danko_2014_Jurkat_forward=GSM1613182_Jurkat_plus.bw.bedGraph.mod.bedGraph
danko_2014_Jurkat_reverse=GSM1613182_Jurkat_minus.bw.bedGraph.mod.bedGraph

Li2013_MCF7_E2_1_2_forward=GSM1115996_Groseq-MCF7-E2-rep1_2.forward.bedGraph
Li2013_MCF7_E2_1_2_reverse=GSM1115996_Groseq-MCF7-E2-rep1_2.reverse.bedGraph
Li2013_MCF7_EtoH_1_2_forward=GSM1115997_Groseq-MCF7-EtoH-rep1_2.forward.bedGraph
Li2013_MCF7_EtoH_1_2_reverse=GSM1115997_Groseq-MCF7-EtoH-rep1_2.reverse.bedGraph
Li2013_MCF7_siSMC3_E2_forward=GSM1115999_Groseq-MCF7-siSMC3-E2.plusstrand.bigWig.bedGraph
Li2013_MCF7_siSMC3_E2_reverse=GSM1115999_Groseq-MCF7-siSMC3-E2.minusstrand.bigWig.bedGraph.plus.bedGraph
Li2013_MCF7_siSMC3_EtoH_forward=GSM1116000_Groseq-MCF7-siSMC3-EtoH.plusstrand.bigWig.bedGraph
Li2013_MCF7_siSMC3_EtoH_reverse=GSM1116000_Groseq-MCF7-siSMC3-EtoH.minusstrand.bigWig.bedGraph.plus.bedGraph



#================================================================
#Final Directory and bedgraph directories
bedgraph_directory=$Li2013
forward=$Li2013_MCF7_siSMC3_EtoH_forward
reverse=$Li2013_MCF7_siSMC3_EtoH_reverse

#================================================================
#out and temp directory
EMG_out=/Users/azofeifa/Lab/gro_seq_files/
Core2014_out=${EMG_out}Core_2014/
Danko_2014_out=${EMG_out}Danko_2014/
HCT116_out=${EMG_out}HCT116/
Li2013_out=${EMG_out}Li2013/
Liu2013_out=${EMG_out}Liu2013/
Meng2014b_out=${EMG_out}Meng2014b/
Zhang2015_out=${EMG_out}Zhang2015/
#================================================================
#FINAL PARAMETERS
out_directory=$Li2013_out





#================================================================
#calling command
cmd="mpirun -np $PBS_NUM_NODES -hosts ${hostlist}"
$cmd $src $config_file -i ${bedgraph_directory}$forward -j ${bedgraph_directory}$reverse -f ${interval_directory}$interval_file   -o $out_directory -log_out $tmp_log_directory 
#================================================================