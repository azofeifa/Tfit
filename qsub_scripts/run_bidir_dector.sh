#PBS -S /bin/bash

#PBS -N Allen2014_Nutlin30_30

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
config_file=/Users/azofeifa/Lab/EMG/cpp_config_files/bidir_config.txt
tmp_log_directory=/Users/azofeifa/EMG_log_files/
hg18=/Users/azofeifa/Lab/genome_files/hg18_no_anntation.bed
hg19=/Users/azofeifa/Lab/genome_files/hg19_no_anntation.bed
#================================================================
#bedgraph directories
pubgro=/projects/dowellLab/groseq/pubgro/
bedgraph_directory=/Users/azofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/ #allen,HCT116
Core2014=${pubgro}Core2014/bedGraph_from_bw/
Danko2014=${pubgro}Danko2014/bedGraph_from_bw/
Li2013=${pubgro}Li2013/forFStitch/
Liu2013=${pubgro}Liu2013/bedGraph_from_bw/
Puc2015=${pubgro}Puc2015/bedGraphs_from_bedGraph/
Meng2014b=${pubgro}Meng2014b/
Zhang2015=${pubgro}Zhang2015/
Allen2014=/Users/azofeifa/Lab/gro_seq_files/Allen2014/bed_graph_files/
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

Liu2013_siBRD4_forward=siBRD4-1-re1_2-HEK293T-forward.bedGraph
Liu2013_siBRD4_reverse=siBRD4-1-re1_2-HEK293T-reverse.bedGraph
Liu2013_siCTL_forward=siCTL-re1_2-HEK293T-forward.bedGraph
Liu2013_siCTL_reverse=siCTL-re1_2-HEK293T-reverse.bedGraph
Liu2013_siJMJD6_forward=siJMJD6-1-re1_2-HEK293-forward.bedGraph
Liu2013_siJMJD6_reverse=siJMJD6-1-re1_2-HEK293-reverse.bedGraph

Puc2015_si_control_vehicle_forward=siControl_1h_vehicle_hg18_forward.bedGraph
Puc2015_si_control_vehicle_reverse=siControl_1h_vehicle_hg18_reverse.bedGraph
Puc2015_siMRE11_dht_forward=siMRE11_1h_dht_hg18_forward.bedGraph
Puc2015_siMRE11_dht_reverse=siMRE11_1h_dht_hg18_reverse.bedGraph
Puc2015_siMRE11_vehicle_forward=siMRE11_1h_vehicle_hg18_forward.bedGraph
Puc2015_siMRE11_vehicle_reverse=siMRE11_1h_vehicle_hg18_reverse.bedGraph
Puc2015_siTOP1_vehicle_forward=siTOP1_1h_vehicle_hg18_forward.bedGraph
Puc2015_siTOP1_vehicle_reverse=siTOP1_1h_vehicle_hg18_reverse.bedGraph
Puc2015_siTOP1_dht_forward=siTOP1_1h_dht_hg18_forward.bedGraph
Puc2015_siTOP1_dht_reverse=siTOP1_1h_dht_hg18_reverse.bedGraph

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
Core2014_out=${EMG_out}Core_2014/EMG_out_files/
Danko_2014_out=${EMG_out}Danko_2014/EMG_out_files/
Allen2014_out=${EMG_out}Allen2014/EMG_out_files/
Li2013_out=${EMG_out}Li2013/EMG_out_files/
Liu2013_out=${EMG_out}Liu2013/EMG_out_files/
Puc2015_out=${EMG_out}Puc2015/EMG_out_files/
Meng2014b_out=${EMG_out}Meng2014b/EMG_out_files/
Zhang2015_out=${EMG_out}Zhang2015/EMG_out_files/
#================================================================
#FINAL PARAMETERS

bedgraph_directory=$Allen2014
forward=$Allen2014_Nutlin30_30_forward
reverse=$Allen2014_Nutlin30_30_reverse
out_directory=$Allen2014_out
genome=$hg19



#================================================================
#calling command
cmd="mpirun -np $PBS_NUM_NODES -hosts ${hostlist}"
$cmd $src $config_file -i ${bedgraph_directory}$forward -j ${bedgraph_directory}$reverse -f ${interval_directory}$interval_file   -o $out_directory -log_out $tmp_log_directory -N $PBS_JOBNAME -nf $genome
#================================================================