#Name the job
#PBS -N EMG_formatting

#PBS -l nodes=1:ppn=1

### Allocate the amount of memory needed
#PBS -l pmem=1gb

### Set your expected walltime
#PBS -l walltime=12:00:00

### Setting to mail when the job is complete
#PBS -e /Users/azofeifa/qsub_errors/EMG/                                                                                              
#PBS -o /Users/azofeifa/qsub_stdo/EMG/  

### Set your email address
#PBS -m a
#PBS -M jgazofeifa@gmail.com


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
config_file=/Users/azofeifa/Lab/EMG/cpp_config_files/format_config.txt
EMG_out_directory=/Users/azofeifa/Lab/gro_seq_files/HCT116/EMG_out_files/EMG_formmated_files/
#interval_directory=/Users/azofeifa/Lab/gro_seq_files/HCT116/interval_files/
interval_directory=/Users/azofeifa/ENCODE/HCT116/Sp1/peak_files/
bedgraph_directory=/Users/azofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/
out_directory=/Users/azofeifa/Lab/gro_seq_files/HCT116/EMG_out_files/EMG_formmated_files/

interval_file=SL12239_Peaks.bed.broadPeak
forward_bedgraph=DMSO2_3.pos.BedGraph
reverse_bedgraph=DMSO2_3.neg.BedGraph
out_file=DMSO2_3_Sp1_rep1.tsv

$src $config_file -i ${interval_directory}$interval_file -j ${bedgraph_directory}$forward_bedgraph -k ${bedgraph_directory}$reverse_bedgraph -o ${out_directory}$out_file