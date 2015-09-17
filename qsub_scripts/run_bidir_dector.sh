#Name the job
#PBS -N DMSO2_3

#PBS -l nodes=7:ppn=32,walltime=24:00:00,mem=100gb
### Allocate the amount of memory needed

### Set your expected walltime
###PBS -q shortmem
### Setting to mail when the job is complete
#PBS -e /Users/azofeifa/qsub_errors/EMG/                                        
#PBS -o /Users/azofeifa/qsub_stdo/EMG/  

### Set your email address
#PBS -m ae
#PBS -M jgazofeifa@gmail.com

#PBS -W x=nmatchpolicy=exactnode

### Choose your shell 
#PBS -S /bin/sh
### Pass enviroment variables to the job
#PBS -V



### ===================
### what machine?
### ===================

module load gcc_4.9.2 ###load most recent gcc compiler
module load openmpi_1.6.4
cd $PBS_O_WORKDIR

src=/Users/azofeifa/Lab/EMG/CPP_src/EMGU
config_file=/Users/azofeifa/Lab/EMG/cpp_config_files/bidir_config.txt
bedgraph_directory=/Users/azofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/
interval_directory=/Users/azofeifa/Lab/gro_seq_files/HCT116/FStitch/
out_directory=/Users/azofeifa/Lab/gro_seq_files/HCT116/EMG_out_files/DMSO2_3_model_fits/

forward_bedgraph=DMSO2_3.pos.BedGraph
reverse_bedgraph=DMSO2_3.neg.BedGraph
interval_file=DMSO2_3_Nutlin2_3_merged_FStitch.bed
# -- OpenMP environment variables --
OMP_NUM_THREADS=$PBS_NUM_PPN
export OMP_NUM_THREADS
# -- program invocation here --
chrom=chr1

mpirun --npernode 1 $src $config_file -i ${bedgraph_directory}$forward_bedgraph -j ${bedgraph_directory}$reverse_bedgraph -f ${interval_directory}$interval_file -chr $chrom   -o $out_directory 