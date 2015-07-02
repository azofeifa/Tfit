#Name the job
#PBS -N EMG_Formatter
### Choose the queue ('short' = <24hrs, 'long' = >24hrs)
#PBS -q short

### Specify the number of nodes/cores
#PBS -l nodes=1:ppn=1

### Allocate the amount of memory needed
#PBS -l mem=64mb

### Set your expected walltime
#PBS -l walltime=05:00:00

### Setting to mail when the job is complete
#PBS -e /Users/azofeifa/qsub_errors/EMG/                                                                                              
#PBS -o /Users/azofeifa/qsub_stdo/EMG/  

### Set your email address
#PBS -M jgazofeifa@gmail.com

### Choose your shell 
#PBS -S /bin/sh
### Pass enviroment variables to the job
#PBS -V
###Load Modules
module load matplotlib_1.3.1
module load numpy_1.9.2
module load scipy_0.12.0

root=/Users/azofeifa/Lab/
src=/Users/azofeifa/Lab/EMG/
ref=${root}genome_files/RefSeqHG19.txt
ffs=${root}gro_seq_files/HCT116/FStitch/DMSO2_3.sorted.fiveprime.pos_segs_IGV.bed
rfs=${root}gro_seq_files/HCT116/FStitch/DMSO2_3.sorted.fiveprime.neg_segs_IGV.bed
fbg=${root}gro_seq_files/HCT116/bed_graph_files/DMSO2_3.pos.BedGraph
rbg=${root}gro_seq_files/HCT116/bed_graph_files/DMSO2_3.neg.BedGraph
wo=${root}gro_seq_files/HCT116/EMG_out_files
pad=100



python ${src}python_src/ formatData FStitchSingleIsoform  -ref $ref -ffs $ffs -rfs $rfs -fbg $fbg -rbg $rbg -pad $pad -wo $wo
