#PBS -S /bin/bash

#PBS -N Danko2013_merged
#PBS -m ae
#PBS -M joseph.azofeifa@colorado.edu

#PBS -e /Users/azofeifa/qsub_errors/EMG/
#PBS -o /Users/azofeifa/qsub_stdo/EMG/

#PBS -l walltime=12:00:00
#PBS -l nodes=20:ppn=64
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
#FINAL PARAMETERS

bedgraph_directory=/projects/dowellLab/groseq/pubgro/Danko2013/mapped/bowtie2/sortedbam/forFstitch/
forward=mergeGRODanko2013.fiveprime.pos.BedGraph
reverse=mergeGRODanko2013.fiveprime.neg.BedGraph
out_directory=/projects/dowellLab/TFIT/Danko2013/EMG_out_files/
genome=$hg19



#================================================================
#calling command
cmd="mpirun -np $PBS_NUM_NODES -hosts ${hostlist}"
$cmd $src $config_file -i ${bedgraph_directory}$forward -j ${bedgraph_directory}$reverse -f ${interval_directory}$interval_file   -o $out_directory -log_out $tmp_log_directory -N $PBS_JOBNAME -nf $genome -rounds 65
#================================================================