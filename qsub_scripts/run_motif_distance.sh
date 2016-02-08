#PBS -S /bin/bash

#PBS -N motif_search
#PBS -m ae
#PBS -M joseph.azofeifa@colorado.edu

#PBS -e /Users/azofeifa/qsub_errors/EMG/
#PBS -o /Users/azofeifa/qsub_stdo/EMG/

#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
module load scipy_0.14.0
module load matplotlib_1.3.1
module load networkx_1.9.1
module load numpy_1.9.2
root=/Users/azofeifa/FIMO_OUT/
query=/Users/azofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/Allen2014_DMSO2_3-1_bidirectional_hits_intervals.bed
out=/Users/azofeifa/TF_predictions/Allen2014_DMSO2_3_motif_2000
pad=2000

#root=/Users/azofeifa/FIMO_OUT/
#query=/Users/azofeifa/temp_ENCODE_files/wgEncodeUwDnaseHct116PkRep1.narrowPeak
#out=/Users/azofeifa/DNAse_motif_overlap.bed

#python /Users/azofeifa/Lab/EMG/TF_predictions/get_DNAse_overlaps.py $root $query $out
python /Users/azofeifa/Lab/EMG/TF_predictions/assign_TF.py $root $query $out $pad