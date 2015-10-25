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

root=/Users/azofeifa/FIMO_OUT/
query=/Users/azofeifa/Lab/gro_seq_files/Yang2013/EMG_out_files/SRR892016_Yang2013-2_bidirectional_hits_intervals.bed
out=/Users/azofeifa/TF_predictions/SRR892016_Yang2013-2

python /Users/azofeifa/Lab/EMG/analysis/assign_TF.py $root $query $out