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

python /Users/azofeifa/Lab/EMG/analysis/assign_TF.py