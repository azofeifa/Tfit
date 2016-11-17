#PBS -S /bin/bash

#PBS -N Tfit
#PBS -m ae


#PBS -l walltime=36:00:00
#PBS -l nodes=1:ppn=64
#PBS -l mem=100gb

hostlist=$( cat $PBS_NODEFILE | sort | uniq | tr '\n' ',' | sed -e 's/,$//' )
# -- OpenMP environment variables --
OMP_NUM_THREADS=64
export OMP_NUM_THREADS
module load gcc_4.9.2
module load mpich_3.1.4



#================================================================
#paths to config and src
src=/path/to/src/Tfit
config_file=/path/to/config_file.txt


#================================================================
#calling command
cmd="mpirun -np $PBS_NUM_NODES -hosts ${hostlist}"
$cmd $src bidir -config $config_file -ij $bg_file  -tss $TSS -bct $bct  -o $out_directory -log_out $tmp_log_directory -N $name
#================================================================