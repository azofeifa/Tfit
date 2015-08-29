#Name the job
#PBS -N bi_detect

#PBS -l nodes=20:ppn=32

### Allocate the amount of memory needed
#PBS -l pmem=1gb

### Set your expected walltime
#PBS -l walltime=12:00:00

### Setting to mail when the job is complete
#PBS -e /Users/azofeifa/qsub_errors/EMG/                                                                                              
#PBS -o /Users/azofeifa/qsub_stdo/EMG/  

### Set your email address
#PBS -m ae
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
config_file=/Users/azofeifa/Lab/EMG/cpp_config_files/bidir_config.txt
bedgraph_directory=/Users/azofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/
out_directory=/Users/azofeifa/Lab/gro_seq_files/HCT116/EMG_out_files/

forward_bedgraph=DMSO2_3.pos.BedGraph
reverse_bedgraph=DMSO2_3.neg.BedGraph
mpi=/opt/openmpi/1.6.4/bin/mpirun
nodes=20
$mpi -np $nodes $src $config_file -i ${bedgraph_directory}$forward_bedgraph -j ${bedgraph_directory}$reverse_bedgraph -o $out_directory