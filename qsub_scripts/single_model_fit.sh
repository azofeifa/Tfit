#Name the job
#PBS -N NU_FIT

#PBS -l nodes=1:ppn=32

### Allocate the amount of memory needed
#PBS -l pmem=9gb

### Set your expected walltime
#PBS -l walltime=12:00:00
###PBS -q shortmem
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
    module load openmpi_1.6.4
fi

src=/Users/azofeifa/Lab/EMG/CPP_src/EMGU
config_file=/Users/azofeifa/Lab/EMG/cpp_config_files/single_config.txt
bedgraph_directory=/Users/azofeifa/Lab/ChIP/HCT116/bedgraph_files/
genome_directory=/Users/azofeifa/Lab/ChIP/HCT116/genome_files/
out_directory=/Users/azofeifa/Lab/ChIP/HCT116/EMG_out_files/

single_bedgraph=tPol_II_DMSO_150bp_genomeCovGraphPDMMR.BedGraph
hg18=hg18_refseq.bed
mpi=/opt/openmpi/1.6.4/bin/mpirun
nodes=5
$mpi $src $config_file -i ${bedgraph_directory}$single_bedgraph -j  ${genome_directory}$hg18  -o $out_directory