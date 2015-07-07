#Name the job
#PBS -N EMG_Formatter
### Choose the queue ('short' = <24hrs, 'long' = >24hrs)
#PBS -q long

### Specify the number of nodes/cores
#PBS -l nodes=1:ppn=16

### Allocate the amount of memory needed
#PBS -l mem=10gb

### Set your expected walltime
#PBS -l walltime=24:00:00

### Setting to mail when the job is complete
#PBS -e /Users/azofeifa/qsub_errors/EMG/                                                                                              
#PBS -o /Users/azofeifa/qsub_stdo/EMG/  

### Set your email address
#PBS -M jgazofeifa@gmail.com

### Choose your shell 
#PBS -S /bin/sh
### Pass enviroment variables to the job
#PBS -V

###Job Array for running across different chromosomes
#PBS -t 1-23

### ===================
### what machine?
### ===================
vieques_pando=true ###unix compute clusters
mac=false ###macOS
if [ "$vieques_pando" = true ] ; then ###load modules 
    module load matplotlib_1.3.1
    module load numpy_1.9.2
    module load scipy_0.12.0
fi

if [ "$vieques_pando" = true ] ; then
    root=/Users/azofeifa/Lab/
    src=/Users/azofeifa/Lab/EMG/
else
    src=/Users/joeyazo/Desktop/Lab/EMG/
    root=/Users/joeyazo/Desktop/Lab/
fi
### ====================
### EMG MODULE TYPE
### ====================
format=false
runModel=true

if [ "$format" = true ] ; then

    format_option=FStitchMerged

    echo "EMG: formatting option"
    ref=${root}genome_files/RefSeqHG19.txt
    ffs=${root}gro_seq_files/HCT116/FStitch/DMSO2_3.sorted.fiveprime.pos_segs_IGV.bed
    rfs=${root}gro_seq_files/HCT116/FStitch/DMSO2_3.sorted.fiveprime.neg_segs_IGV.bed
    fbg=${root}gro_seq_files/HCT116/bed_graph_files/DMSO2_3.pos.BedGraph
    rbg=${root}gro_seq_files/HCT116/bed_graph_files/DMSO2_3.neg.BedGraph
    wo=${root}gro_seq_files/HCT116/EMG_out_files
    pad=100



    python ${src}python_src/ formatData $format_option -ref $ref -ffs $ffs -rfs $rfs -fbg $fbg -rbg $rbg -pad $pad -wo $wo
fi

if [ "$runModel" = true ] ; then 
    echo "EMG: model option"
    formatted_file=${root}gro_seq_files/HCT116/EMG_out_files/out_format_file.tsv
    
    wo=${root}gro_seq_files/HCT116/EMG_out_files/
    k=3
    it=16
    bins=300
    sc=chr$PBS_ARRAYID ###specific chromosome
    bic=0 ###perform model selection?
    st=100 ###standardize, numerical stability
    mc=0.0001 ###EM convergence threshold
    mt=300 ###EM max number of iterations, before we give up
    mu=0
    python ${src}python_src/ runModel -i $formatted_file -wo $wo -k $k -it $it -b $bins -sc $sc -bic $bic -st $st -mc $mc -mt $mt -mu $mu
fi