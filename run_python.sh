
root=/Users/joeyazo/Desktop/Lab/gro_seq_files/
ref=${root}RefSeqHG19.txt
ffs=${root}HCT116/FStitch/DMSO2_3.sorted.fiveprime.pos_segs_IGV.bed
rfs=${root}HCT116/FStitch/DMSO2_3.sorted.fiveprime.neg_segs_IGV.bed
fbg=${root}HCT116/bed_graph_files/DMSO2_3.pos.BedGraph
rbg=${root}HCT116/bed_graph_files/DMSO2_3.neg.BedGraph
wo=${root}
pad=100


python python_src/ formatData FStitchSingleIsoform  -ref $ref -ffs $ffs -rfs $rfs -fbg $fbg -rbg $rbg -pad $pad -wo $wo