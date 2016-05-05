src=/Users/joazofeifa/Lab/EMG/CPP_src/EMGU
config_file=/Users/joazofeifa/Lab/EMG/CPP_src/config_file.txt
ij=/Users/joazofeifa/Lab/gro_seq_files/subsampled/0.462_rep_0_Allen2014_DMSO.bedgraph
NP=3


mpirun -np $NP $src bidir -config $config_file -ij $ij