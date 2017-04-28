#! /bin/bash
#===========================================================
install_location=/Users/joazofeifa/Lab/ #whereever you installed MDS
#change this accordingly
#===========================================================
#input parameters and necessary paths

src=${install_location}/Tfit/src/Tfit #path to src

input_bedgraph=${install_location}/Tfit/examples/test_ij_joint.BedGraph
config_file=${install_location}/Tfit/examples/config_file.txt
scores=${install_location}/Tfit/examples/scores.bedgraph
echo '-------------------------Running Unit Tests-------------------------'

echo 'Installation Location: ' $install_location ' is this correct?' 


mpirun -np 1 $src bidir -config $config_file -ij $input_bedgraph -o ${install_location}/Tfit/examples/ -log_out ${install_location}/Tfit/examples/  -pad 2000 -bct 0.1







