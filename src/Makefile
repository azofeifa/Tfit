#======================================================================================/
#**    Author: Joey Azofeifa                                               	       /
#**    Date:   7/8/2015	                                                               /
#**    Description: Makefile for Tfit (Exponentially Modified Gaussian Mixture Model)  /
#======================================================================================/

CXX             = mpic++
CXXFLAGS        = -Wno-unused-variable -Wno-non-virtual-dtor -std=c++11 -fopenmp -Wno-write-strings -Wno-literal-suffix 
EXEC            = ${PWD}/Tfit
ARCH  		= getconf LONG_BIT

CPP_FLAGS_32  	= -D32_BIT
CPP_FLAGS_64  	= -D64_BIT

GCCVERSION 	= $(shell ${CXX} -dumpversion)
NU_FIT: main.o load.o split.o model.o across_segments.o template_matching.o \
	read_in_parameters.o model_selection.o error_stdo_logging.o \
	MPI_comm.o  density_profiler.o bootstrap.o bidir_main.o model_main.o select_main.o FDR.o BIC.o
	@printf "linking           : "
	@${CXX} ${CXXFLAGS}  ${PWD}/main.o ${PWD}/load.o ${PWD}/model_selection.o \
	${PWD}/split.o ${PWD}/model.o ${PWD}/across_segments.o  \
	${PWD}/template_matching.o ${PWD}/read_in_parameters.o  \
	${PWD}/MPI_comm.o   \
	${PWD}/bootstrap.o ${PWD}/density_profiler.o \
	${PWD}/bidir_main.o ${PWD}/model_main.o ${PWD}/BIC.o ${PWD}/FDR.o  \
	${PWD}/select_main.o  ${PWD}/error_stdo_logging.o -o ${EXEC} -lmpi
	@cp ${PWD}/Tfit ${PWD}/EMGU
	@printf "done\n"
	@echo "========================================="
	@printf "successfully compiled \n\n"

main.o: 
	@printf "\n=========================================\n"
	@echo "GCC version: "${GCCVERSION}
	@printf "main              : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/main.cpp 
	@printf "done\n"
error_stdo_logging.o:
	@printf "error_stdo_logging: "
	@${CXX} -c ${CXXFLAGS} ${PWD}/error_stdo_logging.cpp 
	@printf "done\n"

BIC.o:
	@printf "BIC               : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/BIC.cpp
	@printf "done\n"


select_main.o:
	@printf "select_main       : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/select_main.cpp 
	@printf "done\n"

model_main.o:
	@printf "model_main        : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/model_main.cpp 
	@printf "done\n"

bidir_main.o:
	@printf "bidir_main        : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/bidir_main.cpp 
	@printf "done\n"

bootstrap.o:
	@printf "bootstrap         : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/bootstrap.cpp 
	@printf "done\n"

MPI_comm.o: 
	@printf "MPI_comm          : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/MPI_comm.cpp 
	@printf "done\n"

load.o: 
	@printf "load              : "	
	@${CXX} -c ${CXXFLAGS} ${PWD}/load.cpp 
	@printf "done\n"
model_selection.o:
	@printf "model_selection   : "	
	@${CXX} -c ${CXXFLAGS} ${PWD}/model_selection.cpp 
	@printf "done\n"
split.o:
	@printf "split             : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/split.cpp 
	@printf "done\n"

model.o:	
	@printf "model             : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/model.cpp 
	@printf "done\n"
across_segments.o:
	@printf "across_segments   : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/across_segments.cpp 
	@printf "done\n"

template_matching.o:
	@printf "template_matching : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/template_matching.cpp 
	@printf "done\n"
read_in_parameters.o:
	@printf "read_in_parameters: "
	@${CXX} -c ${CXXFLAGS} ${PWD}/read_in_parameters.cpp 
	@printf "done\n"
model_single.o:
	@printf "model_single      : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/model_single.cpp 
	@printf "done\n"
density_profiler.o:
	@printf "density_profiler  : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/density_profiler.cpp 
	@printf "done\n"
FDR.o:
	@printf "FDR               : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/FDR.cpp
	@printf "done\n"

clean:
	@rm -f *.o
	@rm -f EMGU
	@rm -f Tfit	
	@printf "\nsuccessfully removed binaries\n\n"