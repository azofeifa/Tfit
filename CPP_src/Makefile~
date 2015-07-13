#***************************************************************************/
#**    Author: Joey Azofeifa                                               */
#**    Date:   7/8/2015	                                                   */
#**    Description: Makefile for Exponential-Gaussian-Uniform Mixture Model*/
#***************************************************************************/

CXX             = g++
CXXFLAGS        = -Wno-unused-variable -Wno-non-virtual-dtor -std=c++11  -fopenmp -Wno-write-strings

EXEC            = ${PWD}/EMGU
ARCH  			= getconf LONG_BIT

CPP_FLAGS_32  	= -D32_BIT
CPP_FLAGS_64  	= -D64_BIT

GCCVERSION 		= $(shell gcc --version | grep ^gcc | sed 's/^.* //g')

NU_FIT: main.o load.o split.o model.o across_segments.o template_matching.o read_in_parameters.o
	@${CXX} ${CXXFLAGS} ${PWD}/main.o ${PWD}/load.o \
	${PWD}/split.o ${PWD}/model.o ${PWD}/across_segments.o \
	${PWD}/template_matching.o ${PWD}/read_in_parameters.o -o ${EXEC}
	@echo "========================="
	@echo Sucessfully Compiled

main.o: 
	@echo "========================="
	@${CXX} -c ${CXXFLAGS} ${PWD}/main.cpp 
	@echo Compiled main
load.o:
	@${CXX} -c ${CXXFLAGS} ${PWD}/load.cpp 
	@echo Compiled load
split.o:
	@${CXX} -c ${CXXFLAGS} ${PWD}/split.cpp 
	@echo Compiled split

model.o:
	@${CXX} -c ${CXXFLAGS} ${PWD}/model.cpp 
	@echo Compiled model
across_segments.o:
	@${CXX} -c ${CXXFLAGS} ${PWD}/across_segments.cpp 
	@echo Compiled across_segments

template_matching.o:
	@${CXX} -c ${CXXFLAGS} ${PWD}/template_matching.cpp 
	@echo Compiled template_matching
read_in_parameters.o:
	@${CXX} -c ${CXXFLAGS} ${PWD}/read_in_parameters.cpp 
	@echo Compiled read_in_parameters




clean:
	@rm -f *.o


		

