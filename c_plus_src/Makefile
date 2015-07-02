#***************************************************************************/
#**    Author: Joey Azofeifa                                               */
#**    Date:   1/2/2015	                                                   */
#**    Description: Makefile for FastReadStitcher                          */
#***************************************************************************/

CXX             = g++-4.9.0
CXXFLAGS        = -Wno-unused-variable -Wno-non-virtual-dtor -std=c++11  -fopenmp -Wno-write-strings

EXEC            = ${PWD}/EMGU
ARCH  			= getconf LONG_BIT

CPP_FLAGS_32  	= -D32_BIT
CPP_FLAGS_64  	= -D64_BIT

GCCVERSION 		= $(shell gcc --version | grep ^gcc | sed 's/^.* //g')

NU_FIT: main.o read_in_parameters.o split.o load.o pdfs.o run_fit.o write.o
	@${CXX} ${CXXFLAGS} ${PWD}/main.o  ${PWD}/read_in_parameters.o ${PWD}/split.o ${PWD}/load.o ${PWD}/pdfs.o ${PWD}/run_fit.o ${PWD}/write.o -o ${EXEC}
	@echo "========================================="
	@echo Sucessfully Compiled

main.o: 
		@echo "========================="
		@${CXX} -c ${CXXFLAGS} ${PWD}/main.cpp 
		@echo Compiled main

split.o: 
		@@${CXX} -c ${CXXFLAGS} ${PWD}/split.cpp
		@echo Compiled split
read_in_parameters.o:
		@@${CXX} -c ${CXXFLAGS} ${PWD}/read_in_parameters.cpp
		@echo Compiled read_in_parameters
load.o:
		@@${CXX} -c ${CXXFLAGS} ${PWD}/load.cpp
		@echo Compiled load
pdfs.o:
		@@${CXX} -c ${CXXFLAGS} ${PWD}/pdfs.cpp
		@echo Compiled pdfs
run_fit.o:
		@@${CXX} -c ${CXXFLAGS} ${PWD}/run_fit.cpp
		@echo Compiled run_fit
write.o:
		@${CXX} -c ${CXXFLAGS} ${PWD}/write.cpp 
		@echo Compiled write




clean:
	@rm -f *.o


		

