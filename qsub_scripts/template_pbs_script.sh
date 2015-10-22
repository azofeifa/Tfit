###Name the job
#PBS -N counting_meta
### Specify the number of nodes/cores
#PBS -l nodes=1:ppn=1

### Allocate the amount of memory needed
#PBS -l mem=10gb

### Set your expected walltime
#PBS -l walltime=12:00:00

### Setting to mail when the job is complete
#PBS -e /Users/azofeifa/qsub_errors/                                                                                              
#PBS -o /Users/azofeifa/qsub_stdo/  

### Set your email address
#PBS -m ae
#PBS -M jgazofeifa@gmail.com



### Choose your shell 
#PBS -S /bin/sh
### Pass enviroment variables to the job
#PBS -V

### now call your program

src=/path/to/your/script.py
param1=/input/file/path/ #taken as input to script.py
param2=5 #taken as input to script.py
###...
paramN=1000 #taken as input to script.py

python $src $param1 $param2 $paramN


