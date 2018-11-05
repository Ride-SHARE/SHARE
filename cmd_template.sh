#!/bin/sh
### Set the project name, your department dc by default
#PBS -P cse
### Request email when job begins and ends
#PBS -m ea
### Specify email address to use for notification.
#PBS -M $mcs172525@cse.iitd.ac.in
####
#PBS -l select=2:ncpus=4
### Specify "wallclock time" required for this job, hhh:mm:ss
#PBS -l walltime=12:00:00

#PBS -l software=
# After job starts, must goto working directory. 
# $PBS_O_WORKDIR is the directory from where the job is fired. 
echo "==============================="
echo $PBS_JOBID
cat $PBS_NODEFILE
echo "==============================="
cd ~/combined/
#job  

