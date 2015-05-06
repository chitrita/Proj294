#!/bin/sh
### Set the job name
#PBS -N job_name_here 
### Redirect stdout and stderr by first telling Torque to redirect do /dev/null and then redirecting yourself via exec. This is the way the IT recommends.
#PBS -e localhost:/dev/null
#PBS -o localhost:/dev/null
### Set the queue to which to send
#PBS -q yosef
### Limit the resources used
#PBS -l nodes=1:ppn=2
### Change the walltime and cpu time limit from their default (the default is currently an hour)
### The format is:
### #PBS -l walltime=HH:MM:SS
### #PBS -l cput=HH:MM:SS
#PBS -l walltime=10:00:00
#PBS -l cput=10:00:00
### Move all your environment variables to the job
#PBS -V
### Set the working path to be used by the job
###PBS -d <directory_goes_here> 


### Switch to the working directory; by default Torque launches processes from your home directory.
### Jobs should only be run from /work; Torque returns results via NFS.
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

### Assign stderr and stdout
exec 2> queueErr.txt > queueLog.txt

### Run some informational commands.
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo This jobs runs on the following processors:
echo `cat $PBS_NODEFILE`

### Define number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS cpus

