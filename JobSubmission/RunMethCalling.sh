#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p sq # submit to the serial queue
#SBATCH --time=12:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-193495 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --error=JobSubmissionLogs/RunMethCalling.err # error file
#SBATCH --output=JobSubmissionLogs/RunMethCalling.log # output file
#SBATCH --job-name=RunMethCalling
#SBATCH -D /gpfs/ts0/projects/Research_Project-193495/ONTMethCalling/ # set working directory to

## Run from github repo
echo ${HOME}
export SCRIPTSDIR=$(pwd)/Scripts


module load Anaconda3/5.2.0
source activate ont_meth
source ./Config/SmokingPilot.txt

## align reads
sh ${SCRIPTSDIR}/Alignment.sh

## calculate sequencing statistics
module load BEDTools
sh ${SCRIPTSDIR}/CalcCoverage.sh

sh ${SCRIPTSDIR}/DNAMcallingNanopolish.sh