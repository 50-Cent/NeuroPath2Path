#!/bin/bash
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH -A USER_CLUSTER
#SBATCH --time=28:00:00
#SBATCH --output=OUTPUT.log
#SBATCH --job-name=YOUR_JOB
 



module load matlab

matlab -nodisplay -nosplash -r "parpool('local',5);neuroClassifyInterclass(0.1,15,5,5,9);exit;" 
