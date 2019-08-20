#!/bin/bash
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=7
#SBATCH -A viva-hpc
#SBATCH --time=28:00:00
#SBATCH --output=IraClassify27510Log.log
#SBATCH --job-name=IraLog5
 



module load matlab

matlab -nodisplay -nosplash -r "parpool('local',6);neuroClassifyIntraclass(0.1,7,5,5,10);exit;" 
