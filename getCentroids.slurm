#!/bin/bash
#SBATCH --job-name=getCentroids
#SBATCH --ntasks=94             
#SBATCH --mem=470gb                    
#SBATCH --time=72:00:00   
#SBATCH --partition=standard
#SBATCH --account=benquist
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mdpillet@email.arizona.edu
 
module load R
cd /xdisk/benquist/mich/Scripts/
Rscript getCentroids.R