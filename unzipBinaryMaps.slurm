#!/bin/bash
#SBATCH --job-name=unzipBinaryMaps
#SBATCH --ntasks=1             
#SBATCH --mem=4gb                    
#SBATCH --time=06:00:00   
#SBATCH --partition=standard
#SBATCH --account=benquist
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mdpillet@email.arizona.edu
 
cd /xdisk/benquist/mich/
 
./unzipBinaryMaps.sh