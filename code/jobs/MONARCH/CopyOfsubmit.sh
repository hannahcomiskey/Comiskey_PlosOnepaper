#!/bin/env bash
#SBATCH --job-name=myFirstR
#SBATCH --time=00:10:00
#SBATCH --mem=3000
#SBATCH --ntasks=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=my.username@monash.edu
#SBATCH --output=output.txt
module load  R/4.3.3
R --vanilla < job_simdata.R
# or R --vanilla < example.r > output.txt
