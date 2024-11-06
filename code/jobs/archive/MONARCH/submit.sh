#!/bin/env bash
#SBATCH --job-name=FactorCovFPSource
#SBATCH --time=00:10:00
#SBATCH --mem=3000
#SBATCH --ntasks=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hannah.comiskey@monash.edu
#SBATCH --output=STAN_model_zih_deltak_factorcov_Q_obsdata.RDS
module load  R/4.3.3
R --vanilla < job_file.R
