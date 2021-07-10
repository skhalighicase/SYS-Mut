#!/bin/bash
#SBATCH -J Job_
#SBATCH --time=100:00:00
#SBATCH --nodes=1
#SBATCH -n 8
#SBATCH --output=Fail.out
#SBATCH --mem=16gb
#SBATCH -A vxv89
#SBATCH --mail-user=sxk1241@case.edu
#SBATCH --mail-type=FAIL
#SBATCH -o serial-R.out%j # capture jobid in output file name

module load openmpi/2.0.1
module load gcc
module load R/3.5.0
module load JAGS/4.3.0

cd 
R CMD BATCH 

