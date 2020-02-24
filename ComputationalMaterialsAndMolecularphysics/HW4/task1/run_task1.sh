#!/usr/bin/env bash
#SBATCH -p hebbe
#SBATCH -A SNIC2020-5-51 # Project
#SBATCH -J TDDFT-1 # Name of the job
#SBATCH -N 1 # Use 1 node
#SBATCH -n 1 # Use 10 cores on that node
#SBATCH -t 50:00:00 # Maximum time - expected time is 30-40 hrs
#SBATCH -o stdout # stdout goes this file
#SBATCH -e stderr # stderr goes to this file
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ericlin@chalmers.se


module purge
module load intel/2019a GPAW ASE

echo ------------------ Task 1 start ------------------
mpirun gpaw-python converge.py &&
mpirun gpaw-python spectrum.py
echo ------------------ Task 1 done -------------------
