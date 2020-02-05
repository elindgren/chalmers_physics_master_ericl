#!/usr/bin/env bash
#SBATCH -p hebbe
#SBATCH -A SNIC2020-5-51 # Project
#SBATCH -J task1_test # Name of the job
#SBATCH -N 1 # Use 1 node
#SBATCH -n 20 # Use 20 cores on that node
#SBATCH -t 20:00:00 # Maximum time - expected time is 15 hrs
#SBATCH -o stdout_task1_test # stdout goes this file
#SBATCH -e stderr_task1_test # stderr goes to this file

module purge
module load intel/2019a GPAW ASE

echo ------------------ Task 1 start ------------------
mpirun gpaw-python task1.py
echo ------------------ Task 1 done -------------------
