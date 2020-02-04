#!/usr/bin/env bash
#SBATCH -p hebbe
#SBATCH -A SNIC2020-5-51 # Project
#SBATCH -J task8_ericlin # Name of the job
#SBATCH -N 1 # Use 1 node
#SBATCH -n 1 # Use 1 cores on that node
#SBATCH -t 15:00:00 # Maximum time
#SBATCH -o stdout_task8 # stdout goes this file
#SBATCH -e stderr_task8 # stderr goes to this file

module purge
module load intel/2019a GPAW ASE

echo ------------------ Task 8 start ------------------
mpirun -np 1 gpaw-python task10.py
echo ------------------ Task 8 done -------------------
