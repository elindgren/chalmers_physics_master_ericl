#!/usr/bin/env bash
#SBATCH -p hebbe
#SBATCH -A SNIC2020-5-51 # Project
#SBATCH -J SiBSandDOS # Name of the job
#SBATCH -N 1 # Use 1 node
#SBATCH -n 20 # Use 10 cores on that node
#SBATCH -t 50:00:00 # Maximum time - expected time is 30-40 hrs
#SBATCH -o outSi # stdout goes this file
#SBATCH -e errSi # stderr goes to this file
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ericlin@chalmers.se


module purge
module load intel/2019a GPAW ASE

echo ------------------ Task 7 start ------------------
mpirun gpaw-python task7.py
echo ------------------ Task 7 done -------------------
