#!/usr/bin/env bash
#SBATCH -p hebbe
#SBATCH -A SNIC2020-5-51 # Project
#SBATCH -J NaGA_ericlin # Name of the job
#SBATCH -N 1 # Use 1 node
#SBATCH -n 1 # Use only 1 core on that node
#SBATCH -t 10:00:00 # Maximum time
#SBATCH -o stdout_GA # stdout goes this file
#SBATCH -e stderr_GA # stderr goes to this file

module purge
module load intel/2019a GPAW ASE

cd Na8
echo ------------------ Na8 Start ------------------
mpirun -np 1 gpaw-python ../../Na-clusters-GA-search/ga.py
echo ------------------ Na8 Done -------------------
