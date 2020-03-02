#!/usr/bin/env bash
#SBATCH -p hebbe
#SBATCH -A SNIC2020-5-51 # Project
#SBATCH -J EDOS_nanoparticles # Name of the job
#SBATCH -N 1 # Use 1 node
#SBATCH -n 1 # Use 1 cores on that node
#SBATCH -t 50:00:00 # Maximum time - expected time is 30-40 hrs
#SBATCH -o out_edos # stdout goes this file
#SBATCH -e err_edos # stderr goes to this file
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ericlin@chalmers.se


module purge
module load intel/2019a GPAW ASE

echo ------------------ Task 5 start ------------------
mpirun -np=1 gpaw-python task5.py
echo ------------------ Task 5 done -------------------
