#!/usr/bin/env bash
#SBATCH -p hebbe
#SBATCH -A SNIC2020-5-51 # Project
#SBATCH -J task8_ericlin # Name of the job
#SBATCH -N 1 # Use 1 node
#SBATCH -n 10 # Use 10 cores on that node
#SBATCH -t 15:00:00 # Maximum time
#SBATCH -o stdout_task8 # stdout goes this file
#SBATCH -e stderr_task8 # stderr goes to this file

module purge
# module load intel/2019a GPAW ASE
module spider  GPAW/19.8.1-Python-3.7.2-libvdwxc
module load GCC/8.2.0-2.31.1  OpenMPI/3.1.3
module load  GPAW/19.8.1-Python-3.7.2-libvdwxc

echo ------------------ Task 8 start ------------------
mpirun -np 1 gpaw-python task8.py
echo ------------------ Task 8 done -------------------
