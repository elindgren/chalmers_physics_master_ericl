#!/usr/bin/env bash
#SBATCH -p hebbe
#SBATCH -A SNIC2020-5-51 # Project
#SBATCH -J task1_test # Name of the job
#SBATCH -N 1 # Use 1 node
#SBATCH -n 6 # Use 20 cores on that node
#SBATCH -t 1:00:00 # Maximum time - expected time is 15 hrs
#SBATCH -o stdout_task1 # stdout goes this file
#SBATCH -e stderr_task1 # stderr goes to this file
#SBATCH --mail-type=ALL, TIME_LIMIT, TIME_LIMIT_90


module purge
module load intel/2019a GPAW ASE

echo ------------------ Task 1 start ------------------
rm mdTask1.traj
rm out_mdTask1.txt
rm log_mdTask1.txt
rm stdout_task1
rm stderr_task1

mpirun gpaw-python task1.py
echo ------------------ Task 1 done -------------------
