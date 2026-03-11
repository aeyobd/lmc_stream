#!/bin/bash
#SBATCH -A durham
#SBATCH -p cosma5
#SBATCH -t 10:0:0
#SBATCH --tasks 1
#SBATCH -n 1
#SBATCH --cpus-per-task 16
#SBATCH --mem-per-cpu 2048


set -xe 

export OMP_NUM_THREADS=16
project_potential.jl ../potential_mw_halo_evo.ini -T times.txt -k 1 --limits 400 -n 1025 -o projected_mw_halo.hdf5 
