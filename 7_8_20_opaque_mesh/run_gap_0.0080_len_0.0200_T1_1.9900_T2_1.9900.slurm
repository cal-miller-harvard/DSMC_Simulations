#!/bin/bash
#SBATCH -n 24 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-08:00 # Runtime in minutes
#SBATCH -p shared # Partition to submit to
#SBATCH --mem-per-cpu 1024 # Memory per cpu in MB (see also ?mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o logs/dsmc_job_%j.out # Standard out goes to this file
#SBATCH -e logs/dsmc_job_%j.err # Standard err goes to this filehostname
        
module load intel/19.0.5-fasrc01 openmpi/4.0.2-fasrc01 fftw/3.3.8-fasrc01 cmake/3.12.1-fasrc01 Anaconda3/2019.10 python/3.7.7-fasrc01
module list
        
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export JULIA_NUM_THREADS=$SLURM_CPUS_ON_NODE
        
echo "running...."
        
julia RunCells3He.jl --T1 1.9900 --T2 1.9900 -l 0.0080 -L 0.0200
