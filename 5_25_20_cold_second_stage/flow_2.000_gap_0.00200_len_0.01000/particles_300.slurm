#!/bin/bash
#SBATCH -n 8 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-08:00 # Runtime in minutes
#SBATCH -p shared # Partition to submit to
#SBATCH --mem-per-cpu 2048 # Memory per cpu in MB (see also ?mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o particles_omega_300_job_%j.out # Standard out goes to this file
#SBATCH -e particles_omega_300_job_%j.err # Standard err goes to this filehostname

# module load intel/19.0.5-fasrc01 openmpi/4.0.2-fasrc01 fftw/3.3.8-fasrc01 cmake/3.12.1-fasrc01 Anaconda3/2019.10 python/3.7.7-fasrc01
# module list

# export OMP_PROC_BIND=spread
# export OMP_PLACES=threads
# export JULIA_NUM_THREADS=$SLURM_CPUS_ON_NODE

cd data
pwd
echo "running...."

julia /n/home03/calmiller/DSMC_Simulations/ParticleTracing/ParticleTracing.jl -z 0.035 -T 2.0 -n 100000 ./cell.510001.surfs ./DS2FF.500000.DAT --omega 300 --stats ./stats_omega_300.csv --exitstats ./exitstats_omega_300.csv

