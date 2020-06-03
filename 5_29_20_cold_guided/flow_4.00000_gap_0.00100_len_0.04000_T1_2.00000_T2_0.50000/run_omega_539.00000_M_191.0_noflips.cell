#!/bin/bash
#SBATCH -n 16 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-08:00 # Runtime in minutes
#SBATCH -p shared # Partition to submit to
#SBATCH --mem-per-cpu 1024 # Memory per cpu in MB
#SBATCH --open-mode=append
#SBATCH -o data/particles_noflip_omega_539.00000_M_191.0_zmax_0.10609_pflip_0.10000_job_%j.out # Standard out goes to this file
#SBATCH -e data/particles_noflip_omega_539.00000_M_191.0_zmax_0.10609_pflip_0.10000_job_%j.err # Standard err goes to this filehostname

module load intel/19.0.5-fasrc01 openmpi/4.0.2-fasrc01 fftw/3.3.8-fasrc01 cmake/3.12.1-fasrc01 Anaconda3/2019.10 python/3.7.7-fasrc01
module list

export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export JULIA_NUM_THREADS=$SLURM_CPUS_ON_NODE

cd data
pwd
echo "running...."

julia /n/home03/calmiller/DSMC_Simulations/ParticleTracing/ParticleTracing.jl -z 0.035 -T 2.0 -n 1000000 ./cell.surfs ./DS2FF.DAT --omega 539.00000 --pflip 0.00000 -m 4.00000 -M 191.00000 --sigma 1.30000E-18 --zmin -100.0 --zmax 100.0 --stats ./stats_noflip_omega_539.00000_M_191.0_zmax_0.10609_pflip_0.10000.csv --exitstats ./exitstats_noflip_omega_539.00000_M_191.0_zmax_0.10609_pflip_0.10000.csv --saveall 0
