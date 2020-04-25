import os


# Parameters:
# 4 sccm
# 5, 15, 25 mm gap
# 10 to 100 mm second stage length

flows = [4] # sccm
gaps = [5E-3, 15E-3, 25E-3] # m
lengths = [10E-3, 20E-3, 30E-3, 40E-3, 50E-3, 60E-3, 70E-3, 80E-3, 90E-3, 100E-2] # m
alphas = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900]

nparts = 100000

for flow in flows:
    for gap in gaps:
        for length in lengths:
            directory = "/n/home03/calmiller/DSMC_Simulations/4_24_20_long_second_stage/flow_{:.3f}_gap_{:.3f}_len_{:.3f}/data".format(flow, gap, length)
            print("starting {}".format(directory))
            os.chdir(directory)
            for alpha in alphas:
                with open(r"particles.slurm", "w") as f:
                    f.write("""#!/bin/bash
#SBATCH -n 4 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 2:00:00 # Runtime in minutes
#SBATCH -p shared # Partition to submit to
#SBATCH --mem-per-cpu 2048 # Memory per cpu in MB (see also ?mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o particles_alpha_{}_job_%j.out # Standard out goes to this file
#SBATCH -e particles_alpha_{}_job_%j.err # Standard err goes to this filehostname

module load intel/19.0.5-fasrc01 openmpi/4.0.2-fasrc01 fftw/3.3.8-fasrc01 cmake/3.12.1-fasrc01 Anaconda3/2019.10 python/3.7.7-fasrc01
module list

export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export JULIA_NUM_THREADS=threads

echo "running...."

/n/home03/calmiller/programs/julia /n/home03/calmiller/DSMC_Simulations/ParticleTracing/ParticleTracing.jl -z 0.035 -T 2.0 -n {} ./cell.510001.surfs ./DS2FF.500000.DAT --alpha {}

""".format(nparts, alpha, alpha, alpha))
            os.system("sbatch particles.slurm")
            os.chdir("/n/home03/calmiller/DSMC_Simulations/4_24_20_long_second_stage")
