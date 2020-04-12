import os


# Parameters:
# 0.5 to 10 sccm
# 0 to 4 mm gap
# 3/8" to 1/2" second stage length

nparts = 100000

flows = [0.5, 1, 2, 3, 4, 6, 8, 10] # sccm
gaps = [0, 1E-3, 2E-3, 3E-3, 4E-3] # m
lengths =[9.525E-3, 12.7E-3] # m
for flow in flows:
    for gap in gaps:
        for length in lengths:
            directory = "flow_{:.3f}_gap_{:.3f}_len_{:.3f}/data".format(flow, gap, length)
            print("starting {}".format(directory))
            os.chdir(directory)
            with open(r"particles.slurm", "w") as f:
                f.write("""#!/bin/bash
#SBATCH -n 4 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 2:00:00 # Runtime in minutes
#SBATCH -p shared # Partition to submit to
#SBATCH --mem-per-cpu 2048 # Memory per cpu in MB (see also ?mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o particles_%j.out # Standard out goes to this file
#SBATCH -e particles_%j.err # Standard err goes to this filehostname

module load intel/19.0.5-fasrc01 openmpi/4.0.2-fasrc01 fftw/3.3.8-fasrc01 cmake/3.12.1-fasrc01 Anaconda3/2019.10 python/3.7.7-fasrc01
module list

export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export JULIA_NUM_THREADS=threads

echo "running...."

/n/home03/calmiller/programs/julia /n/home03/calmiller/DSMC_Simulations/ParticleTracing/ParticleTracing.jl -n {} ./cell.510001.surfs ./DS2FF.500000.DAT

""".format(nparts))
            os.system("sbatch particles.slurm")
            os.chdir("/n/home03/calmiller/DSMC_Simulations/4_3_20_flow_gap_length")
