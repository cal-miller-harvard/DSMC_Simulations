import os


# Parameters:
# 4 sccm
# 2, 4, 6, 8, 10 mm gap
# 10, 20, 30 mm second stage length

flows = [0.5,1,2,4,8,16,32,64,128,200] # sccm
omegas = [0]
nparts = 100000
gap = 0
length = 0

for flow in flows:
    directory = "flow_{:.3f}_gap_{:.3f}_len_{:.3f}".format(flow, gap, length)
    os.system("cp -rf template "+directory)
    os.chdir(directory)
    nsim = 5E11*flow
    for omega in omegas:
        with open(r"particles_{}.slurm".format(omega), "w") as f:
            f.write("""#!/bin/bash
#SBATCH -n 16 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 2:00:00 # Runtime in minutes
#SBATCH -p shared # Partition to submit to
#SBATCH --mem-per-cpu 2048 # Memory per cpu in MB (see also ?mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o particles_omega_{}_job_%j.out # Standard out goes to this file
#SBATCH -e particles_omega_{}_job_%j.err # Standard err goes to this filehostname

module load intel/19.0.5-fasrc01 openmpi/4.0.2-fasrc01 fftw/3.3.8-fasrc01 cmake/3.12.1-fasrc01 Anaconda3/2019.10 python/3.7.7-fasrc01
module list

export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export JULIA_NUM_THREADS=threads

cd data
pwd
echo "running...."

/n/home03/calmiller/programs/julia /n/home03/calmiller/DSMC_Simulations/4_8_20_Yuiki_cell/ParticleTracingYuiki.jl -z 0.049 -T 4.0 -n {} ./cell.510001.surfs ./DS2FF.500000.DAT --omega {} --stats ./stats_omega_{}.csv --exitstats ./exitstats_omega_{}.csv

""".format(omega, omega, nparts, omega, omega, omega))
    with open(r"in.cell", "w") as f:
        f.write("""# 2d axial simulation of CBGB cell - Yuiki's T design

# Input density: 1 sccm ~ 2E21?
variable RHO equal "{:.3f} * 2E21"

global		        fnum {:.3e} cellmax 10000 gridcut 0.0 comm/sort yes
seed	    	    12345
dimension   	    2
units               si
timestep 	        1E-6

species		        he4.species He4
mixture		        He4 He4 nrho $(v_RHO) vstream 26.2 0 0 temp 4.0
collide             vss He4 He4.vss #kk

# read_restart data/cell.restart.*

boundary	        o ao p
create_box  	    0 0.12 0 0.03 -0.1 0.1
create_grid 	    20 20 1 
balance_grid        rcb cell

read_surf           data.cell
group       inlet surf id 1
fix		    in emit/surf He4 inlet

surf_collide	    diffuse diffuse 4.0 1.0
surf_modify         all collide diffuse

stats		    1000
stats_style	    wall step cpu np nattempt ncoll nscoll nscheck

# Run until convergence
label loop
variable a loop 40
run 		    10000
adapt_grid all refine particle 16 4
next a
jump in.cell loop

compute temp thermal/grid all He4 temp
compute rhov grid all He4 nrho massrho u v w
fix out ave/grid all 1000 100 100000 c_temp[*] c_rhov[*]
dump out grid all 100000 data/DS2FF.*.DAT xc yc f_out[*]

# Record statistics
label loop2
variable b loop 11
run 		    10000

next b
jump in.cell loop2

# Save surfaces
dump surfs surf all 1 data/cell.*.surfs id v1x v1y v2x v2y
run 1

""".format(flow, nsim))
    with open(r"run.slurm", "w") as f:
        f.write("""#!/bin/bash
#SBATCH -n 16 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 12:00:00 # Runtime in minutes
#SBATCH -p shared # Partition to submit to
#SBATCH --mem-per-cpu 1024 # Memory per cpu in MB (see also ?mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o spa_%j.out # Standard out goes to this file
#SBATCH -e spa_%j.err # Standard err goes to this filehostname

module load intel/19.0.5-fasrc01 openmpi/4.0.2-fasrc01 fftw/3.3.8-fasrc01 cmake/3.12.1-fasrc01 Anaconda3/2019.10 python/3.7.7-fasrc01
module list

export OMP_PROC_BIND=spread
export OMP_PLACES=threads

echo "running...."

# mpirun ~/programs/sparta/spa -kokkos off < in.cell

"""+("".join(["sbatch particles_{}.slurm\n".format(omega) for omega in omegas])))
    os.system("sbatch run.slurm")
    os.chdir("..")
