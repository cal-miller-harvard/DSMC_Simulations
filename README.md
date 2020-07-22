# DSMC_Simulations
SPARTA DSMC Simulations + Monte Carlo particle simulations for simulation of cryogenic buffer gas beams. Developed for the [PolyEDM Experiment](http://www.polyedm.com/).

## SPARTA

Which command runs SPARTA (e.g. `spa_kokkos_omp`, `spa_mpi`, ...) depends on installation. Run a SPARTA script with `sparta_command < path/to/script`. See [here](https://sparta.sandia.gov/doc/Manual.html) for SPARTA documentation and [here](https://github.com/cal-miller-harvard/DSMC_Simulations/blob/master/7_8_20_opaque_mesh/flow_4.00000_gap_0.00200_len_0.02000_T1_2.00000_T2_0.50000/in.cell) for an example SPARTA script.

The Julia script [DrawCell.jl](https://github.com/cal-miller-harvard/DSMC_Simulations/blob/master/7_16_20_methanol/DrawCell.jl) includes functions for generating 2D geometry for input to SPARTA. See a [run script](https://github.com/cal-miller-harvard/DSMC_Simulations/blob/master/7_16_20_methanol/RunCells.jl) for an example of its use.

## Particle Simulations

The particle simulations are run using the Julia program [ParticleTracing.jl](https://github.com/cal-miller-harvard/DSMC_Simulations/blob/master/ParticleTracing/ParticleTracing.jl). It requires as input a list of surfaces with the `--geom` flag (generated with `dump surfs surf all 100000 data/cell.surfs id v1x v1y v2x v2y`) and a list of temperatures, densities, and velocities with the `--flow` flag (generated with `compute temp thermal/grid all He temp; compute rhov grid all He nrho massrho u v w; fix out ave/grid all 1000 100 100000 c_temp[*] c_rhov[*]; dump out grid all 100000 data/DS2FF.DAT xc yc f_out[*]`) from a SPARTA simulation. Other command line options can be found by running `julia ParticleTracing.jl --help`. A list of the final positions and velocities of the simulated particles (e.g. [here](https://github.com/cal-miller-harvard/DSMC_Simulations/blob/master/7_8_20_opaque_mesh/flow_4.00000_gap_0.00200_len_0.02000_T1_2.00000_T2_0.50000/data/particles_omega_0.00000_M_191.0_zmax_0.08609_pflip_0.10000_job_63848593.out)) will be printed (and can be redirected to a file with `> path/to/file`). Statistics on where particles collide can be output with the `--stats` and `--exitstats` flags.

## Run Scripts

See [RunCells.jl](https://github.com/cal-miller-harvard/DSMC_Simulations/blob/master/7_8_20_opaque_mesh/RunCells.jl) for an example script which automatically runs a DSMC simulation and then particle simulations on the cluster. This is the recommended way to run the code, since it simplifies running parameter scans. The script is straightforward to modify for different simulations or parameters.