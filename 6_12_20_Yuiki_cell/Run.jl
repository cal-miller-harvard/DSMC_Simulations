using Printf

flows = [2.0, 5.0, 10.0, 30.0, 70.0]

for flow in flows
    fname = "run_flow_$(flow).slurm"
    open(fname, "w") do f
        write(f,@sprintf("""#!/bin/bash
        #SBATCH -n 32 # Number of cores requested
        #SBATCH -N 1 # Ensure that all cores are on one machine
        #SBATCH -t 0-08:00 # Runtime in minutes
        #SBATCH -p shared # Partition to submit to
        #SBATCH --mem-per-cpu 1024 # Memory per cpu in MB (see also ?mem-per-cpu)
        #SBATCH --open-mode=append
        #SBATCH -o logs/dsmc_job_%%j.out # Standard out goes to this file
        #SBATCH -e logs/dsmc_job_%%j.err # Standard err goes to this filehostname

        module load intel/19.0.5-fasrc01 openmpi/4.0.2-fasrc01 fftw/3.3.8-fasrc01 cmake/3.12.1-fasrc01 Anaconda3/2019.10 python/3.7.7-fasrc01
        module list

        export OMP_PROC_BIND=spread
        export OMP_PLACES=threads
        export JULIA_NUM_THREADS=\$SLURM_CPUS_ON_NODE

        echo "running...."

        julia RunCells.jl --T1 %.4f --T2 %.4f -l %.4f -L %.4f --flow %.4f
        """, 4.0, 4.0, 0.0, 0.0, flow))
    end
    run(`sbatch $fname`)
end
