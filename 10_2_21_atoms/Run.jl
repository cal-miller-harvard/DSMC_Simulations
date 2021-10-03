using Printf

T1 = 2.0
l = 0.0025
L = 0.04
T2s = [0.7, 2.8]
flows = [1, 2, 4, 6]

mkpath("logs")

for flow in flows
    for T2 in T2s
        fname = @sprintf("run_m_%.4f_T2_%.4f_flow_.%.4fslurm", m, T2, flow)
        open(fname, "w") do f
            write(f,@sprintf("""#!/bin/bash
            #SBATCH -n 24 # Number of cores requested
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
        
            julia RunCells.jl --T1 %.4f --T2 %.4f -l %.4f -L %.4f --flow %.4f""", T1, T2, l, L, flow))
        end
        run(`sbatch $fname`)
    end
end
