using Printf

n = 1
m = 1
t = 1
T1s = [2.75]
T2s = [2.75]
ls = [0.00275]
Ls = [0.02]
meshs = [true, false]

mkpath("logs")

for mesh in meshs
    for i in 1:(n*m)
        for j in 1:t
            T1 = T1s[j]
            T2 = T2s[j]
            l = ls[1 + mod(i, n)]
            L = Ls[1 + convert(Int64, (i - 1E-8) ÷ n)]
            if mesh
                fname = @sprintf("run_gap_%.4f_len_%.4f_T1_%.4f_T2_%.4f.slurm", l, L, T1, T2)
            else
                fname = @sprintf("run_gap_%.4f_len_%.4f_T1_%.4f_T2_%.4f_no_mesh.slurm", l, L, T1, T2)
            end
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
            
                julia RunCells.jl --T1 %.4f --T2 %.4f -l %.4f -L %.4f --mesh
                """, T1, T2, l, L, mesh))
            end
            run(`sbatch $fname`)
        end
    end
end