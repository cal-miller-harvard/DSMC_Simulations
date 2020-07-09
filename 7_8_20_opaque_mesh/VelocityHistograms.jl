using Plots
using CSV
using DataFrames
using Printf
using Statistics

function read_file(fname)
    data = DataFrame(CSV.read(fname,
        header = ["idx", "x", "y", "z", "xnext", "ynext", "znext", "vx", "vy", "vz", "collides", "time"], 
        types = [Int64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Int64, Float64],
        skipto=4, ignorerepeated=true, delim=' ', silencewarnings=true))
    dropmissing!(data)
    return data
end

function get_extracted(df)
    # zmax = maximum(df.z)
    # return df[(df.znext .> zmax) .& (df.vz .> 0) .& (df.z .> 0.12), :]
    return df
end

function make_hist!(hist, label, fname, np=1E6, bins=[i for i in 0:3:90])
    data = read_file(fname)
    exited = get_extracted(data)
    extr = length(exited.vz)/np
    stephist!(hist, v(exited), label=label*@sprintf(", extr: %.2E",extr), bins=bins, normed=false)
end

v(df) = sqrt.(df.vx.^2 .+ df.vy.^2 .+ df.vz.^2)


hist = plot(xlabel="v (m/s)", ylabel="n (per 10^6)")
make_hist!(hist, "YbOH, 0.5K, 1mm gap", "/home/cal/Documents/DSMC_Simulations/5_29_20_cold_guided/flow_4.00000_gap_0.00100_len_0.04000_T1_2.00000_T2_0.50000/data/particles_exit_omega_0.00000_M_191.0_zmax_0.10609_pflip_0.10000_job_58251244.out")
make_hist!(hist, "YbOH, 0.5K, 3mm gap", "/home/cal/Documents/DSMC_Simulations/5_29_20_cold_guided/flow_4.00000_gap_0.00300_len_0.01200_T1_2.00000_T2_2.00000/data/particles_exit_omega_0.00000_M_191.0_zmax_0.07809_pflip_0.10000_job_58355349.out")
make_hist!(hist, "YbOH, 2K, 3mm gap", "/home/cal/Documents/DSMC_Simulations/5_29_20_cold_guided/flow_4.00000_gap_0.00300_len_0.01200_T1_2.00000_T2_2.00000/data/particles_exit_omega_0.00000_M_191.0_zmax_0.07809_pflip_0.10000_job_58355349.out")

make_hist!(hist, "CaOH, 0.5K, 1mm gap", "/home/cal/Documents/DSMC_Simulations/5_29_20_cold_guided/flow_4.00000_gap_0.00100_len_0.04000_T1_2.00000_T2_0.50000/data/particles_exit_omega_0.00000_M_57.0_zmax_0.10609_pflip_0.10000_job_58370013.out")
make_hist!(hist, "CaOH, 0.5K, 3mm gap", "/home/cal/Documents/DSMC_Simulations/5_29_20_cold_guided/flow_4.00000_gap_0.00300_len_0.04000_T1_2.00000_T2_0.50000/data/particles_exit_omega_0.00000_M_57.0_zmax_0.10609_pflip_0.10000_job_58370014.out")
make_hist!(hist, "CaOH, 2K, 3mm gap", "/home/cal/Documents/DSMC_Simulations/5_29_20_cold_guided/flow_4.00000_gap_0.00300_len_0.01200_T1_2.00000_T2_2.00000/data/particles_exit_omega_0.00000_M_57.0_zmax_0.07809_pflip_0.10000_job_58355350.out")
savefig("summary/hist.pdf")
display(hist)
