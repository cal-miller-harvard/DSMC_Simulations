using CSV, StatsPlots, Printf, DataFrames, Plots, Statistics

program_dir = "/home/cal/Documents/DSMC_Simulations/4_8_20_Yuiki_cell"
data_dir = "/home/cal/Documents/DSMC_Simulations/4_8_20_Yuiki_cell"

"""
    read_batch()

Reads in particle simulation data for a set of simulations generated by SPARTA and ParticleTracing.jl
"""
function read_batch()
    dirs = readdir()
    dir_pattern = r"flow_(\d+\.?\d+)_gap_(\d+\.?\d+)_len_(\d+\.?\d+)"
    part_pattern = r"particles_omega_(\d+)_job_(\d+).out"

    data = DataFrame()

    for dir in dirs
        md = match(dir_pattern, dir)
        if md === nothing
            println("no match for: ", dir)
        else
            flow = parse(Float64, md.captures[1])
            @printf("flow: %.3f\n", flow)
            files = readdir(dir)
            for f in files
                md = match(part_pattern, f)
                if !(md === nothing)
                    part_file = f
                    omega = parse(Float64, md.captures[1])
                    valid = false
                    linenum = 0
                    fname = dir*"/"*part_file
                    for line in eachline(fname)
                        if occursin("idx x y z", line)
                            valid = true
                            break
                        end
                        linenum += 1
                        if linenum > 100
                            break
                        end
                    end
                    if valid
                        println(fname)
                        table = CSV.read(fname,
                        header = ["idx", "x", "y", "z", "xnext", "ynext", "znext", "vx", "vy", "vz", "collides", "time"], 
                        types = [Int64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Int64, Float64],
                        skipto=linenum+2, ignorerepeated=true,delim=' ')
                        table[!, :flow] .= flow
                        table[!, :omega] .= omega
                        append!(data, dropmissing(table))
                    else
                        println("particle data not valid.")
                    end
                end
            end
        end
    end
    return data
end


# Loading
cd(data_dir)
data = read_batch()
cd(program_dir)

# Plotting
npts = 200000
flows = sort(unique(data.flow))
extractions = zero(flows)
medians = zero(flows)
vmax = maximum(data.vz)

plotdir = data_dir*"/summary"
plotpath = mkpath(plotdir)
cd(plotdir)

plt = plot(xlabel="vz (m/s)", ylabel="probability", legend=:topleft)
bins = 0:vmax/100:vmax
for (i, f) in enumerate(flows)
    vzs = data[data.flow .== f,:].vz
    extractions[i] = length(vzs)/npts
    medians[i] = median(vzs)
    stephist!(vzs, label="$f sccm", bins=bins, normalize=:probability)
end
display(plt)
savefig("vz.pdf")

plt = plot(flows, extractions, xlabel="flow (sccm)", ylabel="extraction")
display(plt)
savefig("extraction.pdf")

plt = plot(flows, medians, xlabel="flow (sccm)", ylabel="median vz (m/s)")
display(plt)
savefig("medianvz.pdf")

cd(program_dir)
