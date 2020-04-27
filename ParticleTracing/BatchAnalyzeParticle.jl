using CSV, StatsPlots, Printf, DataFrames, Plots

program_dir = "/home/cal/Documents/DSMC_Simulations/ParticleTracing"
data_dir = "/home/cal/Documents/DSMC_Simulations/4_24_20_long_second_stage"


"""
    read_batch()

Reads in particle simulation data for a set of simulations generated by SPARTA and ParticleTracing.jl
"""
function read_batch()
    dirs = readdir()
    dir_pattern = r"flow_(\d+\.?\d+)_gap_(\d+\.?\d+)_len_(\d+\.?\d+)"
    part_pattern = r"particles_alpha_(\d+)_job_(\d+).out"
    cell_pattern = r"cell.(\d+).surfs"

    data = DataFrame()

    for dir in dirs
        md = match(dir_pattern, dir)
        if md === nothing
            println("no match for: ", dir)
        else
            flow = parse(Float64, md.captures[1])
            gap = parse(Float64, md.captures[2])
            len = parse(Float64, md.captures[3])
            @printf("flow: %.3f gap: %.3f length: %.3f\n", flow, gap, len)
            subdirs = readdir(dir)
            if "data" in subdirs
                files = readdir(dir*"/data")
                part_file = nothing
                geom_exists = false
                max_x = 10.0E-3
                for f in files
                    mc = match(cell_pattern, f)
                    if !(mc === nothing)
                        fname = dir*"/data/"*f
                        box = Matrix(CSV.read(fname, header = ["min","max"], skipto=6, limit=2,ignorerepeated=true,delim=' ',silencewarnings=true))
                        max_x = box[1,2]
                        geom_exists = true
                    end
                end
                for f in files
                    md = match(part_pattern, f)
                    if !(md === nothing) && geom_exists
                        part_file = f
                        omega = parse(Float64, md.captures[1])
                        valid = false
                        linenum = 0
                        fname = dir*"/data/"*part_file
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
                            table = CSV.read(fname,
                            header = ["idx", "x", "y", "z", "xnext", "ynext", "znext", "vx", "vy", "vz", "collides", "time"], 
                            types = [Int64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Int64, Float64],
                            skipto=linenum+2, ignorerepeated=true,delim=' ')
                            table[!, :flow] .= flow
                            table[!, :gap] .= gap
                            table[!, :len] .= len
                            table[!, :omega] .= omega
                            table[!, :maxx] .= max_x
                            append!(data, table)
                        else
                            println("particle data not valid.")
                        end
                    end
                end
            else
                println("data subdirectory not found.")
            end
        end
    end
    return data
end


# Loading
# data_dir = "/home/cal/Documents/DSMC_Simulations/4_3_20_flow_gap_length"
cd(data_dir)
data = read_batch()
cd(program_dir)

# Analysis
lens = unique(data.len)
lens = [i for i in lens if i < 0.1] 
flows = unique(data.flow)
gaps = unique(data.gap)
omegas = unique(data.omega)
omegas = [i for i in omegas if i < 650] 
extracted = zeros(length(lens), length(flows), length(gaps), length(omegas))
slow = zeros(length(lens), length(flows), length(gaps), length(omegas))

for (i, len) in enumerate(lens)
    for (j, flow) in enumerate(flows)
        for (k, gap) in enumerate(gaps)
            for (l, omega) in enumerate(omegas)
                part = data[(data.len .== len) .& (data.flow .== flow) .& (data.gap .== gap) .& (data.omega .== omega), :]
                extr = part[part.znext .> part.maxx, :]
                extracted[i,j,k,l] = nrow(extr)/nrow(part)
                if(isnan(extracted[i,j,k,l]))
                    extracted[i,j,k,l] = 0
                end
                slowed = extr[(extr.vz .< 15.0) .& (extr.vz .> 0.0), :]
                slow[i,j,k,l] = nrow(slowed)/nrow(part)
                if(isnan(slow[i,j,k,l]))
                    slow[i,j,k,l] = 0
                end
            end
        end
    end
end

plotdir = data_dir*"/summary"
plotpath = mkpath(plotdir)
cd(plotdir)

colors = cgrad(:inferno, scale = :log)

extractionplots = Array{Plots.Plot}(UndefInitializer(), (length(flows),length(gaps)))
for (i, flow) in enumerate(flows)
    for (j, gap) in enumerate(gaps)
        extractionplots[i,j] = heatmap(omegas, lens*1000, extracted[:,i,j,:], xlabel="omega (1/s)", ylabel="Second stage length (mm)", title=@sprintf("Extraction for flow of %.3f sccm and gap of %.3f mm", flow, 1000*gap),clims=(minimum(extracted),maximum(extracted)),c = colors)
        savefig(@sprintf("extraction_flow_%.3f_gap_%.3f.pdf", flow, gap))
    end
end

slowplots = Array{Plots.Plot}(UndefInitializer(), (length(flows),length(gaps)))
for (i, flow) in enumerate(flows)
    for (j, gap) in enumerate(gaps)
        slowplots[i,j] = heatmap(omegas, lens*1000, slow[:,i,j,:], xlabel="omega (1/s)", ylabel="Second stage length (mm)", title=@sprintf("Fraction < 15 m/s for flow of %.3f sccm and gap of %.3f mm", flow, 1000*gap),clims=(minimum(slow),maximum(slow)),c = colors)
        savefig(@sprintf("slowed_flow_%.3f_gap_%.3f.pdf", flow, gap))
    end
end

extr_0 = data[(data.znext .> data.maxx) .& (data.omega .== 0), :]
extr_300 = data[(data.znext .> data.maxx) .& (data.omega .== 300), :]
gaps = unique(data.gap)
vzbins = [0:maximum(data.vz)/20:maximum(data.vz)+0.001;]
lenbins = [0:maximum(data.len)/length(unique(data.len)):maximum(data.len)+0.001;]
colors = cgrad(:inferno)
sx = 600
sy = 600

histograms = Array{Plots.Plot}(UndefInitializer(), length(gaps))
for (i, gap) in enumerate(gaps)
    d0 = extr_0[extr_0.gap .== gap, :]
    d300 = extr_300[extr_300.gap .== gap, :]
    histograms[i] = plot(
        histogram2d(d0.vz, d0.len, clims=(0, 250), colorbar=:bottom, bins = (vzbins, lenbins), c=colors,size=(sx,sy),
        title=@sprintf("omega: 0, gap: %.0f mm", 1000*gap)),
        histogram2d(d300.vz, d300.len, clims=(0, 250), colorbar=:bottom, bins = (vzbins, lenbins), c=colors,size=(sx,sy),
        title=@sprintf("omega: 300, gap: %.0f mm", 1000*gap)),
        layout=@layout grid(2,1)
    )
    savefig(@sprintf("hist_gap_%.3f.pdf", gap))
end
# plot(histograms..., layout=@layout grid(1, length(histograms)))
# savefig("vz_len_histogram.pdf")

cd(program_dir)
