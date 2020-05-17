using Distributed
using Random
using LinearAlgebra
using NearestNeighbors
using CSV
using DataFrames
using Printf
using ArgParse
using OnlineStats

import Base: convert

# Constants
const MASS_PARTICLE = 174 + 16 + 1    # AMU
const MASS_BUFFER_GAS = 4             # AMU
const MASS_REDUCED = MASS_PARTICLE * MASS_BUFFER_GAS /
    (MASS_PARTICLE + MASS_BUFFER_GAS) # AMU
const kB = 8314.46                    # AMU m^2 / (s^2 K)
const σ_BUFFER_GAS_PARTICLE = 100E-20 # m^2

# Data structures and functions for tracking trajectory statistics
struct TrajStats
    v::OnlineStat
    t::OnlineStat
end

struct StatsArray
    stats::Matrix{TrajStats}
    minr::Float64
    maxr::Float64
    rbins::Int
    minz::Float64
    maxz::Float64
    zbins::Int
    rstep::Float64
    zstep::Float64
end

@inline TrajStats() = TrajStats(CovMatrix(2), Variance())

@inline function StatsArray(minr, maxr, rbins, minz, maxz, zbins)
    stats = Array{TrajStats}(undef, rbins, zbins)
    for i in 1:rbins
        for j in 1:zbins
            stats[i,j] = TrajStats()
        end
    end
    return StatsArray(stats, minr, maxr, rbins, minz, maxz, zbins, rbins/(maxr-minr), zbins/(maxz-minz))
end

@inline function updateStats!(s::StatsArray, x, v, t)
    r = sqrt(x[1]^2+x[2]^2)
    ridx = min(s.rbins, 1+floor(Int, s.rstep*(r-s.minr)))
    zidx = min(s.zbins, 1+floor(Int, s.zstep*(x[3]-s.minz)))
    fit!(s.stats[ridx, zidx].v, [(-x[2]*v[1]+x[1]*v[2])/sqrt(x[1]^2+x[2]^2),v[3]])
    fit!(s.stats[ridx, zidx].t, t)
    return s.stats[ridx, zidx]
end

@inline function merge!(a::StatsArray, b::StatsArray)
    for i in 1:a.rbins
        for j in 1:a.zbins
            OnlineStats.merge!(a.stats[i,j].v, b.stats[i,j].v)
            OnlineStats.merge!(a.stats[i,j].t, b.stats[i,j].t)
        end
    end
    return a
end

@inline function convert(::Type{Matrix}, s::StatsArray)
    # r, z, n, t, tvar, vr, vz, vrcov, vzcov, vrvzcov
    M = Array{Float64}(undef, s.rbins*s.zbins, 10)
    idx = 1
    for i in 1:s.rbins
        for j in 1:s.zbins
            stats = s.stats[i,j]
            M[idx,1] = s.minr+(i-0.5)/s.rstep
            M[idx,2] = s.minz+(j-0.5)/s.zstep
            M[idx,3] = stats.t.n
            M[idx,4] = stats.t.μ
            M[idx,5] = stats.t.σ2
            M[idx,6] = stats.v.b[1]
            M[idx,7] = stats.v.b[2]
            M[idx,8] = stats.v.A[1,1]
            M[idx,9] = stats.v.A[2,2]
            M[idx,10] = stats.v.A[1,2]
            idx += 1
        end
    end
    return M
end

@inline function convert(::Type{DataFrame}, s::StatsArray)
    m = convert(Matrix, s)
    df = DataFrame(m, [:r, :z, :n, :t, :tvar, :vr, :vz, :vrvar, :vzvar, :vrvzcov])
    return df
end

"""
    collide!(v, vgx, vgy, vgz, T)

Accepts as input the velocity of a particle v, the mean velocity of a buffer gas atom vgx, vgy, vgz, and the buffer gas temperature T. Computes the velocity of the particle after they undergo a collision, treating the particles as hard spheres and assuming a random scattering parameter and buffer gas atom velocity (assuming the particles are moving slower than the buffer gas atoms). Follows Appendices B and C of Boyd 2017. Note that v and vg are modified.
"""
@inline function collide!(v::Vector, vgx::Number, vgy::Number, vgz::Number, T::Number)
    vg = sqrt(abs(-2 * kB * T / MASS_BUFFER_GAS * log(1-Random.rand())))
    θv = π * Random.rand()
    φv = 2 * π * Random.rand()
    vgx += Random.randn() * vg * sin(θv) * cos(φv)
    vgy += Random.randn() * vg * sin(θv) * sin(φv)
    vgz += Random.randn() * vg * cos(θv)
    cosχ = 2*Random.rand() - 1
    sinχ = sqrt(1 - cosχ^2)
    θ = 2 * π * Random.rand()
    g = sqrt((v[1] - vgx)^2 + (v[2] - vgy)^2 + (v[3] - vgz)^2)
    v[1] = MASS_PARTICLE * v[1] + MASS_BUFFER_GAS * (vgx + g * cosχ)
    v[2] = MASS_PARTICLE * v[2] + MASS_BUFFER_GAS * (vgy + g * sinχ * cos(θ))
    v[3] = MASS_PARTICLE * v[3] + MASS_BUFFER_GAS * (vgz + g * sinχ * sin(θ))
    v .= v ./ (MASS_PARTICLE + MASS_BUFFER_GAS)
end

"""
    freePropagate!(xnext, x, v, d, ω)

Updates xnext by propagating a particle at x with velocity v a distance d in a harmonic potential with frequency ω
"""
@inline function freePropagate!(xnext::Vector, x::Vector, v::Vector, d::Number, ω::Number)
    if ω != 0
        t = min(d/LinearAlgebra.norm(v), 1.0)
        sint = sin(sqrt(2)*ω*t)
        cost = cos(sqrt(2)*ω*t)
        xnext[1] = x[1]*cost + v[1]*sint/(sqrt(2)*ω)
        xnext[2] = x[2]*cost + v[2]*sint/(sqrt(2)*ω)
        xnext[3] = x[3] + v[3]*t
        v[1] = v[1]*cost-2*x[1]*ω*sint
        v[2] = v[2]*cost-2*x[2]*ω*sint
    else
        xnext .= x .+ d .* LinearAlgebra.normalize(v)
    end
end

"""
    freePath(vrel, T, ρ)

Accepts as input the velocity of a particle relative to a buffer gas atom v, the buffer gas temperature T and the buffer gas density ρ. Draws a distance the particle travels before it hits a buffer gas atom from an exponential distribution, accounting for a velocity-dependent mean free path. Note that this assumes that the gas properties don't change significantly over a mean free path.
"""
@inline function freePath(vrel::Number, T::Number, ρ::Number)
    return -log(Random.rand()) * sqrt(vrel^2/(3*kB*T/MASS_BUFFER_GAS + vrel^2))/(ρ*σ_BUFFER_GAS_PARTICLE)
end


"""
    getIntersection(x1, x2, y2, x3, y3, x4, y4)

Returns whether the line segments ((x1, y1), (x2, y2)) and ((x3, y3), (x4, y4)) intersect. See "Faster Line Segment Intersection" from Graphics Gems III ed. David Kirk.
"""
@inline function getIntersection(x1::Number, y1::Number, x2::Number, y2::Number, x3::Number, y3::Number, x4::Number, y4::Number)
    denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    num = x4*(y1 - y3) + x1*(y3-y4) + x3*(y4-y1)
    if denom > 0
        if num < 0 || num > denom
            return false
        end
        num = x3*(y2 - y1) + x2*(y1 - y3) + x1*(y3-y2)
        if num < 0 || num > denom
            return false
        end
    else
        if num > 0 || num < denom
            return false
        end
        num = x3*(y2 - y1) + x2*(y1 - y3) + x1*(y3-y2)
        if num > 0 || num < denom
            return false
        end
    end
    return true
end


"""
    propagate(xinit, vin, interp!, getCollision, omega=0.0)

Accepts as input the position of a particle xinit, its velocity v, the function interp!, which takes as input a position and a vector and updates the vector to describe the gas [x, y, vgx, vgy, vgz, T, ρ], and the function getCollision(x1, x2), which returns whether if the segment from x1 to x2 intersects geometry. Computes the path of the particle until it getCollision returns true. Returns a vector of simulation results with elements x, y, z, xnext, ynext, znext, vx, vy, vz, collides, time.
"""
@inline function propagate(xinit::Vector, vin::Vector, interp!::Function, getCollision::Function, ω=0.0, stats=nothing)
    x = deepcopy(xinit)
    props = zeros(8) # x, y, vgx, vgy, vgz, T, ρ, dmin
    interp!(props, x)
    xnext = deepcopy(x)
    v = deepcopy(vin)
    time = 0.0
    collides = 0
    
    if LinearAlgebra.norm(v) < 1E-6
        collide!(v, props[3], props[4], props[5], props[6])
        collides += 1
    end

    while true
        interp!(props, x)
        vrel = sqrt((v[1] - props[3])^2 + (v[2] - props[4])^2 + (v[3] - props[5])^2)
        dist = freePath(vrel, props[6], props[7])
        freePropagate!(xnext, x, v, dist, ω)
        if getCollision(x, xnext) != 0
            xnext .= min.(max.(xnext, -1000),1000)
            return (x[1], x[2], x[3], xnext[1], xnext[2], xnext[3], v[1], v[2], v[3], collides, time)
        else
            time += dist / LinearAlgebra.norm(v)
            collides += 1
        end
        x .= xnext
        collide!(v, props[3], props[4], props[5], props[6])
        if !isnothing(stats)
            updateStats!(stats, x, v, time)
        end
    end
end

interps = 0
"""
    SimulateParticles(geomFile, gridFile, nParticles, generateParticle)

    Runs nParticles simulations using the SPARTA outputs with paths geomFile and gridFile to define the surfaces and buffer gas properties. Generates each particle with position and velocity returned by generateParticle. Returns a nParticles by 11 matrix of simulation results, with columns x, y, z, xnext, ynext, znext, vx, vy, vz, collides, time.
"""
function SimulateParticles(
    geomFile::AbstractString, 
    gridFile::AbstractString, 
    nParticles::Integer,
    generateParticle::Function,
    print_stuff=true,
    ω=0.0,
    saveall=0,
    savestats=nothing,
    saveexitstats=nothing,
    rbins=100,
    zbins=100
    )

    bounds = Matrix(CSV.read(geomFile, header = ["min","max"], skipto=6, limit=2,ignorerepeated=true,delim=' '))
    geom = Matrix(CSV.read(geomFile, header=["ID","x1","y1","x2","y2"],skipto=10, ignorerepeated=true,delim=' '))
    griddf = CSV.read(gridFile, header=["x","y","T","ρ","ρm","vx","vy","vz"],skipto=10,ignorerepeated=true,delim=' ')

    # Reorder columns and only include grid cells with data
    DataFrames.select!(griddf, [:x, :y, :vx, :vy, :vz, :T, :ρ])
    grids = griddf[griddf.T .> 0, :]
    griddf[!, :dmin] .= 0
    grids = Matrix(griddf)

    # Make a tree for efficient nearest neighbor search
    kdtree = KDTree(transpose(Matrix(grids[:,1:2])); leafsize=10)

    # For each point, find the 100 nearest neighbors and compute the distance of the closest point at which any of the parameters varies by more than 10%. Add that distance as a column to grid.

    err = 0.2
    for i in 1:size(grids)[1]
        idxs, dists = knn(kdtree, grids[i,[1,2]], 100, true)
        for (j, idx) in enumerate(idxs)
            for k in [3,4,5,6,7]
                if !(err*grids[i,k] < grids[idx,k] < (1+err)*grids[i,k])
                    grids[i, 8] = dists[j]
                    break
                end
            end
        end
    end

    @printf(stderr, "Mean dmin: %f\n", sum(grids[:,8]/size(grids)[1]))

    """
        interpolate!(props, x)
    
    Updates the gas properties props with the data from point x.
    """
    @inline function interpolate!(props::Vector, x::Vector)
        # x, y, vgx, vgy, vgz, T, ρ, dmin
        if sqrt((x[3] - props[1])^2 + (sqrt(x[1]^2 + x[2]^2) - props[2])^2) > props[8]
            global interps
            interps += 1
            interp = view(grids, knn(kdtree, [x[3], sqrt(x[1]^2 + x[2]^2)], 1)[1][1], :)
            props[1] = interp[1]            # x 
            props[2] = interp[2]            # y
            θ = atan(x[2],x[1])
            props[3] = interp[4] * cos(θ)   # vgx
            props[4] = interp[4] * sin(θ)   # vgy
            props[5] = interp[3]            # vgz
            props[6] = interp[6]            # T
            props[7] = interp[7]            # ρ
            props[8] = interp[8]
        end
    end

    """
        getCollision(x1, x2)

    Checks whether the line between x1 and x2 intersects geometry or the boundary of the simulation region.
    """
    @inline function getCollision(x1::Vector, x2::Vector)
        r1 = sqrt(x1[1]^2 + x1[2]^2)
        r2 = sqrt(x2[1]^2 + x2[2]^2)
        if x2[3] < bounds[1,1] || x2[3] > bounds[1,2] || r2 > bounds[2,2]
            return 2
        end
        for i in 1:size(geom)[1]
            if getIntersection(geom[i,2],geom[i,3],geom[i,4],geom[i,5],x1[3],r1,x2[3],r2)
                return 1
            end
        end
        return 0
    end

    output_dim = length(propagate(zeros(3), zeros(3), interpolate!, (x,y)->true))
    outputs = zeros(nParticles, output_dim)
    if !isnothing(savestats)
        allstats = StatsArray(bounds[2,1], bounds[2,2], rbins, bounds[1,1], bounds[1,2], zbins)
    end
    if !isnothing(saveexitstats)
        boundstats = StatsArray(bounds[2,1], bounds[2,2], rbins, bounds[1,1], bounds[1,2], zbins)

    end
    Threads.@threads for i in 1:nParticles
        stats = StatsArray(bounds[2,1], bounds[2,2], rbins, bounds[1,1], bounds[1,2], zbins)
        xpart, vpart = generateParticle()
        outputs[i,:] .= propagate(xpart, vpart, interpolate!, getCollision, ω, stats)
        colltype = getCollision(outputs[i,[1,2,3]], outputs[i,[4,5,6]])
        if !isnothing(savestats)
            merge!(allstats, stats)
        end
        if colltype == 2 && !isnothing(saveexitstats)
            merge!(boundstats, stats)
        end
        if print_stuff && (saveall != 0 || colltype == 2)
            println(@sprintf("%d %e %e %e %e %e %e %e %e %e %d %e", i, 
            outputs[i,1], outputs[i,2], outputs[i,3], outputs[i,4], outputs[i,5], outputs[i,6], outputs[i,7], outputs[i,8], outputs[i,9], outputs[i,10], outputs[i,11]))
        end
    end
    
    return outputs, boundstats, allstats
end

"""
    parse_commandline()

Parses command-line arguments.
"""
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "geom"
            help = "the file generated by SPARTA containing a list of surfaces"
            required = true
        "flow"
            help = "the file generated by SPARTA containing a dump of flow properties at each grid point"
            required = true
        "-n"
            help = "the number of particles to simulate."
            arg_type = Int
            default = 10000
        "-z"
            help = "axial position at which to spawn particles (m)"
            arg_type = Float64
            default = 0.035
        "-r"
            help = "radial position at which to spawn particles (m)"
            arg_type = Float64
            default = 0.0
        "--vz"
            help = "mean axial velocity at which to spawn particles (m/s)"
            arg_type = Float64
            default = 0.0
        "--vr"
            help = "mean radial velocity at which to spawn particles (m/s)"
            arg_type = Float64
            default = 0.0
        "-T"
            help = "temperature at which to spawn particles (K)"
            arg_type = Float64
            default = 0.0
        "-m"
            help = "mass of buffer gas atoms (AMU); not currently supported"
            arg_type = Float64
            default = 4.0
        "-M"
            help = "mass of particles (AMU); not currently supported"
            arg_type = Float64
            default = 191.0
        "--sigma"
            help = "collision cross section between buffer gas and particle (m^2); not currently supported"
            arg_type = Float64
            default = 100E-20
        "--omega"
            help = "sqrt(v/m) for a harmonic trap of V(x,y,z) = v*(x^2+y^2)"
            arg_type = Float64
            default = 0.0
        "--saveall"
            help = "save all particles or just those that leave the cell"
            arg_type = Int
            default = 0
        "--stats"
            help = "filename to save average properties of trajectories"
            arg_type = String
            default = nothing
        "--exitstats"
            help = "filename to save average properties of trajectories where the particle hits the simulation boundary"
            arg_type = String
            default = nothing
    end

    return parse_args(s)
end
"""
    main()

The main function starts a particle simulation based on the command line arguments, printing outputs to stdout and timing information to stderr.
"""
function main()
    args = parse_commandline()

    nParticles = args["n"]

    println("idx x y z xnext ynext znext vx vy vz collides time")

    # Define particle generation
    boltzmann = sqrt(kB*args["T"]/MASS_PARTICLE)
    generateParticle() = (
        [args["r"], 0.0, args["z"]],
        [args["vr"] + Random.randn() * boltzmann, Random.randn() * boltzmann, args["vz"] + Random.randn() * boltzmann])

    # Set simulation parameters and run simulation
    nthreads = Threads.nthreads()
    @printf(stderr, "Threads: %d\n", nthreads)
    start = time()
    outputs, boundstats, allstats = SimulateParticles(
        args["geom"],
        args["flow"],
        nParticles,
        generateParticle,
        true,
        args["omega"],
        args["saveall"],
        !isnothing(args["stats"]),
        !isnothing(args["exitstats"]))
    runtime = time() - start
    if !isnothing(args["stats"])
        CSV.write(args["stats"], convert(DataFrame, allstats))
    end
    if !isnothing(args["exitstats"])
        CSV.write(args["exitstats"], convert(DataFrame, boundstats))
    end

    # Compute and display timing statistics
    @printf(stderr, "Time: %.3e\n",runtime)
    @printf(stderr, "Time per particle: %.3e\n", runtime/nParticles)
    @printf(stderr, "Time per collision: %.3e\n", runtime/sum(outputs[:,10]))
    @printf(stderr, "Interpolates: %.3e\n",interps)
    @printf(stderr, "Collides: %.3e\n",sum(outputs[:,10]))

    return allstats
end

allstats = main()
