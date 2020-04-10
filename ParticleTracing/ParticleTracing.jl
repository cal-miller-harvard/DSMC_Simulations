using Distributed
using Random, LinearAlgebra, NearestNeighbors
using CSV, DataFrames, Printf

# Constants
const MASS_PARTICLE = 174 + 16 + 1    # AMU
const MASS_BUFFER_GAS = 4             # AMU
const MASS_REDUCED = MASS_PARTICLE * MASS_BUFFER_GAS /
    (MASS_PARTICLE + MASS_BUFFER_GAS) # AMU
const kB = 8314.46                    # AMU m^2 / (s^2 K)
const σ_BUFFER_GAS_PARTICLE = 100E-20 # m^2

"""
    collide!(v, vgx, vgy, vgz, T)

Accepts as input the velocity of a particle v, the mean velocity of a buffer gas atom vgx, vgy, vgz, and the buffer gas temperature T. Computes the velocity of the particle after they undergo a collision, treating the particles as hard spheres and assuming a random scattering parameter and buffer gas atom velocity (assuming the particles are moving slower than the buffer gas atoms). Follows Appendices B and C of Boyd 2017. Note that v and vg are modified.
"""
@inline function collide!(v::Vector, vgx::Number, vgy::Number, vgz::Number, T::Number)
    vg = sqrt(-2 * kB * T / MASS_BUFFER_GAS * log(1-Random.rand()))
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
    freePropagate(vrel, T, ρ)

Accepts as input the velocity of a particle relative to a buffer gas atom v, the buffer gas temperature T and the buffer gas density ρ. Draws a distance the particle travels before it hits a buffer gas atom from an exponential distribution, accounting for a velocity-dependent mean free path. Note that this assumes that the gas properties don't change significantly over a mean free path.
"""
@inline function freePropagate(vrel::Number, T::Number, ρ::Number)
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
    propagate(xinit, vin, interp!, getCollision)

Accepts as input the position of a particle xinit, its velocity v, the function interp!, which takes as input a position and a vector and updates the vector to describe the gas [x, y, vgx, vgy, vgz, T, ρ], and the function getCollision(x1, x2), which returns whether if the segment from x1 to x2 intersects geometry. Computes the path of the particle until it getCollision returns true. Returns a vector of simulation results with elements x, y, z, xnext, ynext, znext, vx, vy, vz, collides, time.
"""
@inline function propagate(xinit::Vector, vin::Vector, interp!::Function, getCollision::Function)
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
        dist = freePropagate(vrel, props[6], props[7])
        xnext .= x .+ dist .* LinearAlgebra.normalize(v)
        if getCollision(x, xnext)
            return (x[1], x[2], x[3], xnext[1], xnext[2], xnext[3], v[1], v[2], v[3], collides, time)
        else
            time += dist / LinearAlgebra.norm(v)
            collides += 1
        end
        x .= xnext
        collide!(v, props[3], props[4], props[5], props[6])
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
    generateParticle::Function
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

    @printf("Mean dmin: %f\n", sum(grids[:,8]/size(grids)[1]))

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
            return true
        end
        for i in 1:size(geom)[1]
            if getIntersection(geom[i,2],geom[i,3],geom[i,4],geom[i,5],x1[3],r1,x2[3],r2)
                return true
            end
        end
        return false
    end

    output_dim = length(propagate(zeros(3), zeros(3), interpolate!, (x,y)->true))
    outputs = zeros(nParticles, output_dim)
    Threads.@threads for i in 1:nParticles
    # for i in 1:nParticles
        xpart, vpart = generateParticle()
        outputs[i,:] .= propagate(xpart, vpart, interpolate!, getCollision)
    end
    
    return outputs
end

# Set simulation parameters
nParticles = 100000
# dir="/home/cal/Documents/DSMC_Simulations/4_3_20_flow_gap_length/flow_1.000_gap_0.003_len_0.010/data/"
dir="/home/cal/Documents/DSMC_Simulations/4_3_20_flow_gap_length/flow_4.000_gap_0.003_len_0.010/data/"
generateParticle() = ([0.0, 0.0, 0.035], zeros(3))

# Run simulation
start = time()
outputs = SimulateParticles(
    dir*"cell.510001.surfs",
    dir*"DS2FF.500000.DAT",
    nParticles,
    generateParticle)
runtime = time() - start

# Write output file
CSV.write(dir*"particles.out",
convert(DataFrame,outputs),
delim = ' ',
header = ["x", "y", "z", "xnext", "ynext", "znext", "vx", "vy", "vz", "collides", "time"])

# Compute timing statistics
nthreads = Threads.nthreads()
@printf("Time: %.3e\n",runtime)
@printf("Time per particle: %.3e\n", runtime/nParticles)
@printf("Time per collision: %.3e\n", runtime/sum(outputs[:,10]))
@printf("Interpolates: %.3e\n",interps)
@printf("Collides: %.3e\n",sum(outputs[:,10]))