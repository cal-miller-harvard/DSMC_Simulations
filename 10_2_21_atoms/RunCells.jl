using ArgParse, Luxor, Printf

include("DrawCell.jl")

function runsim(lgap, lstage, T1, T2, pflip, ismesh, flow)
    # DSMC Parameters
    zmax = 0.35 # m
    rmax = 0.03 # m
    timestep = 2E-7 # s
    nperstep = 1
    fnum = (3.0/nperstep)*3E10*flow # physical particles per simulated particle

    # Cell parameters
    gap = lgap*1000mm
    stage_len = lstage*1000mm
    stage_id = 12.7mm
    stage_od = 22.23mm
    raperture = 4.4mm
    facethickness = 0.6mm

    # Particle simulation parameters
    n_particles = 1000000
    omegas = [0.0]
    zmin = 65.09E-3
    zend = zmin+lstage+1E-3
    zmaxs = [zend]

    σs = [pi * (327E-12).^2] #pi .* (140E-12 .+ 1E-12 .* [227, 231, 242]).^2 # collision cross section (m^2)
    Ms = [23.0] #[23.0, 40.0, 173.0] # mass of molecule (AMU)

    # Paths
    PROG_PATH = pwd()
    TEMPLATE_PATH = "./template"
    RUN_PATH = @sprintf("./flow_%.5f_T2_%.5f", flow, T2)
    SPARTA_CMD = `mpirun /n/home03/calmiller/programs/sparta/spa -kokkos off`

    mkpath(RUN_PATH)
    for f in readdir(TEMPLATE_PATH)
        cp(TEMPLATE_PATH*"/"*f, RUN_PATH*"/"*f, force=true)
    end
    mkpath(RUN_PATH*"/data")

    # Define cell geometry
    Drawing(round(Int32,zmax*1000mm), round(Int32,rmax*1000mm), RUN_PATH*"/cell.pdf")

    stage = Point(gap, 0) + make_stage(stage_len, stage_id, stage_od, raperture, facethickness)
    mesh = Point(BoundingBox(stage).corner2.x, 0.14mm) .+ make_mesh(:y, 4, 0.675mm, 0.675mm, 1.055mm)

    if ismesh
        polygons = [stage, mesh...]
    else
        polygons = [stage]
    end
    
    # Draw and save cell image
    sethue("black")
    setline(0.1mm)
    for p in polygons
        poly(p, :stroke, close=true)
    end
    finish()

    try
        cmd = "firefox"
        arg = RUN_PATH*"/cell.pdf"
        run(`$cmd $arg`)
    catch
        println("firefox not found")
    end

    toSPARTA(polygons, RUN_PATH*"/data.stages")

    # if T1 > 1.9 && T2 > 1.9
    #     he = "he4"
    # else
    #     he = "he3"
    # end

    he = "he3"

    open(RUN_PATH*"/in.cell", "w") do f
        write(f,@sprintf("""# 2d axial simulation of CBGB cell

        echo screen

        # variables to be set in RunCells.jl
        variable FNUM equal "%.3E"
        variable ZMAX equal "%.3E"
        variable RMAX equal "%.3E"
        variable TIMESTEP equal "%.3E"
        variable NPERSTEP equal "%.3E"
        variable T1 equal "%.3E"
        variable T2 equal "%.3E"
        
        global		        fnum \${FNUM} cellmax 10000
        seed	    	    12345
        dimension   	    2
        units               si
        timestep 	        \${TIMESTEP}
        
        species		        %s.species He
        mixture		        He He nrho 1 vstream 26.2 0 0 temp \${T1}
        collide             vss He %s.vss #kk
        
        # read_restart data/cell.restart.*
        
        boundary	        o ao p
        create_box  	    0 \${ZMAX} 0 \${RMAX} -0.1 0.1
        create_grid 	    20 20 1 
        balance_grid        rcb cell
        
        read_surf           data.cell
        group       inlet surf id 12
        fix		    in emit/surf He inlet n \${NPERSTEP} perspecies no

        read_surf           data.diffuser
        
        read_surf           data.stages trans 65.09E-3 0 0 group stages
        
        surf_collide	    diffuse diffuse \${T1} 1.0
        surf_collide	    cold diffuse \${T2} 1.0
        surf_modify         all collide diffuse
        surf_modify         stages collide cold
        
        stats		    1000
        stats_style	    wall step cpu np nattempt ncoll nscoll nscheck
        
        # Run until convergence
        label loop
        variable a loop 125
        run 		    20000
        adapt_grid all refine particle 16 4
        balance_grid rcb part
        next a
        jump in.cell loop
        
        compute temp thermal/grid all He temp
        compute rhov grid all He nrho massrho u v w
        fix out ave/grid all 1000 100 100000 c_temp[*] c_rhov[*]

        run 1
        # Save statistics and surfaces
        dump out grid all 100000 data/DS2FF.DAT xc yc f_out[*]
        dump surfs surf all 100000 data/cell.surfs id v1x v1y v2x v2y
        
        # Record statistics
        label loop2
        variable b loop 10
        run 		    10000
        
        next b
        jump in.cell loop2
        write_restart data/restart.slurm""", fnum, zmax, rmax, timestep, nperstep, T1, T2, he, he))
    end

    cd(RUN_PATH)
    # run(pipeline(SPARTA_CMD, stdin="in.cell"), wait=true)

    if he == "he3"
        m = 3.0
    else
        m = 4.0
    end

    for (i, omega) in enumerate(omegas)
        for (j, M) in enumerate(Ms)
            fname = @sprintf("run_omega_%.5f_M_%.1f.cell",omega, M)
            open(fname, "w") do f
                write(f,@sprintf("""#!/bin/bash
                #SBATCH -n 16 # Number of cores requested
                #SBATCH -N 1 # Ensure that all cores are on one machine
                #SBATCH -t 0-08:00 # Runtime in minutes
                #SBATCH -p shared # Partition to submit to
                #SBATCH --mem-per-cpu 1024 # Memory per cpu in MB
                #SBATCH --open-mode=append
                #SBATCH -o data/particles_omega_%.5f_M_%.1f_zmax_%.5f_pflip_%.5f_job_%%j.out # Standard out goes to this file
                #SBATCH -e data/particles_omega_%.5f_M_%.1f_zmax_%.5f_pflip_%.5f_job_%%j.err # Standard err goes to this filehostname

                module load intel/19.0.5-fasrc01 openmpi/4.0.2-fasrc01 fftw/3.3.8-fasrc01 cmake/3.12.1-fasrc01 Anaconda3/2019.10 python/3.7.7-fasrc01
                module list

                export OMP_PROC_BIND=spread
                export OMP_PLACES=threads
                export JULIA_NUM_THREADS=\$SLURM_CPUS_ON_NODE

                cd data
                pwd
                echo "running...."

                julia /n/home03/calmiller/DSMC_Simulations/ParticleTracing/ParticleTracing.jl -z 0.035 -T 2.0 -n %d ./cell.surfs ./DS2FF.DAT --omega %.5f --pflip %.5f -m %.5f -M %.5f --sigma %.5E --zmin %.5f --zmax %.5f --saveall 0""", omega, M, zmaxs[i], pflip, omega, M, zmaxs[i], pflip, n_particles, omega, pflip, m, M, σs[j], zmin, zmaxs[i]))
            end
            run(`sbatch $fname`)
        end
    end
    cd(PROG_PATH)
end

"""
    parse_commandline()

Parses command-line arguments.
"""
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "-l"
            help = "gap between first and second stage (m)"
            arg_type = Float64
            default = 0.001
        "-L"
            help = "length of second stage (m)"
            arg_type = Float64
            default = 0.040
        "--T1"
            help = "temperature of first stage (K)"
            arg_type = Float64
            default = 2.0
        "--T2"
            help = "temperature of second stage (K)"
            arg_type = Float64
            default = 2.0
        "--pflip"
            help = "collision spin flip probability"
            arg_type = Float64
            default = 0.1
        "--flow"
            help = "flow rate (sccm)"
            arg_type = Float64
            default = 4.0
        "--mesh"
            help = "is there a mesh?"
            arg_type = Bool
            default = true
    end

    return parse_args(s)
end

args = parse_commandline()
runsim(args["l"], args["L"], args["T1"], args["T2"], args["pflip"], args["mesh"], args["flow"])
