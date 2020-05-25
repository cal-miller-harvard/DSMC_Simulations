using ArgParse

include("DrawCell.jl")

function runsim(lwallarg)
    # DSMC Parameters
    flow = 4 # sccm
    zmax = 0.35 # m
    rmax = 0.03 # m
    timestep = 2E-7 # s
    nperstep = 3
    fnum = 3E10*flow # physical particles per simulated particle

    # Cell parameters
    nstages = 10
    wallthickness = 2mm
    wallperiod = 10mm
    lwall = lwallarg*1000mm
    # gap = wallperiod - lwall
    gap = lwallarg*1000mm

    stage_len = 10mm
    stage_id = 12.7mm
    stage_od = 22.23mm
    raperture = 4.5mm
    facethickness = 0.6mm

    # Paths
    PROG_PATH = pwd()
    TEMPLATE_PATH = "./template"
    RUN_PATH = @sprintf("./flow_2.000_gap_%.5f_len_0.000", gap/(1000mm))
    SPARTA_CMD = `mpirun /n/home03/calmiller/programs/sparta/spa -kokkos off`    

    mkpath(RUN_PATH)
    for f in readdir(TEMPLATE_PATH)
       cp(TEMPLATE_PATH*"/"*f, RUN_PATH*"/"*f, force=true)
    end
    mkpath(RUN_PATH*"/data")

    # Define cell geometry
    Drawing(round(Int32,zmax*1000mm), round(Int32,rmax*1000mm), RUN_PATH*"/cell.pdf")

    # side = Point(gap, stage_id) .+ make_mesh(:x, nstages, wallthickness, lwall, wallperiod)
    # stage = Point(gap+maximum([a.corner2.x for a in BoundingBox.(side)]), 0) + make_stage(stage_len, stage_id, stage_od, raperture, facethickness)
    # mesh = Point(BoundingBox(stage).corner2.x, 0.36mm) .+ make_mesh(:y, 4, 0.48mm, 0.48mm, 1.21mm)

    stage = Point(gap, 0) + make_stage(stage_len, stage_id, stage_od, raperture, facethickness)
    mesh = Point(BoundingBox(stage).corner2.x, 0.36mm) .+ make_mesh(:y, 4, 0.48mm, 0.48mm, 1.21mm)

    # polygons = [stage, side..., mesh...]
    polygons = [stage, mesh...]

    # 
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


    open(RUN_PATH*"/in.cell", "w") do f
        write(f,@sprintf("""# 2d axial simulation of CBGB cell

        echo screen

        # variables to be set in RunCells.jl
        variable FNUM equal "%.3E"
        variable ZMAX equal "%.3E"
        variable RMAX equal "%.3E"
        variable TIMESTEP equal "%.3E"
        variable NPERSTEP equal "%.3E"
        
        global		        fnum \${FNUM} cellmax 10000
        seed	    	    12345
        dimension   	    2
        units               si
        timestep 	        \${TIMESTEP}
        
        species		        he4.species He4
        mixture		        He4 He4 nrho 1 vstream 26.2 0 0 temp 2.0
        collide             vss He4 He4.vss #kk
        
        # read_restart data/cell.restart.*
        
        boundary	        o ao p
        create_box  	    0 \${ZMAX} 0 \${RMAX} -0.1 0.1
        create_grid 	    20 20 1 
        balance_grid        rcb cell
        
        read_surf           data.cell
        group       inlet surf id 12
        fix		    in emit/surf He4 inlet n \${NPERSTEP} perspecies no

        read_surf           data.diffuser
        
        read_surf           data.stages trans 65.09E-3 0 0 group stages
        
        surf_collide	    diffuse diffuse 2.0 1.0
        surf_collide	    cold diffuse 0.5 1.0
        surf_modify         all collide diffuse
        surf_modify         stages collide cold
        
        stats		    1000
        stats_style	    wall step cpu np nattempt ncoll nscoll nscheck
        
        # Run until convergence
        label loop
        variable a loop 40
        run 		    10000
        adapt_grid all refine particle 16 4
        balance_grid rcb part
        next a
        jump in.cell loop
        
        compute temp thermal/grid all He4 temp
        compute rhov grid all He4 nrho massrho u v w
        fix out ave/grid all 1000 100 100000 c_temp[*] c_rhov[*]
        dump out grid all 100000 data/DS2FF.*.DAT xc yc f_out[*]
        
        # Record statistics
        label loop2
        variable b loop 11
        run 		    10000
        
        next b
        jump in.cell loop2
        
        # Save surfaces
        dump surfs surf all 1 data/cell.*.surfs id v1x v1y v2x v2y
        run 1""", fnum, zmax, rmax, timestep, nperstep))
    end

    println("cat $RUN_PATH/in.cell")

    cd(RUN_PATH)
    run(pipeline(SPARTA_CMD, stdin="in.cell"), wait=true)
    run(`sbatch particles_0.slurm`)
    run(`sbatch particles_300.slurm`)
    run(`sbatch particles_600.slurm`)
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
            default = 0.003
    end

    return parse_args(s)
end

args = parse_commandline()
runsim(args["l"])
