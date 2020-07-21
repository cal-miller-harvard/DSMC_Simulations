using ArgParse, Luxor, Printf

include("DrawCell.jl")

const SCCM = 4.41941E17 # particles per second

function runsim(lgap, lstage, T1, T2, pflip)
    # Methanol
    methanol_flow = 0.08 # sccm
    methanol_T = 250 # K 
    pstick = 0.5 # methanol sticking probability

    # DSMC Parameters
    flow = 4.0 # sccm
    zmax = 0.35 # m
    rmax = 0.03 # m
    timestep = 2E-7 # s
    methanolnperstep = 1
    nperstep = 1
    methanolnevery = round(Int64, flow/methanol_flow * nperstep/methanolnperstep)
    fnum = flow * SCCM * timestep / nperstep


    # Cell parameters
    gap = lgap*1000mm
    stage_len = lstage*1000mm
    stage_id = 12.7mm
    stage_od = 22.23mm
    raperture = 4.5mm
    facethickness = 0.6mm

    # Paths
    PROG_PATH = pwd()
    TEMPLATE_PATH = "./template"
    RUN_PATH = @sprintf("./flow_%.5f_gap_%.5f_len_%.5f_T1_%.5f_T2_%.5f_methanol_%.5f_T_%.5f", flow, gap/(1000mm), lstage, T1, T2, methanol_flow, methanol_T)
    SPARTA_CMD = `mpirun /n/home03/calmiller/programs/sparta/spa -kokkos off`
    # SPARTA_CMD = `/usr/local/bin/spa_kokkos_omp`

    mkpath(RUN_PATH)
    for f in readdir(TEMPLATE_PATH)
        cp(TEMPLATE_PATH*"/"*f, RUN_PATH*"/"*f, force=true)
    end
    mkpath(RUN_PATH*"/data")

    # Define cell geometry
    Drawing(round(Int32,zmax*1000mm), round(Int32,rmax*1000mm), RUN_PATH*"/cell.pdf")

    # First Stage
    max_len = 5mm
    first = divide_polygon(1mm*[
        Point(0.79, 0),
        Point(0.79, 1.59),
        Point(3.18, 1.59),
        Point(3.18, 12.7),
        Point(10.27, 12.7),
        Point(10.27, 21.59),
        Point(56.52, 21.59),
        Point(56.52, 3.5),
        Point(58.1, 3.5),
        Point(58.1, 25.4),
        Point(0, 25.4),
        Point(0, 0)
    ], max_len)

    # Diffuser
    diffuser = divide_polygon(1mm*[
        Point(7.14, 0),
        Point(7.14, 1.59),
        Point(9.53, 1.59),
        Point(9.53, 9.53),
        Point(6.35, 9.53),
        Point(6.35, 0)
    ], max_len)

    # Second Stage
    stage = Point(BoundingBox(first).corner2.x + gap, 0) + make_stage(stage_len, stage_id, stage_od, raperture, facethickness)
    mesh = Point(BoundingBox(stage).corner2.x, 0.1mm) .+ make_mesh(:y, 7, 0.425mm, 0.425mm, 0.624mm)
    lfirst = BoundingBox(first).corner2.x/(1000*mm)

    polygons = [first, diffuser, stage, mesh...]
    He_inlet = 1
    methanol_inlet = length(first) + 1

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

    if T1 > 1.9 && T2 > 1.9
        he = "he4"
    else
        he = "he3"
    end

    open(RUN_PATH*"/methanol.surf", "w") do f
        write(f, @sprintf("""methanol --> NULL
        R S %.3E""", pstick))
    end

    open(RUN_PATH*"/in.cell", "w") do f
        write(f,@sprintf("""# 2d axial simulation of CBGB cell

        echo screen

        # variables to be set in RunCells.jl
        variable FNUM equal "%.3E"
        variable ZMAX equal "%.3E"
        variable RMAX equal "%.3E"
        variable TIMESTEP equal "%.3E"
        variable NPERSTEP equal "%.3E"
        variable METHANOLPERSTEP equal "%.3E"
        variable METHANOLNEVERY equal "%d"
        variable T1 equal "%.3E"
        variable T2 equal "%.3E"
        variable TM equal "%.3E"
        variable FIRSTSTAGE equal "%.3E"
        
        global		        fnum \${FNUM} cellmax 10000
        seed	    	    12345
        dimension   	    2
        units               si
        timestep 	        \${TIMESTEP}
        
        species		        %s.species He
        mixture		        He He nrho 1 vstream 26.2 0 0 temp \${T1}
        species		        methanol.species methanol
        mixture		        methanol methanol nrho 1 vstream 26.2 0 0 temp \${TM}

        collide             vss all %s.vss #kk

        
        # read_restart data/cell.restart.*
        
        boundary	        o ao p
        create_box  	    0 \${ZMAX} 0 \${RMAX} -0.1 0.1
        create_grid 	    20 20 1 
        balance_grid        rcb cell
        
        read_surf   data.stages
        group       inlet surf id 1
        group       methanolinlet surf id %d
        region      stagesr cylinder z 0 0 100 \${FIRSTSTAGE} 100
        group       stages surf region stagesr one

        group       hot surf id %d:%d

        fix		    in emit/surf He inlet n \${NPERSTEP} perspecies no
        fix         in2 emit/surf methanol methanolinlet n \${METHANOLPERSTEP} perspecies no nevery \${METHANOLPERSTEP}
        
        surf_collide	    diffuse diffuse \${T1} 1.0
        surf_collide	    cold diffuse \${T2} 1.0
        surf_collide        hot diffuse \${TM} 1.0

        surf_react          stick prob methanol.surf

        surf_modify         all collide diffuse react stick
        surf_modify         stages collide cold
        surf_modify         hot collide hot react none
        
        stats		    1000
        stats_style	    wall step cpu np nattempt ncoll nscoll nscheck
        
        # Run until convergence
        label loop
        variable a loop 50
        run 		    20000
        adapt_grid all refine particle 16 4
        balance_grid rcb part
        next a
        jump in.cell loop
        
        compute temp thermal/grid all He temp
        compute rhov grid all He nrho massrho u v w
        fix out ave/grid all 1000 100 100000 c_temp[*] c_rhov[*]

        compute methanoltemp thermal/grid all methanol temp
        compute methanolrhov grid all methanol nrho massrho u v w
        fix methanolout ave/grid all 1000 100 100000 c_methanoltemp[*] c_methanolrhov[*]
        
        compute collisions react/surf all stick r:methanol p:methanol
        fix collisions ave/surf all 1000 100 100000 c_collisions[*]

        run 100000
        # Save statistics and surfaces
        dump out grid all 100000 data/DS2FF.DAT xc yc f_out[*]
        dump out2 grid all 100000 data/DS2FF.methanol.DAT xc yc f_methanolout[*]
        dump surfs surf all 100000 data/cell.surfs id v1x v1y v2x v2y f_collisions[*]

        dump_modify out append no format float %.5f
        dump_modify out2 append no format float %.5f
        dump_modify surfs append no format float %.5f
        
        # Record statistics
        label loop2
        variable b loop 10
        run 		    10000
        
        next b
        jump in.cell loop2
        write_restart data/restart.slurm""", fnum, zmax, rmax, timestep, nperstep, methanolnperstep, methanolnevery, T1, T2, methanol_T, lfirst + 0.0005, he, he, methanol_inlet, methanol_inlet, methanol_inlet+2))
    end

    cd(RUN_PATH)
    run(pipeline(SPARTA_CMD, stdin="in.cell"), wait=true)
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
            default = 0.002
        "-L"
            help = "length of second stage (m)"
            arg_type = Float64
            default = 0.020
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
    end

    return parse_args(s)
end

args = parse_commandline()
runsim(args["l"], args["L"], args["T1"], args["T2"], args["pflip"])
