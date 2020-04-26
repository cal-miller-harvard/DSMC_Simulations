import os


# Parameters:
# 4 sccm
# 5, 15, 25 mm gap
# 100 mm second stage length

flows = [4] # sccm
gaps = [5E-3, 15E-3, 25E-3] # m
lengths =[100E-3] # m
for flow in flows:
    for gap in gaps:
        for length in lengths:
            directory = "flow_{:.3f}_gap_{:.3f}_len_{:.3f}".format(flow, gap, length)
            os.system("cp -rf template "+directory)
            os.chdir(directory)
            with open(r"in.cell", "w") as f:
                f.write("""# 2d axial simulation of CBGB cell

# Gap between 1st and second stage (mm)
variable GAP equal {:.3f} # 2.7E-3
# Length of second stage (mm)
variable LEN equal {:.3f} # 11.1E-3
# Input density: 1 sccm ~ 2E21?
variable RHO equal "{:.3f} * 2E21"

global		        fnum 1.0E12 cellmax 10000 gridcut 0.0 comm/sort yes
seed	    	    12345
dimension   	    2
units               si
timestep 	        1E-6

species		        he4.species He4
mixture		        He4 He4 nrho $(v_RHO) vstream 26.2 0 0 temp 2.0
collide             vss He4 He4.vss #kk

# read_restart data/cell.restart.*

boundary	        o ao p
create_box  	    0 0.20 0 0.03 -0.1 0.1
create_grid 	    20 20 1 
balance_grid        rcb cell

read_surf           data.cell
group       inlet surf id 12
fix		    in emit/surf He4 inlet

read_surf           data.diffuser group diffuser
read_surf           data.secondstage atrans $(v_GAP - 2.7E-3) 0 0 group secondstage

read_surf           data.mesh atrans $(v_GAP - 2.7E-3) 0 0 group mesh
read_surf           data.mesh atrans $(v_GAP - 2.7E-3) 1.21E-3 0 group mesh
read_surf           data.mesh atrans $(v_GAP - 2.7E-3) 2.42E-3 0 group mesh
read_surf           data.mesh atrans $(v_GAP - 2.7E-3) 3.63E-3 0 group mesh

surf_collide	    diffuse diffuse 2.0 1.0
surf_modify         all collide diffuse

region mesh block 78E-3 82E-3 0 5E-3 -1 1
group mesh grid region mesh one
adapt_grid mesh refine random 1 0 iterate 2


region ssbody cylinder x 0 0 1 $(v_GAP + 65.09E-3 + 1E-3) 1
group ssbody surf region ssbody all
move_surf ssbody trans $(v_LEN-11.11E-3) 0 0 connect yes

stats		    1000
stats_style	    wall step cpu np nattempt ncoll nscoll nscheck

# Run until convergence
label loop
variable a loop 40
run 		    10000
adapt_grid all refine particle 16 4
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
run 1
""".format(gap, length, flow))
            os.system("sbatch run.slurm")
            os.chdir("..")
