# Grabs aligned crystal structure (aligned to x axis for convenience) and removes waters
# Run using
# $> vmd -dispdev none -e 4_import_xtal.tcl
# see http://www.ks.uiuc.edu/Research/vmd/current/ug/node246.html for how to run scripts

set outdir /Volumes/REA_data/AADH/traj_4_aligned

# Get reference structure
### file names
set filename /Users/robert_arbon/Code/AADH/MD/2b_solvate_box/2agy_c36_state0_ws_realigned.pdb
set basename [file tail $filename]
set outname [string map {.pdb .nowat.pdb} $basename ]

### Load data and select protein
set mol [mol new $filename type pdb waitfor all]
set sel [atomselect $mol "protein"]

### Write out data
$sel writepdb $outdir/$outname

quit