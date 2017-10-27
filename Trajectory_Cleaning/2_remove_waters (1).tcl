# Loads trajectory file and removes water from DCD
# Run using
# $> vmd -dispdev none -e remove_waters.tcl
# see http://www.ks.uiuc.edu/Research/vmd/current/ug/node246.html for more info. 

set psf /Users/robert_arbon/Code/AADH/MD/common/2agy_final.psf
set inpdir /Volumes/REA_data/AADH/traj_1_raw
set root 2agy-310k-1atm-prod
set outdir /Volumes/REA_data/AADH/traj_2_no_water

# Load single frame for writing reduced psf
set mol [mol new $psf type psf waitfor all]
set filename ${root}1.1-1.dcd
mol addfile $inpdir/$filename type {dcd} first 0 last 0 step 1 waitfor all top
set prot [atomselect $mol "protein"]

# Write the psf
$prot writepsf  $outdir/2agy_final_nowat.psf
mol delete top

# Load full DCDs and remove waters
set files [glob -directory $inpdir *.dcd]
set mol [mol new $psf type psf waitfor all]

foreach filename $files {

    mol addfile $filename type {dcd} first 0 last -1 step 1 waitfor all top
    set prot [atomselect $mol "protein"]

    set basename [file tail $filename]

    # Write out protein only DCD
    set num [molinfo top get numframes] 
    incr num -1 
    set out [string map {.dcd .nowat.dcd} $basename ] 
    animate write dcd $outdir/$out beg 0 end $num skip 1 waitfor all sel $prot top 
    animate delete all  
}

quit