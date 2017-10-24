# Script to make short movies from
# Globals

set psf 2agy_rxt.psf
set traj tica-dimension-1.dcd
set outdir [file rootname $traj]
file mkdir $outdir

set mol [mol new $psf  waitfor all]
set molid [molinfo $tar get id]
mol addfile $traj type {dcd} first 0 last -1 step 1 waitfor all $molid

