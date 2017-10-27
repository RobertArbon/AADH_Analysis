# Writes pdb/psf for arbitrary selection

package require psfgen

# These are the crystal structures
set inppsf  /Users/robert_arbon/Code/AADH/MD/2b_solvate_box/2agy_c36_state0_ws_realigned.psf
set inppdb /Users/robert_arbon/Code/AADH/MD/2b_solvate_box/2agy_c36_state0_ws_realigned.pdb

# Selection TO REMOVE
set sel "not protein"

set outputdir /Users/robert_arbon/Code/AADH/Analysis/RMSD_Analysis
set outputbasename  2agy_protein

readpsf $inppsf
coordpdb $inppdb
mol load psf $inppsf pdb $inppdb

set selection [atomselect top $sel]
foreach segid [$selection get segid] resid [$selection get resid] {
    delatom $segid $resid
}

writepsf $outputdir/$outputbasename.psf
writepdb $outputdir/$outputbasename.pdb

quit

