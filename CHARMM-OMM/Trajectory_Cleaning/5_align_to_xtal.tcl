# Aligns trajectories to crystal structre
# Run using
# $> vmd -dispdev none -e 5_align_to_xtal.tcl
# see http://www.ks.uiuc.edu/Research/vmd/current/ug/node246.html for how to run scripts

# Globals
set outdir /Volumes/REA_data/AADH/traj_4_aligned
set inpdir /Volumes/REA_Data/AADH/traj_3_combined
set psf /Volumes/REA_Data/AADH/traj_2_no_water/2agy_final_nowat.psf
set pdb /Volumes/REA_Data/AADH/traj_4_aligned/2agy_c36_state0_ws_realigned.nowat.pdb



# Procedures

proc writedcd { molid outdir outname } {

    # Writes and then deletes frames
    puts [concat "Writing frames to " $outdir/$outname " from " $molid]
    set nframes [molinfo $molid get numframes]
    incr nframes -1
    animate write dcd $outdir/$outname beg 0 end $nframes skip 1 $molid

}

proc deleteframes { molid } {

    puts "Deleting frames"

    # Deletes all frames
    set nframes [molinfo $molid get numframes]
    animate delete beg 0 end $nframes skip 1 $molid

}

proc align { refid tarid sel } {

 # refid : reference molecule ID
 # tarid : target molecule ID
 # sel : the selection to align against

    puts "Aligning trajectory"
    set ref [atomselect $refid $sel frame 0]
    set tar [atomselect $tarid $sel]
    set all [atomselect $tarid all]
    set nframes [molinfo $tarid get numframes]
    puts [concat "Number of frames " $nframes ]

    for { set i 0 } { $i < $nframes } { incr i } {

        $tar frame $i
        $all frame $i
        # puts [concat "Pre-align RMSD" [measure rmsd $tar $ref]]
        # puts [concat "Aligning frame " $i ]
        $all move [measure fit $tar $ref]
        # puts [concat "Post-align RMSD" [measure rmsd $tar $ref]]
    }

    return
 }

# Load reference
set ref [mol new $pdb waitfor all]
set refid [molinfo $ref get id]
puts [concat "Added reference: " $pdb]

# Load target
set tar [mol new $psf  waitfor all]
set tarid [molinfo $tar get id]

set files [glob -directory $inpdir *.dcd]

foreach filename $files {

    set tar [mol new $psf  waitfor all]
    set tarid [molinfo $tar get id]
    mol addfile $filename type {dcd} first 0 last -1 step 1 waitfor all $tarid

    puts [concat "Added trajectory: " $filename]

    # Align trajectory
    set sel "noh"
    align $refid $tarid $sel

    # Write out trajectory
    set basename [file tail $filename]
    set outfile [string map {.dcd -aligned.dcd} $basename]
    writedcd $tarid $outdir $outfile
    mol delete $tarid
}

quit


