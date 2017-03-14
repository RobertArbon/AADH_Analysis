proc make_rotation_animated_gif { movie_name } {
	set frame 0
	set rot 5
	for {set i 0} {$i < 360} {incr i $rot} {
		set filename snap.[format "%04d" $frame].rgb
		render snapshot $filename
		incr frame
		rotate y by $rot
	}
	exec convert -loop 1 -delay 10 snap*rgb $movie_name
	file delete {*}[glob snap*rgb]
}

set files [glob -directory /Users/robert_arbon/Code/AADH/Analysis/MSM_Reactants_Only/Clustering/clusters/ *.dcd]
set psf /Users/robert_arbon/Code/AADH/Analysis/MSM_Reactants_Only/2agy_rxt.psf
set dcd /Users/robert_arbon/Code/AADH/Analysis/MSM_Reactants_Only/Clustering/clusters/rmsd_cluster-0.dcd

foreach file $files {
    set mol [mol new $psf waitfor all]
    mol addfile $file type {dcd} first 0 last -1 step 1 waitfor all $mol

    display resetview
    display update on
    display update ui

    mol modstyle 0 $mol Licorice
    mol modcolor 0 $mol Name
    mol drawframes $mol 0 1:100

    set basename [file tail $file]
    set outfile [string map ".dcd  .gif"  $basename]
    make_rotation_animated_gif $outfile

    mol delete $mol

}

quit