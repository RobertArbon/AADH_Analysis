proc make_rotation_animated_gif { movie_name } {
	set frame 0
	set rot 60
	for {set i 0} {$i < 360} {incr i $rot} {
		set filename snap.[format "%04d" $frame].rgb
		render snapshot $filename
		incr frame
		rotate y by $rot
	}
	exec convert -loop 1 -delay 10 snap*rgb $movie_name.gif
	exec convert snap*rgb $movie_name.png
	file delete {*}[glob snap*rgb]
}

# set sel [atomselect top "protein"]

set outfile cluster0

display resetview
display update on
display update ui
display distance -4

mol modstyle 0 top Licorice
mol modcolor 0 top Name
mol drawframes top 0 1:100

make_rotation_animated_gif $outfile
