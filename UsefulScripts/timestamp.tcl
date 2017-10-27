proc enabletrace {} {
  global vmd_frame;
  trace variable vmd_frame([molinfo top]) w drawcounter
}
proc disabletrace {} {
  global vmd_frame;
  trace vdelete vmd_frame([molinfo top]) w drawcounter
}
proc drawcounter { psperframe } {

    draw delete all
    set num [molinfo top get numframes]

    # loop through the frames
    for {set i 0} {$i < $num} {incr i} {

        # go to the given frame
        animate goto $i

        ## Label text color
        draw color white

        ## Sets the 'time' (string with the time and units)
        set time [format "%8.3f ns" [expr ($i + 1) * $psperframe / 1000]]

        ## Finally, draw the text
        draw text {20 30 30} "$time" size 3

    }
}


