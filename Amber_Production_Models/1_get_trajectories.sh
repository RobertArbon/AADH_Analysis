#!/bin/bash
# Creates a directory of symlinks for trajectories and a prmtop file

# Where the data is stored
rootdir=/Volumes/dynamics_of_glassy_aerosols/AADH/Amber/round_1
subdirs=({1..100}ns)

oldtrajname=100ns-production.mdcrd
oldtopname=2agy_final_min.prmtop

newtrajname=2agy
newtopname=2agy.prmtop

trajsize=13822480796
# Make the directory to store the trajectories
newdir=raw_traj
mkdir -p $newdir

# link trajectories
for dir in ${subdirs[*]}
do
    trajfile=$rootdir/$dir/$oldtrajname
    if [ -e $trajfile ]
    then
        actualsize=$(wc -c  < $trajfile)

        if [ $actualsize -eq $trajsize ]
        then
            ln -s $trajfile $newdir/$newtrajname-${dir%ns}.nc
        fi
    fi
done

# link prmtop
ln -s $rootdir/$subdirs[0]/$newtopname $newdir/$newtopname

