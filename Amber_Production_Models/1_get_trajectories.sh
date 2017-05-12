#!/bin/bash
# Creates a directory of symlinks for trajectories and a prmtop file

# Where the data is stored
#rootdir=/Volumes/dynamics_of_glassy_aerosols/AADH/Amber/round_1
rootdir=~/RDSF/Amber/round_1
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
            newfile=$newdir/$newtrajname-${dir%ns}.nc
            if [[ ! -h $newfile ]]
            then
              echo 'Creating link to ' $newfile
              ln -s $trajfile $newfile
            fi
        fi
    fi
done

# link prmtop
#ln -s $rootdir/$subdirs[0]/$oldtopname $newdir/$newtopname
ln -s /home/robert/Research/AADH/MD/common/$oldtopname $newdir/$newtopname
