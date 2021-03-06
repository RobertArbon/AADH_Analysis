#!/usr/bin/env bash
rootdir=. #/Users/robert_arbon/Code/AADH/Analysis/Amber_Production_Models

topfile=2agy.prmtop
mkdir -p proc_traj
inpdir=raw_traj
outdir=proc_traj

for traj in raw_traj/*
do
    bn=$(basename $traj)
    if [[ ${bn#*.} == "nc" ]]
    then
        fname=${bn%.*}

        newfile=$rootdir/$outdir/$fname-as1.nc
        if [[ ! -e $newfile ]]
        then
            echo 'Processing ' $traj 'saving as ' $newfile
            cpptraj -p $rootdir/$inpdir/$topfile -i site1.cpptraj -y $rootdir/$inpdir/$bn -x $newfile \
            --log $outdir/as1.log
        fi

        newfile=$rootdir/$outdir/$fname-as2.nc
        if [[ ! -e $newfile ]]
        then
            echo 'Processing ' $traj 'saving as ' $newfile
            cpptraj -p $rootdir/$inpdir/$topfile -i site2.cpptraj -y $rootdir/$inpdir/$bn -x $newfile \
            --log $outdir/as2.log
        fi
    fi
done
