import MDAnalysis
from MDAnalysis.analysis.align import *
from os.path import basename, join, isfile
from glob import glob
import sys

inp_dir = 'unaligned/'
out_dir = 'aligned'
top = 'act_site.pdb'

traj_list = glob(inp_dir+'/*-as[12].dcd')
# traj_list = [inp_dir+'2agy_310k-1atm-prod{0}.{1}-as{2}.dcd'.format(i,j,k) for i in range(1,11) for j in range(1,11)
#              for k in range(1,3)]



for traj in traj_list:

    trj = MDAnalysis.Universe(top, traj)
    ref = MDAnalysis.Universe(top, traj)

    bname = basename(traj).split('.dcd')[0]
    print(bname)
    aligned_outfile = join(out_dir, '{}-aligned.dcd'.format(bname))

    if not isfile(aligned_outfile):
        # Set reference frame and align trajectory
        print('Writing {}'.format(aligned_outfile))
        ref.trajectory[0]
        rms_fit_trj(trj, ref, filename=aligned_outfile)




