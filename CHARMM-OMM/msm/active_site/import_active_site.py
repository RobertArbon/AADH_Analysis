import MDAnalysis
from os.path import join, isfile
from glob import glob
from multiprocessing import Pool, cpu_count

psf = '../../common/2agy_final.psf'
as_num = 2
selection = 'segid BT{} and (resid 39 or resid 58)'.format(as_num)
label = '-as{}'.format(as_num)
dirname = 'unaligned'

def write_backbone(traj):

    bname = traj.split('/')[-1].split('.dcd')[0]
    outname = join(dirname, bname+label+'.dcd')

    if not isfile(outname):
        print 'Writing ', outname
        u = MDAnalysis.Universe(psf, traj)
        atom_sel = u.select_atoms(selection)
        natoms = atom_sel.n_atoms
        with MDAnalysis.Writer(outname, multiframe=True, n_atoms=atom_sel.n_atoms) as W:
            for ts in u.trajectory:
                W.write(atom_sel)

if __name__ == '__main__':
    traj_list = glob('/Volumes/REA_data/AADH/MD/Trajectories/*.dcd')
    p = Pool(cpu_count())
    p.map(write_backbone, traj_list)