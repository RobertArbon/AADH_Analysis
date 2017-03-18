import mdtraj as md
from glob import glob
import numpy as np
import matplotlib.pyplot as plt


traj_files = glob('/Volumes/REA_Data/AADH/traj_4_aligned/*.dcd')
topology_file = '/Users/robert_arbon/Code/AADH/Analysis/RMSD_Analysis/2agy_protein.pdb'
reference_file = '/Users/robert_arbon/Code/AADH/Analysis/RMSD_Analysis/2agy_protein.pdb'
psf_file = '/Users/robert_arbon/Code/AADH/Analysis/RMSD_Analysis/2agy_protein.psf'

selections = ['protein and backbone', 'not type H*',
              'segid BT1 and (resid 39 or resid 58)', 'segid BT2 and (resid 39 or resid 58)']


def get_rmsd(args):
    traj_file, ref_file, topology, selection = args
    print('Traj_file: {0}\nSelection: {1}'.format(traj_file, selection))
    reference = md.iterload(ref_file, top=topology, chunk=100, stride=10)
    target = md.load(traj_files, top=topology)
    print('Calculating RMSD')
    rmsd = md.rmsd(target=target, reference=reference, frame=0, atom_indices=selection)
    return rmsd


if __name__ == "__main__":

    # Get indices of selections.  We're using MDAnalysis as that has better topology recognition
    u = MDAnalysis.Universe(psf_file, reference_file)
    selection_indices = [u.select_atoms(sel).indices for sel in selections]
    traj_files = traj_files[:4]
    num_traj = len(traj_files)

    for idx, selection in enumerate(selection_indices[2:]):
        args = zip(traj_files,
                   np.repeat(reference_file, num_traj),
                   np.repeat(topology_file, num_traj),
                   np.repeat([selection], num_traj, axis=0))

        # with contextlib.closing(Pool(2)) as pool:
        #     rmsd_list = pool.map(get_rmsd, args)

        # Plot chart
        print args
        rmsd_list = get_rmsd(args[0])

        fig, axes = plt.subplots(nrows=num_traj, ncols=1, sharex=True, sharey=True)
        axes = charts.chart_trajectory(rmsd_list, row_max=num_traj, axes=axes)
        plt.savefig('{}.png'.format(selections[idx]))





