import mdtraj as md
import numpy as np
from glob import glob
from os.path import join, basename

dihed_idx = np.load('../Common/high_variance_dihedrals_active_site.npy')
traj_paths = glob('../Data/proc_traj/*as1.nc')
top_path = '../Data/proc_traj/2agy-as1.prmtop'
dest_dir = '../Data/proc_traj/high_var_dihedrals'

for traj_path in traj_paths:

    # load trajectory
    traj = md.load(traj_path, top=top_path, stride=10)
    dt = traj.timestep/1000

    # Create new file name
    fname = basename(traj_path)
    dest_path = join(dest_dir, fname.replace('.nc', '-dihedrals-dt-{}ns.npy'.format(dt)))


    # Calculate dihedrals
    diheds = md.compute_dihedrals(traj, dihed_idx)

    # Calculate sinccos
    diheds = np.concatenate([np.sin(diheds), np.cos(diheds)], axis=1)

    # Save
    print('Saving: {}'.format(dest_path))
    np.save(dest_path, diheds)