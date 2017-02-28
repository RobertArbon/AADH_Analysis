import mdtraj as md
from glob import glob
import MDAnalysis

traj_files = glob('/Volumes/REA_Data/AADH/traj_4_aligned/*.dcd')
topology_file = '/Users/robert_arbon/Code/AADH/Analysis/RMSD_Analysis/2agy_protein.pdb'
reference_file = '/Users/robert_arbon/Code/AADH/Analysis/RMSD_Analysis/2agy_protein.pdb'
psf_file = '/Users/robert_arbon/Code/AADH/Analysis/RMSD_Analysis/2agy_protein.psf'

selections = ['protein and backbone', 'not type H*',
              'segid BT1 and (resid 39 or resid 58)', 'segid BT2 and (resid 39 or resid 58)',
              'around 15.0 segid BT1 and resid 39','around 15.0 segid BT2 and resid 39']


# Get indices of selections.  We're using MDAnalysis as that has better topology recognition
u = MDAnalysis.Universe(psf_file, reference_file)
selection_indices = [u.select_atoms(sel).indices for sel in selections]
for sel in selection_indices:
    print(len(sel))
# Load trajectories and reference

target = md.load(traj_files[0], top=topology_file)
reference = md.load(reference_file)

