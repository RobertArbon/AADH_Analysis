#
#   LOADS DATA
#
#


import os
from msmbuilder.io import GenericParser, save_meta, gather_metadata
#
# Set up symlinks for easy referencing
#


if not os.path.exists('traj'):
    os.symlink('/Volumes/REA_Data/AADH/traj_5_rxts', 'traj')
if not os.path.exists('topology'):
    os.symlink('/Users/robert_arbon/Code/AADH/Analysis/MSM_Reactants_Only/2agy_rxt.pdb', 'topology.pdb')

#
# Helper functions
#


def identity(x):
    return x


#
# File name parsing and metadata
#


re_pattern = '(\w+)-([0-9]{3})k-([0-9])atm-prod([0-9]+\.[0-9]+).*BT([0-9]+)*'
captured_group_names = ['PDB', 'Temp', 'Pressure', 'Prod_Round', 'Act_Site']
captured_group_transforms = [identity, float, float, identity, int]
time_step = 1
file_type = 'dcd'

#
#  Gather and save the metadata
#


parser = GenericParser(re_pattern,
                       group_names=captured_group_names,
                       group_transforms=captured_group_transforms,
                       top_fn='topology.pdb', step_ps=time_step)

meta = gather_metadata(os.path.join('traj', "*.{}".format(file_type)), parser)
save_meta(meta)