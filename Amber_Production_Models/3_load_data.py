from msmbuilder.io import GenericParser, save_meta, gather_metadata
from os.path import join
#
# File name parsing and metadata
#
def identity(x):
    return x

re_pattern = '(\w+)-([0-9]+)-as([0-9]+)*'
captured_group_names = ['PDB', 'Traj_Num', 'Act_Site']
captured_group_transforms = [identity, int, int]
time_step = 10 #10 ps
file_type = 'nc'

#
#  Gather and save the metadata
#


parser = GenericParser(re_pattern,
                       group_names=captured_group_names,
                       group_transforms=captured_group_transforms,
                       top_fn='proc_traj/2agy-as1.prmtop', step_ps=time_step)

meta = gather_metadata(join('proc_traj', "*.{}".format(file_type)), parser)
save_meta(meta)