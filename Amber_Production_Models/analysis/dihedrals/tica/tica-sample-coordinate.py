"""Sample tICA coordinates

msmbuilder autogenerated template version 2
created 2017-05-23T16:38:49.125959
please cite msmbuilder in any publications


"""

import mdtraj as md

from msmbuilder.io.sampling import sample_dimension
from msmbuilder.io import load_trajs, save_generic, preload_top, backup

## Load
meta, ttrajs = load_trajs('ttrajs')

## Sample
inds = sample_dimension(ttrajs,
                        dimension=0,
                        n_frames=200, scheme='random')

save_generic(inds, "tica-dimension-0-inds.pickl")

## Make trajectory
top = preload_top(meta)

# Use loc because sample_dimension is nice
traj = md.join(
    md.load_frame(meta.loc[traj_i]['traj_fn'], index=frame_i, top=top)
    for traj_i, frame_i in inds
)

## Save
traj_fn = "tica-dimension-0.xtc"
backup(traj_fn)
traj.save(traj_fn)