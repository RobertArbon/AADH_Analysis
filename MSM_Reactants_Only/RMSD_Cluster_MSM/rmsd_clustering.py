from msmbuilder.io import load_meta, preload_tops, save_generic, itertrajs, backup, save_trajs
import mdtraj as md
from msmbuilder.cluster import LandmarkAgglomerative
from msmbuilder.featurizer import RMSDFeaturizer
import matplotlib
matplotlib.use('Agg')
from matplotlib.pylab import plt
from multiprocessing import Pool
from utilities import msmb_feat
import numpy as np
import seaborn as sns
import pandas as pd
import sys
from utilities import to_dataframe


# load trajectories
meta = load_meta()
tops = preload_tops(meta)
nframes = int(np.max(meta['nframes'].unique()[0]))
totframes = meta['nframes'].sum()
print(totframes)
ref = md.load('topology.pdb')

def traj_load(irow):
    i, row = irow
    traj = md.load(row['traj_fn'], top=tops[row['top_fn']])
    return i, traj

traj_dict = dict(map(traj_load, meta.iterrows()))
trajs = [traj for traj in traj_dict.values()]

# cluster
print('Attempting to cluster')
num_clusters=20
cluster = LandmarkAgglomerative(n_clusters=num_clusters, n_landmarks=int(totframes/100), linkage='ward', metric='rmsd')
cluster.fit(trajs)

#
# print('Fitting cluster labels')
# ctraj = {}
# for k, v in traj_dict.items():
#     v = cluster.partial_predict(v)
#     diff = nframes-v.shape[0]
#     v = np.append(v, np.zeros(diff)-1)
#     ctraj[k] = v

# Convert to DF for plotting and sampling.
# df = to_dataframe(ctraj, nframes, dt=1)

print('Fitting cluster labels for MSM')
ctraj = {}
for k, v in traj_dict.items():
    ctraj[k] = cluster.partial_predict(v)



# Save dataframe
save_generic(df, 'clusters/rmsd_cluster_trajectory.pickl')
save_trajs(ctraj, 'ftraj', meta)
