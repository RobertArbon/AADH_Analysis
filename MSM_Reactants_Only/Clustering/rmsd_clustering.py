from msmbuilder.io import load_meta, preload_tops, save_generic, itertrajs, backup
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
ref = md.load('topology.pdb')

def traj_load(irow):
    i, row = irow
    traj = md.load(row['traj_fn'], top=tops[row['top_fn']])
    return i, traj

traj_dict = dict(map(traj_load, meta.iterrows()))
trajs = [traj for traj in traj_dict.values()]

# cluster
num_clusters=10
cluster = LandmarkAgglomerative(n_clusters=num_clusters, n_landmarks=200, linkage='ward', metric='rmsd')
cluster.fit(trajs)

ctraj = {}
for k, v in traj_dict.items():
    v = cluster.partial_predict(v)
    diff = nframes-v.shape[0]
    v = np.append(v, np.zeros(diff)-1)
    ctraj[k] = v

# Convert to DF for plotting and sampling.
df = to_dataframe(ctraj, nframes, dt=1)

# Plot trajectories
sample = df.sample(frac=0.1, axis=0)
sample.sort_values(by=['Prod_ID', 'Site_ID', 'Time_ps'], inplace=True)
g = sns.FacetGrid(sample, col='Prod_ID',hue='Site_ID', col_wrap=10)
g.map(plt.scatter, 'Time_ps', 'Trajectory', alpha=0.5)
g.set(ylim=(-1.5,10))
g.fig.tight_layout()
plt.savefig('figures/rmsd_cluster_trajectory.pdf')

# Save dataframe
save_generic(df, 'clusters/rmsd_cluster_trajectory.pickl')

# Sampling
for i in range(num_clusters):
    df_smp = df.ix[df['Trajectory']==i, ['Key', 'Time_ps']].sample(100)
    inds = zip(df_smp['Key'], df_smp['Time_ps'])

    # Use loc because sample_dimension is nice
    traj = md.join(
        md.load_frame(meta.loc[traj_i]['traj_fn'], index=frame_i, top=meta.loc[traj_i]['top_fn'])
        for traj_i, frame_i in inds
    )

    # Original trajectories include both BT1 and BT2 so need to superpose
    traj.superpose(reference=ref)

    # Save
    traj_fn = "clusters/rmsd_cluster-{}.dcd".format(i)
    backup(traj_fn)
    traj.save(traj_fn)