# The pipeline for constructing an MSM
from msmbuilder.io import load_trajs, save_generic, load_meta, preload_tops, save_trajs
from msmbuilder.cluster import LandmarkAgglomerative
from msmbuilder.msm import MarkovStateModel
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.pylab import plt
import sys
from sklearn.pipeline import Pipeline
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV, ShuffleSplit
from time import time
import pandas as pd
import mdtraj as md


# Load data
meta = load_meta()
tops = preload_tops(meta)
totframes = meta['nframes'].sum()

def traj_load(irow):
    i, row = irow
    traj = md.load(row['traj_fn'], top=tops[row['top_fn']])
    return i, traj

traj_dict = dict(map(traj_load, meta.iterrows()))
trajs = [traj for traj in traj_dict.values() if traj.n_frames > 1000]
print(len(trajs))
num_clust = 20
cluster = LandmarkAgglomerative(n_clusters=num_clust, n_landmarks=int(totframes/100), linkage='ward', metric='rmsd')
ctrajs = cluster.fit_transform(trajs)

# print('Fitting cluster labels for MSM')
# ctraj = {}
# count = 0
# for k, v in traj_dict.items():
#     print(k, count)
#     count +=1
#     ctraj[k] = cluster.partial_predict(v)
#
# ctrajs = [traj for traj in ctraj.values() if traj.shape[0] > 1000]

print('Fitting MSM')
lag = 4000
msm = MarkovStateModel(lag_time=lag, n_timescales=50)
msm.fit(ctrajs)

# save_trajs(ctraj, 'results/nclusters-{0}-ctraj'.format(num_clust), meta)
save_generic(cluster, 'results/clusterer-nclusters-{0}.pickle'.format(num_clust))
save_generic(msm, 'results/msm-lag-{0}-nclusters-{1}.pickl'.format(lag, num_clust))
