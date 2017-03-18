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
all_trajs = [traj for traj in traj_dict.values() ]

cluster = LandmarkAgglomerative(n_clusters=200, n_landmarks=int(totframes/100), linkage='ward', metric='rmsd')
cluster.fit(all_trajs)

all_ctrajs_dict = {}
count = 0
for k, v in traj_dict.items():
    count += 1
    print(count)
    all_ctrajs_dict[k] = cluster.transform(v)

long_ctrajs = [traj for traj in all_ctrajs_dict.values() if traj.shape[0] > 1000]
all_ctrajs = [traj for traj in all_ctrajs_dict.values()]

save_generic(cluster, 'cluster-200')
save_trajs(all_ctrajs_dict, 'ctraj-200', meta)

if __name__ == "__main__":

    lags = np.concatenate((np.arange(200, 1000, 200),np.arange(1000, 5000, 500)))
    all_msms = []

    for lag in lags:
        print('Fitting lag {}'.format(lag))
        if lag > 1000:
            trajs = long_ctrajs
        else:
            trajs = all_ctrajs
        msm = MarkovStateModel(lag_time=int(lag), n_timescales=100)
        msm.fit(trajs)
        all_msms.append(msm)
    save_generic(all_msms, 'rmsd_msms.pickl')