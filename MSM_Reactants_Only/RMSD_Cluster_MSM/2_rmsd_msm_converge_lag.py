# The pipeline for constructing an MSM
from msmbuilder.io import load_trajs, save_generic, load_meta, preload_tops, save_trajs, load_generic
from msmbuilder.cluster import LandmarkAgglomerative
from msmbuilder.msm import MarkovStateModel
import numpy as np
import mdtraj as md
from multiprocessing import Pool
from os.path import isdir, isfile
import sys

def clust(args):
    k, v, cluster = args
    print(k)
    ctraj = cluster.transform(v)
    return k, ctraj

if __name__ == "__main__":

    # Load data
    meta = load_meta()
    tops = preload_tops(meta)
    totframes = meta['nframes'].sum()

    ctraj_path = 'ctraj-200'
    if isdir(ctraj_path):
        meta, all_ctrajs_dict = load_trajs(ctraj_path)
    else:

        def traj_load(irow):
            i, row = irow
            traj = md.load(row['traj_fn'], top=tops[row['top_fn']])
            return i, traj


        traj_dict = dict(map(traj_load, meta.iterrows()))
        all_trajs = [traj for traj in traj_dict.values()]

        cluster = LandmarkAgglomerative(n_clusters=200, n_landmarks=int(totframes /200), linkage='ward', metric='rmsd')
        cluster.fit(all_trajs)
        # TODO will this work?
        args = [(k,v,cluster) for k, v in traj_dict.items()]

        with Pool() as pool:
            all_ctrajs_dict = dict(pool.imap_unordered(clust, args))

        save_generic(cluster, 'cluster-200')
        save_trajs(all_ctrajs_dict, 'ctraj-200', meta)

    long_ctrajs = [np.squeeze(traj) for traj in all_ctrajs_dict.values() if traj.shape[0] > 1000]
    all_ctrajs = [np.squeeze(traj) for traj in all_ctrajs_dict.values()]

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