## Create a MSM using RMSD as a clustering metric with Ward clustering.
from msmbuilder.io import load_trajs, save_generic, load_meta, preload_tops
from msmbuilder.cluster import LandmarkAgglomerative
from sklearn.pipeline import Pipeline
# from sklearn.model_selection import GridSearchCV, ShuffleSplit
import mdtraj as md
from multiprocessing import Pool
import numpy as np
from msmbuilder.msm import MarkovStateModel
import matplotlib.pyplot as plt

# Load data
meta = load_meta()
tops = preload_tops(meta)


def get_timings(meta):
    frames_tot = meta['nframes'].sum()
    n_frames = meta['nframes'].unique()
    assert (len(n_frames) == 1, 'Different trajectory lengths')
    n_frames = n_frames[0]
    dt = meta['step_ps'][0]
    to_ns = dt/1000
    t_max = n_frames*to_ns
    return to_ns, t_max, frames_tot


def traj_load(irow):
    i, row = irow
    traj = md.load(row['traj_fn'], top=tops[row['top_fn']])
    return i, traj

if __name__ == '__main__':

    with Pool() as pool:
        trajs_dct = dict(pool.imap_unordered(traj_load, meta.iterrows()))
    trajs = [traj for traj in trajs_dct.values()]

    to_ns, t_max, frames_tot = get_timings(meta)
    n_clusters = int(np.sqrt(frames_tot))
    print(n_clusters)
    # n_clusters = int(frames_tot/1000)
    clusterer = LandmarkAgglomerative(n_clusters=n_clusters,
                                      n_landmarks=n_clusters//10,
                                      linkage='ward', metric='rmsd',
                                      landmark_strategy='stride',
                                      random_state=None, max_landmarks=None,
                                      ward_predictor='ward')
    ctrajs = clusterer.fit_transform(trajs)

    lags = (np.arange(1,50,1)/to_ns).astype(int)
    n_timescales=50
    timescales = np.zeros((lags.shape[0], n_timescales))
    eigenvalues = np.zeros((lags.shape[0], n_timescales))

    for idx, lag in enumerate(lags):
        msm = MarkovStateModel(lag_time=lag, n_timescales=n_timescales)
        msm.fit_transform(ctrajs)
        timescales[idx] = msm.timescales_
        eigenvalues[idx] = msm.eigenvalues_[1:]



    for idx in range(n_timescales):
        plt.plot(lags*to_ns, timescales.T[idx])
    plt.savefig('figures/rmsd_timescales.png')
    plt.ylim((0,int(np.max(timescales.T[1]))))
    plt.savefig('figures/rmsd_timescales-detail.png')
    plt.clf()

    for idx in range(n_timescales):
        plt.plot(lags*to_ns, eigenvalues.T[idx])
    plt.savefig('figures/rmsd_eigenvalues.png')
    # Make Pipeline


