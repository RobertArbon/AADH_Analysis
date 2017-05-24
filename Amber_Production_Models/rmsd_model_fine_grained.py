## Create a MSM using RMSD as a clustering metric with Ward clustering.
from msmbuilder.io import load_trajs, save_generic, load_meta, preload_tops
from msmbuilder.cluster import LandmarkAgglomerative
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV, ShuffleSplit
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

    # -------------------------------------------------------------------------
    # Load data and set parameters
    # -------------------------------------------------------------------------
    with Pool() as pool:
        trajs_dct = dict(pool.imap_unordered(traj_load, meta.iterrows()))
    trajs = [traj for traj in trajs_dct.values()]

    to_ns, t_max, frames_tot = get_timings(meta)
    n_clusters = int(np.sqrt(frames_tot))
    msm_lag = int(25/to_ns) # from eyeballing rmsd_timescales(-detail).png

    # -------------------------------------------------------------------------
    # Set pipeline
    # -------------------------------------------------------------------------
    clusterer = LandmarkAgglomerative(n_clusters=n_clusters,
                                      n_landmarks=n_clusters//10,
                                      linkage='ward', metric='rmsd',
                                      landmark_strategy='stride',
                                      random_state=None, max_landmarks=None,
                                      ward_predictor='ward')
    msm = MarkovStateModel(lag_time=msm_lag)
    pipe = Pipeline([('cluster', clusterer), ('msm', msm)])

    # -------------------------------------------------------------------------
    # Set param search object
    # -------------------------------------------------------------------------

    params = {'cluster__n_clusters': list((np.logspace(-0.5, 2, 10)*n_clusters).astype(int))}
    print(params)
    cv_iter = ShuffleSplit(n_splits=10, test_size=0.5)
    param_search = GridSearchCV(pipe, param_grid=params, cv=cv_iter)

    # -------------------------------------------------------------------------
    # Search param space and save
    # -------------------------------------------------------------------------

    param_search.fit(trajs)
    save_generic(param_search, 'models/rmsd_model.pickl')

    print('Best score of {0} was achieved with \n {1}'.format(
        param_search.best_score_, param_search.best_params_))




