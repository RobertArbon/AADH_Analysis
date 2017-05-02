# The pipeline for constructing an MSM
from msmbuilder.io import load_trajs, save_generic, load_meta, preload_tops
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
trajs = [traj for traj in traj_dict.values()]


# Make Pipeline
cv_iter = ShuffleSplit(n_splits=5, test_size=0.5)
estimators = [('cluster', LandmarkAgglomerative(n_clusters=2, n_landmarks=int(totframes/200), linkage='ward', metric='rmsd')),
              ('msm', MarkovStateModel())]

params = {'cluster__n_clusters': [200]}

pipe = Pipeline(estimators)
pipe.set_params(msm__lag_time=999)
pipe.set_params(msm__n_timescales=20)

if __name__ == "__main__":

    cvSearch = GridSearchCV(pipe, params, n_jobs=1, verbose=1, cv=cv_iter)

    print("Performing grid search...")
    print("pipeline:", [name for name, _ in pipe.steps])
    print("parameters:")
    print(params)
    t0 = time()
    cvSearch.fit(trajs)
    print("done in %0.3fs" % (time() - t0))
    print()

    print("Best score: %0.3f" % cvSearch.best_score_)
    print("Best parameters set:")
    best_parameters = cvSearch.best_estimator_.get_params()
    for param_name in sorted(params.keys()):
        print("\t%s: %r" % (param_name, best_parameters[param_name]))

    df = pd.DataFrame(cvSearch.cv_results_)
    save_generic(df, 'results/lag999-ncluster200.pickl')