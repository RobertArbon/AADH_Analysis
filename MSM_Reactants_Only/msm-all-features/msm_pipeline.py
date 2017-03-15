# The pipeline for constructing an MSM
from msmbuilder.io import load_trajs, save_generic
from msmbuilder.preprocessing import StandardScaler
from msmbuilder.decomposition import tICA
from msmbuilder.cluster import MiniBatchKMeans
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
import scipy


# Load all features
# features = ['angles', 'dihedrals', 'bonds', 'contacts']
features = ['dihedrals', 'contacts']

all_trajs = []
for feature in features:
    meta, ftraj = load_trajs('pruned_trajectories/{}-ftraj'.format(feature))
    ftraj = [traj for traj in ftraj.values()]
    all_trajs.append(ftraj)
all_trajs = list(map(list, zip(*all_trajs)))

# Concatenate features
# ftraj = num trajectories x np.array(n_frames, n_features)
ftraj = []
for traj in all_trajs:
    tmp = []
    for feat in traj:
        if feat.ndim == 1:
            feat = feat.reshape(-1, 1)
        tmp.append(feat)
    ftraj.append(np.concatenate(tmp, axis=1))

# Make Pipeline
cv_iter = ShuffleSplit(n_splits=5, test_size=0.5)
estimators = [('scale',StandardScaler()),('tica', tICA()), ('cluster', MiniBatchKMeans(random_state=0)),
              ('msm', MarkovStateModel())]
param_grid = {'cluster__n_clusters': list(np.linspace(200, 500, num=2).astype(int)),
              'tica__n_components': list(np.linspace(10, 30, num=2).astype(int)),
              'tica__lag_time': list(np.linspace(200, 500, num=2).astype(int))}

params = {'cluster__n_clusters': scipy.stats.randint(low=200,high=700),
              'tica__n_components':  scipy.stats.randint(low=2,high=40),
              'tica__lag_time':  scipy.stats.randint(low=100,high=800)}


pipe = Pipeline(estimators)
pipe.set_params(msm__lag_time=500)
pipe.set_params(msm__n_timescales=10)

if __name__ == "__main__":

    cvSearch = RandomizedSearchCV(pipe, params, n_jobs=3, verbose=1, cv=cv_iter, n_iter=100)

    print("Performing grid search...")
    print("pipeline:", [name for name, _ in pipe.steps])
    print("parameters:")
    print(params)
    t0 = time()
    cvSearch.fit(ftraj)
    print("done in %0.3fs" % (time() - t0))
    print()

    print("Best score: %0.3f" % cvSearch.best_score_)
    print("Best parameters set:")
    best_parameters = cvSearch.best_estimator_.get_params()
    for param_name in sorted(params.keys()):
        print("\t%s: %r" % (param_name, best_parameters[param_name]))

    df = pd.DataFrame(cvSearch.cv_results_)
    save_generic(df, 'results/20-40-Clusters.pickl')


