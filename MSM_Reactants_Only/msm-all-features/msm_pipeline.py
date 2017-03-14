# The pipeline for constructing an MSM
from msmbuilder.io import load_trajs, save_generic
from msmbuilder.decomposition import tICA
from msmbuilder.cluster import KMeans, MiniBatchKMeans, LandmarkAgglomerative
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

# Load all features
features = ['angles', 'dihedrals', 'bonds', 'contacts']
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
cv_iter = ShuffleSplit(n_splits=5, test_size=0.1)
estimators = [('tica', tICA()), ('cluster', MiniBatchKMeans(random_state=0)), ('msm', MarkovStateModel())]
param_grid = {'cluster__n_clusters': list(np.linspace(200, 500, num=2).astype(int)),
              'tica__n_components': list(np.linspace(10, 30, num=2).astype(int)),
              'tica__lag_time': list(np.linspace(200, 500, num=2).astype(int))}

pipe = Pipeline(estimators)
pipe.set_params(msm__lag_time=500)

if __name__ == "__main__":

    grid_search = GridSearchCV(pipe, param_grid, n_jobs=-1, verbose=1, cv=cv_iter)

    print("Performing grid search...")
    print("pipeline:", [name for name, _ in pipe.steps])
    print("parameters:")
    print(param_grid)
    t0 = time()
    grid_search.fit(ftraj)
    print("done in %0.3fs" % (time() - t0))
    print()

    print("Best score: %0.3f" % grid_search.best_score_)
    print("Best parameters set:")
    best_parameters = grid_search.best_estimator_.get_params()
    for param_name in sorted(param_grid.keys()):
        print("\t%s: %r" % (param_name, best_parameters[param_name]))

    df = pd.DataFrame(grid_search.cv_results_)
    save_generic(df, 'results/grid_search.pickl')
    # # Print
    # timescales = msm.timescales_.T
    # fig, ax = plt.subplots()
    # ax.hlines(timescales, xmin=0, xmax=1, label='Lag {}'.format(lag))
    # plt.legend()
    # ax.set_xlim(0,1)
    # ax.set_ylim(0,10000)
    # plt.savefig('figures/timescales-clusters-{}.png'.format(num_clusters))

