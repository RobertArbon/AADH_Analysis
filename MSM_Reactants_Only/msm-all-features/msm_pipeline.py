# The pipeline for constructing an MSM
from msmbuilder.io import load_trajs
from msmbuilder.decomposition import tICA
from msmbuilder.cluster import KMeans, MiniBatchKMeans
from msmbuilder.msm import MarkovStateModel
import numpy as np
from multiprocessing import Pool
import matplotlib
matplotlib.use('Agg')
from matplotlib.pylab import plt
import sys

if __name__ == "__main__":
    features = ['angles', 'dihedrals', 'bonds']
    all_trajs = []
    for feature in features:
        meta, ftraj = load_trajs('featurized_trajectories/{}-ftraj'.format(feature))
        ftraj = [traj for traj in ftraj.values()]
        all_trajs.append(ftraj)
    all_trajs = map(list, zip(*all_trajs))
    ftraj = []
    for traj in all_trajs:
        ftraj.append(np.concatenate([feat for feat in traj], axis=1))

    # Tica
    print('Peforming tICA decomposition')
    num_comp = 50
    tica_lag = 2500
    decomp = tICA(n_components=num_comp, lag_time=tica_lag, kinetic_mapping=True)
    ttraj = decomp.fit_transform(ftraj)

    # Clustering
    print('Peforming clustering')
    num_clusters = 200
    cluster = MiniBatchKMeans(n_clusters=num_clusters, random_state=0)
    ctraj = cluster.fit_transform(ttraj)

    # MSM
    print('Estimating MSM')
    lags = np.linspace(2000, 4000, 10).astype(int)
    msms = []
    for lag in lags:
        msm = MarkovStateModel(lag_time=lag)
        msm.fit_transform(ctraj)
        msms.append(msm)

    num_ts = 20
    timescales = np.array([msm.timescales_[:num_ts] for msm in msms])
    timescales = timescales.T
    fig, ax = plt.subplots()
    for i in range(num_ts):
        ax.plot(lags, timescales[i])
    ax.set_ylim(0,10000)
    plt.savefig('figures/timescales.png')

