# Choose features
from msmbuilder.io import load_trajs, save_trajs
import numpy as np
from multiprocessing import Pool
import matplotlib
matplotlib.use('Agg')
from matplotlib.pylab import plt
import sys
import seaborn as sns
from sklearn.neighbors.kde import KernelDensity
from scipy.signal import argrelextrema

# Don't prune these:
for feature in ['angles', 'bonds', 'contacts']:
    meta, ftraj = load_trajs('featurized_trajectories/{}-ftraj'.format(feature))
    save_trajs(ftraj, 'pruned_trajectories/{}-ftraj'.format(feature), meta)


# Prune these:

for feature in ['dihedrals']:
    meta, ftraj_dict = load_trajs('featurized_trajectories/{}-ftraj'.format(feature))
    ftraj = np.concatenate([traj for traj in ftraj_dict.values()])
    cos_idx = np.arange(0,ftraj.shape[1]-1, 2).reshape(-1,1)
    variance = ftraj[:, cos_idx].var(axis=0).reshape(-1,1)

    # Do KDE and split the data
    num_splits=3
    bandwidths = np.linspace(.01,.10,num=100)
    x = np.linspace(0,.5,1000).reshape(-1,1)
    for bw in bandwidths:
        kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(variance)
        dens = kde.score_samples(x)
        mi = argrelextrema(dens, np.less)[0]
        if mi.shape[0] == num_splits:
            break

    fig, ax = plt.subplots()
    ax.plot(x[:,0], dens, label="bw = {0:3.2f}".format(bw))
    ax.scatter(x[mi], dens[mi], label='Minima', marker='o')
    stride=100
    ax.set_xticks(x[::stride,0])
    ax.set_xticklabels(["{:2.3f}".format(xlab) for xlab in x[::stride,0]])
    ax.set_ylabel("log {} density".format(feature), fontsize=16)
    ax.set_xlabel("{} variance".format(feature), fontsize=16)
    plt.legend()
    plt.savefig("figures/{}-features-kde.pdf".format(feature))

    # Create mask for min variance features
    min_var = x[mi[0]]
    min_cos_idx = cos_idx[variance < min_var]
    max_cos_idx = cos_idx[variance >= min_var]
    min_sin_idx = min_cos_idx + 1
    max_sin_idx = max_cos_idx + 1
    min_idx = np.sort(np.append(min_cos_idx, min_sin_idx))
    max_idx = np.sort(np.append(max_cos_idx, max_sin_idx))

    # replot variance
    variance = ftraj.var(axis=0)
    fig, ax = plt.subplots()
    ax.bar(min_cos_idx, variance[min_cos_idx], color='g', align='center', label='Low variance', width=2)
    ax.bar(max_cos_idx, variance[max_cos_idx], color='r', align='center', label='High variance', width=2)
    stride=10
    # ax.set_xticks(cos_idx[::stride])    ax.set_ylabel("{} Variance".format(feature), fontsize=16)
    ax.set_xlabel("{} feature index".format(feature), fontsize=16)
    plt.legend()
    plt.savefig("figures/{}-features-high-low.pdf".format(feature))

    # Save pruned trajectory
    ptraj = {}
    for k, v in ftraj_dict.items():
        ptraj[k] = v[:, max_idx]

    save_trajs(ptraj, 'pruned_trajectories/{}-ftraj'.format(feature), meta)

