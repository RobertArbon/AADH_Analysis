from msmbuilder.io import load_meta, preload_tops, save_generic
import mdtraj as md
from msmbuilder.featurizer import RMSDFeaturizer
import matplotlib
matplotlib.use('Agg')
from matplotlib.pylab import plt
from multiprocessing import Pool
from utilities import msmb_feat
import numpy as np
import seaborn as sns
import pandas as pd
import sys

sns.set_palette("Spectral", 50)

meta = load_meta()
tops = preload_tops(meta)
legends = []
with sns.plotting_context("notebook", font_scale=1.5):

    for frame in np.arange(0, 1000, 200):
        print(frame)
        traj = md.join(
            md.load_frame(row['traj_fn'], index=frame, top=tops[row['top_fn']])
            for i, row in meta.iterrows()
        )
        traj = traj.superpose(reference=traj)
        distances = np.empty((traj.n_frames, traj.n_frames))
        for i in range(traj.n_frames):
            distances[i] = md.rmsd(traj, traj, i)

        distances = np.tril(distances, k=0) + np.triu(np.ones(distances.shape), k=0)*(-1)
        dist = distances.reshape(-1,1)
        dist = dist[dist>=0]*10 # to angstroms

        g = sns.distplot(dist, hist=False, kde=True, rug=False, label='t = {} ps'.format(frame))

    for frame in np.arange(1000, 5000, 200):
        print(frame)
        traj = md.join(
            md.load_frame(row['traj_fn'], index=frame, top=tops[row['top_fn']])
            for i, row in meta.iterrows() if row['nframes'] > 1000
        )
        traj = traj.superpose(reference=traj)
        distances = np.empty((traj.n_frames, traj.n_frames))
        for i in range(traj.n_frames):
            distances[i] = md.rmsd(traj, traj, i)

        distances = np.tril(distances, k=0) + np.triu(np.ones(distances.shape), k=0)*(-1)
        dist = distances.reshape(-1,1)
        dist = dist[dist>=0]*10 # to angstroms

        g = sns.distplot(dist, hist=False, kde=True, rug=False, label='t = {} ps'.format(frame))

    handels, labels = g.get_legend_handles_labels()
    handels = [handels[0], handels[-1]]
    labels = [labels[0], labels[-1]]
    plt.legend(handles=handels, labels=labels)
    plt.xlabel('RMSD $\AA$')
    plt.ylabel('Frequency')
    plt.title('KDE of sampling distribution')
    plt.tight_layout()
    plt.savefig('KDE_sampling_dist.png', transparent=True)