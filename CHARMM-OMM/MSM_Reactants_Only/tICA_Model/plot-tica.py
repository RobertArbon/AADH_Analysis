from msmbuilder.io import load_generic, load_trajs
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from utilities import *
import numpy as np
import pandas as pd
import seaborn as sns

# LOAD DATA
# TODO change the default name of the param search
param_search = load_generic('Positions-grid-search-results.pickl')
tica = load_generic('Positions-tica.pickl')
meta, ttrajs = load_trajs('Positions-ttrajs')
txx = np.concatenate(list(ttrajs.values()))
params, n_comb = get_param_combs(param_search)

# # PARAMETER SEARCH PLOT
# fig, axes = plt.subplots(nrows=n_comb, ncols=1)
# axes = plot_param_line(param_search, axes, params)
# plt.tight_layout()
# plt.savefig('Positions-param-results.pdf')
# plt.clf()
#
# # tICA DISTRIBUTION PLOT
# plot_tica_distribution(txx, sample_size=2000, ndims=4)
# plt.savefig('Positions-tica-dist.pdf')
# plt.clf()
#
# # TIME SCALES PLOT
# fig, axes = plt.subplots()
# axes = plot_timescales(axes, meta, tica, units='ps')
# plt.tight_layout()
# plt.savefig('Positions-tica-timescales.pdf')
# plt.clf()

# # PLOT EIGENVECTORS
# ndim = 4
# fig, axes = plt.subplots(nrows=ndim, sharex=True, sharey=True)
# axes = plot_eigenvectors(axes, tica)
# plt.savefig('Positions-tica-eigenvectors.pdf')

# # SAMPLE TRAJECTORY
# for i in range(4):
#     sample_tica_dim(dim=i, n_frames=200, meta=meta, ttrajs=ttrajs)
i = 0
inds = sample_dimension(ttrajs,
                        dimension=i,
                        n_frames=100, scheme='random')

tica_values = np.array([ttrajs[traj_i][frame_i][i] for traj_i, frame_i in inds])
tica_values = (tica_values - tica_values.min())/(tica_values.max() - tica_values.min())

## Make trajectory
top = preload_top(meta)
bfactors = tica_values[:, np.newaxis]*np.ones((tica_values.shape[0], top.n_atoms))


# Use loc because sample_dimension is nice
traj = md.join(
    md.load_frame(meta.loc[traj_i]['traj_fn'], index=frame_i, top=top)
    for traj_i, frame_i in inds
)

## Save
traj_fn = "tica-dimension-{}.pdb".format(i + 1)
backup(traj_fn)
traj.save_pdb(traj_fn, bfactors=bfactors)