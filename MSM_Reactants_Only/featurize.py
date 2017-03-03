#
# FEATURIZES DATA
#

from msmbuilder.featurizer import RawPositionsFeaturizer
from msmbuilder.preprocessing import RobustScaler
import numpy as np
import seaborn as sns
import mdtraj as md
from msmbuilder.io import load_meta, preload_tops, save_trajs, save_generic
from multiprocessing import Pool
import matplotlib
matplotlib.use('Agg')

sns.set_style('ticks')
colors = sns.color_palette('colorblind')

#
# Helper functions
#

## Box and whisker plot
def plot_box(ax, fxx):
    n_feats_plot = min(fxx.shape[1], 100)
    ax.boxplot(fxx[:, :100],
               boxprops={'color': colors[0]},
               whiskerprops={'color': colors[0]},
               capprops={'color': colors[0]},
               medianprops={'color': colors[2]},
               )

    if fxx.shape[1] > 100:
        ax.annotate("(Only showing the first 100 features)",
                    xy=(0.05, 0.95),
                    xycoords='axes fraction',
                    fontsize=14,
                    va='top',
                    )

    ax.set_xlabel("Feature Index", fontsize=16)
    xx = np.arange(0, n_feats_plot, 10)
    ax.set_xticks(xx)
    ax.set_xticklabels([str(x) for x in xx])
    ax.set_xlim((0, n_feats_plot + 1))
    ax.set_ylabel("Feature Value", fontsize=16)


def feat(irow):
    i, row = irow
    traj = md.load(row['traj_fn'], top=tops[row['top_fn']])
    feat_traj = featurizer.partial_transform(traj)
    return i, feat_traj


## Load
meta = load_meta()
tops = preload_tops(meta)
featurizer = RawPositionsFeaturizer()


## Do it in parallel
with Pool() as pool:
    feature_trajs = dict(pool.imap_unordered(feat, meta.iterrows()))

## Save
save_trajs(feature_trajs, 'ftrajs', meta)
save_generic(featurizer, 'featurizer.pickl')