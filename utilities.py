import numpy as np
import seaborn as sns
import mdtraj as md
import re
from itertools import combinations
import matplotlib
matplotlib.use('Agg')
from matplotlib.pylab import plt
import pandas as pd
import seaborn as sns
import mdtraj as md
from msmbuilder.io.sampling import sample_dimension
from msmbuilder.io import load_trajs, save_generic, preload_top, backup
from pyemma.coordinates import featurizer, source
sns.set_style('ticks')
colors = sns.color_palette('colorblind')

## Box and whisker plot
# def plot_box(ax, fxx, feature_name='Feature'):
#     n_feats_plot = min(fxx.shape[1], 100)
#     ax.boxplot(fxx[:, :100],
#                boxprops={'color': colors[0]},
#                whiskerprops={'color': colors[0]},
#                capprops={'color': colors[0]},
#                medianprops={'color': colors[2]},
#                )
#
#     if fxx.shape[1] > 100:
#         ax.annotate("(Only showing the first 100 features)",
#                     xy=(0.05, 0.95),
#                     xycoords='axes fraction',
#                     fontsize=14,
#                     va='top',
#                     )
#
#     ax.set_xlabel("Feature Index", fontsize=16)
#     xx = np.arange(0, n_feats_plot, 10)
#     ax.set_xticks(xx)
#     ax.set_xticklabels([str(x) for x in xx])
#     ax.set_xlim((0, n_feats_plot + 1))
#     ax.set_ylabel("{} Value".format(feature_name), fontsize=16)

def plot_box(axes, fxx, n_feats_plot, nrows,  feature_name='Feature'):
    ndim = fxx.shape[1]
    if nrows == 1:
        axes.boxplot(fxx,
                   boxprops={'color': colors[0]},
                   whiskerprops={'color': colors[0]},
                   capprops={'color': colors[0]},
                   medianprops={'color': colors[2]},
                   )
        axes.set_xlabel("Feature Index", fontsize=16)
        xx = np.arange(0, n_feats_plot, 10)
        axes.set_xticks(xx)
        axes.set_xticklabels([str(x) for x in xx])
        axes.set_xlim((0, n_feats_plot + 1))
        axes.set_ylabel("{} Value".format(feature_name), fontsize=16)

    elif nrows > 1:
        for idx, ax in enumerate(axes):
            start = idx*n_feats_plot
            end = min((idx+1)*n_feats_plot, ndim)




    # axes.set_xlabel("Feature Index", fontsize=16)
    # xx = np.arange(0, n_feats_plot, 10)
    # axes.set_xticks(xx)
    # axes.set_xticklabels([str(x) for x in xx])
    # axes.set_xlim((0, n_feats_plot + 1))
    # axes.set_ylabel("{} Value".format(feature_name), fontsize=16)


def msmb_feat(args):
    """
    featurizes with MSM builder featurizers.
    :param args:
    :return:
    """
    irow, feat, tops = args
    i, row = irow
    traj = md.load(row['traj_fn'], top=tops[row['top_fn']])
    feat_traj = feat.partial_transform(traj)
    return i, feat_traj


def pyemma_feat(args):
    irow, featurizer_name, tops, indices = args
    i, row = irow
    traj, top = row['traj_fn'], tops[row['top_fn']]
    feat = featurizer(top)
    try:
        adder = getattr(feat, featurizer_name)
        adder(indexes=indices, cossin=True)
        feat_traj = np.squeeze(source(traj, features=feat).get_output(), axis=0)
        return i, feat_traj
    except AttributeError:
        print("pyEMMA doesn't have {} as a featurizer".format(featurizer_name))


def plot_param_line(df, axes, params=None):
    """
    Plots parameter scores.  Dimension of axes must match size of number of
    combinations.
    :param df: output from sk-learn CV search algo.
    :return: None
    """
    if params is None:
        params = [x for x in df.columns if re.search('param_(.*)', x)]
    par_combs = list(combinations(params, 2))
    for idx, (p, q) in enumerate(par_combs):
        p_vals = df[p].unique()
        for idx, p_val in enumerate(p_vals):
            cols = sns.color_palette('Blues', n_colors=len(p_vals))
            x = df.ix[df[p] == p_val, q].values
            y = df.ix[df[p] == p_val, 'mean_test_score'].values
            xlab = re.search('param_(.*)', q).group(1)
            ylab = 'score'
            zlab = '{0} = {1:4.2f}'.format(re.search('param_(.*)', p).group(1), p_val)
            if len(par_combs) > 1:
                ax = axes[idx]
            else:
                ax = axes
            ax.plot(x, y, marker='o', label=zlab, c=cols[idx])
            ax.set_xlabel(xlab)
            ax.set_ylabel(ylab)
            ax.legend()
    return axes

def get_param_combs(df):
    """
    gets list of parameters from CV search DF.  Also returns number of combinations
    :param df:
    :return:
    """
    params = [x for x in df.columns if re.search('param_(.*)', x)]
    n_params = len(params)
    n_comb = int(0.5 * n_params * (n_params - 1))
    return params, n_comb


def plot_tica_distribution(txx, sample_size=100, ndims=4):
    """
    plots pairgrid of tICa
    :param txx:
    :param sample_size:
    :param ndims:
    :return:
    """
    straj = txx[np.random.choice(txx.shape[0], size=sample_size, replace=False), :ndims]
    df = pd.DataFrame(straj)
    df.columns = ['tIC {}'.format(int(x) + 1) for x in df.columns]
    g = sns.PairGrid(df, diag_sharey=True)
    g.map_diag(sns.kdeplot, lw=3)
    g.map_lower(sns.kdeplot, shade=True, shade_lowest=False, cmap="Blues_d", )
    g.map_upper(plt.scatter, alpha=0.1)
    return g


## Timescales
def plot_timescales(ax, meta, tica, units='ps'):
    timestep = meta['step_ps'].unique()
    assert len(timestep) == 1, timestep
    timestep = float(timestep[0])  # ps
    to_us = (
        (1.0 / 1000)  # ps -> ns
        * (1.0 / 1000)  # ns -> us
        * (timestep / 1)  # steps -> ps
    )
    if units=='ps':
        fac = 1e6
        lab = 'ps'
    elif units=='ns':
        fac = 1e3
        lab = 'ns'
    elif units=='us':
        fac = 1
        lab = '\mu s'
    else:
        raise ValueError('units should be ps, ns, us')

    ax.hlines(tica.timescales_ * to_us*fac,
              0, 1,
              color=colors[0])
    ax.set_ylabel(r'Timescales / $\mathrm{'+lab+'}$', fontsize=18)
    ax.set_xticks([])
    ax.set_xlim((0, 1))
    return(ax)


def plot_eigenvectors(axes, tica):
    vectors = tica.eigenvectors_.T
    width = 1
    feat_idx = np.arange(vectors.shape[1])
    for idx, ax in enumerate(axes):
        ax.bar(feat_idx, vectors[idx], align='center')
        ax.set_ylabel('EigVec\ncomponent')
    ax.set_xlabel('Feature index')
    return axes


# Sample trajectories
def sample_tica_dim(dim=0, n_frames=200,meta=None, ttrajs=None):

    ## Load
    if (not meta is None) & (not ttrajs is None):

        ## Sample
        # These are apparently ordered according tica value
        inds = sample_dimension(ttrajs,
                                dimension=dim,
                                n_frames=n_frames, scheme='random')

        save_generic(inds, "tica-dimension-{}-inds.pickl".format(dim+1))

        ## Get tica components
        tica_values = np.array([ttrajs[traj_i][frame_i][dim] for traj_i, frame_i in inds])
        tica_values = (tica_values - tica_values.min())/(tica_values.max()-tica_values.min())
        tica_values *= 10
        ## Make trajectory
        top = preload_top(meta)

        # Use loc because sample_dimension is nice
        traj = md.join(
            md.load_frame(meta.loc[traj_i]['traj_fn'], index=frame_i, top=top)
            for traj_i, frame_i in inds
        )

        ## Supperpose

        ## Save
        traj_fn = "tica-dimension-{}.dcd".format(dim+1)
        backup(traj_fn)
        traj.save(traj_fn)
    else:
        raise ValueError('Specify meta data and trajectory objects')


def to_dataframe(traj_dict, nframes, dt):

    id_cols = ['PDB', 'Temp', 'Pres', 'Prod_ID', 'Site_ID']
    df = pd.DataFrame.from_records(data=traj_dict)
    df['Time_ps'] = np.arange(nframes) * dt
    df = pd.melt(df, id_vars=['Time_ps'], var_name='Key', value_name='Trajectory')
    df[id_cols] = pd.DataFrame(df['Key'].tolist())
    # del df['Production_ID']
    return df