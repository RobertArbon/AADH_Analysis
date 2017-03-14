# Plot features
from msmbuilder.io import load_trajs
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.pylab import plt
import sys
import seaborn as sns
colors = sns.color_palette("colorblind", 8)

for feature in ['angles', 'dihedrals', 'bonds', 'contacts']:
    meta, ftraj = load_trajs('featurized_trajectories/{}-ftraj'.format(feature))
    ftraj = np.concatenate([traj for traj in ftraj.values()])

    if feature in ['angles', 'dihedrals']:
        sample = ftraj[np.random.choice(ftraj.shape[0], size=10000), :]
        sample = sample[:,np.arange(0,ftraj.shape[1], 2)]
        print(feature, sample.shape)
    elif feature in ['contacts']:
        sample = ftraj[np.random.choice(ftraj.shape[0], size=10000)]
        print(feature, sample.shape)
    else:
        sample = ftraj[np.random.choice(ftraj.shape[0], size=10000), :]
        print(feature, sample.shape)

    try:
        n_feats_plot = sample.shape[1]
    except IndexError:
        n_feat_plot = 1

    # BOX PLOT
    fig, axes = plt.subplots()

    axes.boxplot(sample,
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
    axes.set_ylabel("{} Value".format(feature), fontsize=16)
    plt.savefig("figures/{}-features-box.pdf".format(feature))

    variance = sample.var(axis=0)
    fig, axes = plt.subplots()
    try:
        index = np.arange(variance.shape[0])
    except IndexError:
        index = np.array([1])

    axes.bar(index, variance, align='center')
    axes.set_xlabel("Feature Index", fontsize=16)
    xx = np.arange(0, n_feats_plot, 10)
    axes.set_xticks(xx)
    axes.set_xticklabels([str(x) for x in xx])
    axes.set_xlim((-0.5, n_feats_plot + 0.5))
    axes.set_ylabel("{} Variance".format(feature), fontsize=16)
    plt.savefig("figures/{}-features-var.pdf".format(feature))
    print('Finished {}'.format(feature))
