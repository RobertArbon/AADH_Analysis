from msmbuilder.io import load_meta, preload_tops, save_generic, itertrajs, backup, load_trajs
import mdtraj as md
from msmbuilder.cluster import LandmarkAgglomerative
import matplotlib
matplotlib.use('Agg')
from matplotlib.pylab import plt
import numpy as np
import seaborn as sns
from utilities import to_dataframe


# load trajectories
feature = 'dihedrals'
meta, traj_dict= load_trajs('pruned_trajectories/{}-ftraj'.format(feature))
trajs = [traj for traj in traj_dict.values()]
nframes = int(np.max(meta['nframes'].unique()[0]))

# cluster
num_clusters=10
cluster = LandmarkAgglomerative(n_clusters=num_clusters, n_landmarks=200, linkage='ward', metric='euclidean')
cluster.fit(trajs)

ctraj = {}
for k, v in traj_dict.items():
    v = v.copy(order='C')
    v = cluster.partial_predict(v)
    diff = nframes-v.shape[0]
    v = np.append(v, np.zeros(diff)-1)
    ctraj[k] = v

# Convert to DF for plotting and sampling.
df = to_dataframe(ctraj, nframes, dt=1)

# Plot trajectories
sample = df.sample(frac=0.1, axis=0)
sample.sort_values(by=['Prod_ID', 'Site_ID', 'Time_ps'], inplace=True)
g = sns.FacetGrid(sample, col='Prod_ID',hue='Site_ID', col_wrap=10)
g.map(plt.scatter, 'Time_ps', 'Trajectory', alpha=0.5)
g.set(ylim=(-0.5,num_clusters))
g.fig.tight_layout()
plt.savefig('figures/{}_cluster_trajectory.pdf'.format(feature))

# Plot  histograms
g = sns.FacetGrid(sample, col='Prod_ID',hue='Site_ID', col_wrap=10)
g = g.map(plt.hist, 'Trajectory', bins=range(num_clusters), histtype='step', lw='5')
g.fig.tight_layout()
plt.savefig('figures/{}_cluster_hist.pdf'.format(feature))


# Save dataframe
save_generic(df, 'clusters/{}_cluster_trajectory.pickl'.format(feature))

# Sampling (only plot 10 random assortments
to_plot = np.random.choice(range(num_clusters), min(10, num_clusters), replace=False)
for i in to_plot:
    num_samples = 100
    df_smp = df.ix[df['Trajectory']==i, ['Key', 'Time_ps']].sample(num_samples)
    inds = zip(df_smp['Key'], df_smp['Time_ps'])

    # Use loc because sample_dimension is nice
    traj = md.join(
        md.load_frame(meta.loc[traj_i]['traj_fn'], index=frame_i, top=meta.loc[traj_i]['top_fn'])
        for traj_i, frame_i in inds
    )

    # Original trajectories include both BT1 and BT2 so need to superpose
    ref = md.load('topology.pdb')
    traj.superpose(reference=ref)

    # Save
    traj_fn = "clusters/{0}_cluster-{1}.dcd".format(feature, i)
    backup(traj_fn)
    traj.save(traj_fn)