
from msmbuilder.io import load_meta, preload_tops, save_generic, itertrajs, backup
import mdtraj as md
from msmbuilder.cluster import LandmarkAgglomerative
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
from utilities import to_dataframe

df = pd.read_pickle('clusters/rmsd_cluster_trajectory.pickl')
num_clusters = np.max((df['Trajectory'].unique())) + 1

def plot_rmsd_trajectory():
    cluster_labels = np.arange(1,num_clusters,2)

    # Plot trajectories
    print('Charting trajectories')
    sample = df.sample(frac=0.1, axis=0)
    sample.sort_values(by=['Prod_ID', 'Site_ID', 'Time_ps'], inplace=True)

    with sns.plotting_context("notebook", font_scale=2):

        g = sns.FacetGrid(sample, col='Prod_ID',hue='Site_ID', col_wrap=10)
        g.map(plt.scatter, 'Time_ps', 'Trajectory', alpha=0.5)
        g.set(ylim=(-0.5,num_clusters), yticks=cluster_labels)
        g.set_ylabels("Cluster")
        g.set_xlabels("Time $ps$")
        g.set_titles("")
        g.fig.subplots_adjust(wspace=0.05, hspace=0.05)
        plt.savefig('figures/rmsd_cluster_trajectory.png', transparent=True)



# Sampling
def sample_clusters():

    meta = load_meta()
    tops = preload_tops(meta)
    print('Sampling trajectories')
    ref = md.load('topology.pdb')
    for i in range(int(num_clusters)):
        print(i)
        df_smp = df.ix[df['Trajectory']==i, ['Key', 'Time_ps']].sample(100)
        inds = zip(df_smp['Key'], df_smp['Time_ps'])

        # Use loc because sample_dimension is nice
        traj = md.join(
            md.load_frame(meta.loc[traj_i]['traj_fn'], index=frame_i, top=meta.loc[traj_i]['top_fn'])
            for traj_i, frame_i in inds
        )

        # Original trajectories include both BT1 and BT2 so need to superpose
        traj.superpose(reference=ref)

        # Save
        traj_fn = "clusters/rmsd_cluster-{}.dcd".format(i)
        backup(traj_fn)
        traj.save(traj_fn)


sample_clusters()