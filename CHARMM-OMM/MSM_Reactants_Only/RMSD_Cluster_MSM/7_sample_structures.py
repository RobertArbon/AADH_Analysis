# The pipeline for constructing an MSM
from msmbuilder.io import load_trajs, backup
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
from matplotlib.pylab import plt
import pandas as pd
import math
import mdtraj as md
import sys

sns.set_style("white")


def to_dataframe(traj_dict):
    all_dfs = []
    for k, v in traj_dict.items():
        ids = np.array(k*len(v)).reshape(len(v), len(k))
        v = v[:,np.newaxis]
        data = np.concatenate((ids, v, np.arange(len(v)).reshape(-1,1)), axis=1)
        id_cols = ['PDB', 'Temp', 'Pres', 'Prod_ID', 'Site_ID', 'Trajectory', 'Frame']
        df = pd.DataFrame(data=data, columns=id_cols)
        df['Temp'] = df['Temp'].astype(float)
        df['Pres'] = df['Pres'].astype(float)
        df['Site_ID'] = df['Site_ID'].astype(int)

        df['Key'] = df[['PDB', 'Temp', 'Pres', 'Prod_ID', 'Site_ID']].apply(tuple, axis=1)
        all_dfs.append(df)

    df = pd.concat(all_dfs)
    df.ix[df['Trajectory'].isnull(), 'Trajectory'] = -1
    df.ix[df['Trajectory']=='nan', 'Trajectory'] = -1
    df['Trajectory'] = df['Trajectory'].astype(float).astype(int)
    df['Frame'] = df['Frame'].astype(int)
    return df


def plot_cluster_trajectory(df):
    cluster_labels = df['Trajectory'].unique()
    ymin = np.min(cluster_labels)-0.5
    ymax = np.max(cluster_labels)+0.5

    # Plot trajectories
    print('Charting trajectories')
    # sample = df.iloc[0:-1:2,:]
    long_trajs = ['{}.1'.format(x) for x in range(1,11)]
    sample = df.ix[df['Prod_ID'].isin(long_trajs),:]

    sample.sort_values(by=['Prod_ID', 'Site_ID', 'Frame'], inplace=True)
    print(sample.head())
    with sns.plotting_context("notebook", font_scale=2):

        g = sns.FacetGrid(sample, col='Prod_ID',hue='Site_ID', col_wrap=5)
        g.map(plt.scatter, 'Frame', 'Trajectory', alpha=0.5, marker='o', s=5)
        g.set(ylim=(ymin,ymax), yticks=cluster_labels)
        g.set_ylabels("Cluster")
        g.set_xlabels("Time $ps$")
        g.set_titles("")
        g.fig.subplots_adjust(wspace=0.05, hspace=0.05)
        plt.savefig('figures/pcca_cluster_trajectory.png', transparent=True)


def sample_clusters(meta, trajs, df):

    clust_id = df['Trajectory'].unique()

    for i in clust_id:
        print(i)

        df_smp = df.ix[df['Trajectory']==i, ['Key', 'Frame']].sample(1000)
        inds = zip(df_smp['Key'], df_smp['Frame'])

        # Use loc because sample_dimension is nice
        traj = md.join(
            md.load_frame(meta.loc[traj_i]['traj_fn'], index=int(frame_i), top=meta.loc[traj_i]['top_fn'])
            for traj_i, frame_i in inds
        )

        # Save
        traj_fn = "/Users/robert_arbon/Code/AADH/Analysis/KR_Comparison/pcca_cluster-{}.dcd".format(i)
        backup(traj_fn)
        traj.save(traj_fn)


if __name__ == "__main__":
    meta, ctraj = load_trajs('pcca-2-traj')

    df = to_dataframe(ctraj)
    # sample_clusters(meta, ctraj, df)

    plot_cluster_trajectory(df)