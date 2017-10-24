import mdtraj as md
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


hod2 = 1.9

froot = 'pcca_cluster-{}.dcd'
all_dfs = []
for j in [2,3]:
    for i in [0,1,-1]:

        fname = froot.format(i)
        ref = md.load('../MSM_Reactants_Only/2agy_rxt.pdb')
        traj = md.load(fname, top='../MSM_Reactants_Only/2agy_rxt.psf')

        p1 = traj.topology.select("(resname TTW and name HI{}) or (resname ASP and name OD1)".format(j))
        p2 = traj.topology.select("(resname TTW and name HI{}) or (resname ASP and name OD2)".format(j))

        pairs = np.array([p1, p2])
        distances = md.compute_distances(traj, pairs)*10 # Nanometers to Angstroms
        df = pd.DataFrame(data=distances, columns=['OD1','OD2'])
        df['PCCA_state'] = i
        df['H-donor'] = 'HI{}'.format(j)
        df = pd.melt(frame=df, id_vars=['PCCA_state', 'H-donor'], var_name='H-acceptor', value_name='Distance')
        all_dfs.append(df)


df = pd.concat(all_dfs)
print(df.head())

with sns.plotting_context('talk', font_scale=1):
    sns.set_style('white')
    cols = sns.color_palette('colorblind', 3)
    grid = sns.FacetGrid(df, col='H-donor', size=5, legend_out=True, )
    grid = grid.map(sns.violinplot, 'PCCA_state', 'Distance', 'H-acceptor', palette=cols[:2])
    for ax in grid.axes.flatten():
        ax.hlines(hod2, xmin=-1, xmax=2.5, colors=cols[-1], label='QM/MM (HI2-OD1/2)')
        ax.legend(loc=2)

    grid.set_ylabels(r'O-H Distance ($\AA$)')
    grid.set_xlabels('PCCA state')

    plt.savefig('OH_distance.png', transparent=True)
