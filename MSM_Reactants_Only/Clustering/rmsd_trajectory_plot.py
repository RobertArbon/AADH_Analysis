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


meta = load_meta()
tops = preload_tops(meta)

ref = md.load('topology.pdb')
feat = RMSDFeaturizer(reference_traj=ref)

args = zip(meta.iterrows(), [feat] * meta.shape[0], [tops] * meta.shape[0])

with Pool() as pool:
    ftrajs = dict(pool.imap_unordered(msmb_feat, args))

# Squeeze and extend short trajectories with zeros
# MSMBuilder does rmsd in nm so multiply by 10 to get angstroms
nframes = int(np.max(meta['nframes'].unique()[0]))
ns_to_ang = 10
rtrajs = {}
for k, v in ftrajs.items():
    v = np.squeeze(v)
    diff = nframes-v.shape[0]
    v = np.append(v, np.zeros(diff)+np.nan)*ns_to_ang
    rtrajs[k] = v

# Create DataFrame
id_cols = ['PDB', 'Temp', 'Pres','Prod_ID', 'Site_ID']
dt = float(meta['step_ps'].unique()[0])
df = pd.DataFrame.from_records(data=rtrajs)
df['Time_ps'] = np.arange(nframes)*dt
df = pd.melt(df, id_vars=['Time_ps'], var_name='Production_ID', value_name='RMSD')
df[id_cols] = pd.DataFrame(df['Production_ID'].tolist())
del df['Production_ID']

df = df.join(df.groupby(id_cols)['RMSD'].rolling(10).mean()
             .reset_index(level=[0,1,2,3,4]), rsuffix='_r')
df.drop(labels=[x+'_r' for x in id_cols], axis=1, inplace=True)

long_trajs = ['{}.1'.format(x) for x in range(1,11)]
# Plot rolling
with sns.plotting_context("notebook", font_scale=2):
    sample = df.ix[df['Prod_ID'].isin(long_trajs),:].sample(frac=0.1, axis=0)
    sample.sort_values(by=['Prod_ID', 'Site_ID', 'Time_ps'], inplace=True)
    g = sns.FacetGrid(sample, col='Prod_ID',hue='Site_ID', col_wrap=5)
    g.map(plt.plot, 'Time_ps', 'RMSD_r')
    g.set_ylabels("RMSD $\AA$")
    g.set_xlabels("")
    g.set_titles("")
    g.fig.subplots_adjust(wspace=0.05, hspace=0.05)
    plt.savefig('rmsd_trajectory_long.png', transparent=True)

# Save dataframe
save_generic(df, 'rmsd_trajectory.pickl')
