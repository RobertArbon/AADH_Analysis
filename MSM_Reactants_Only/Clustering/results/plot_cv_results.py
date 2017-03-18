from msmbuilder.io import load_generic
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.pylab import plt
import sys
import seaborn as sns
colors = sns.color_palette("colorblind", 8)
import pandas as pd
from glob import glob

files = glob("*.pickl")
all_dfs = []
for file in files:
    all_dfs.append(load_generic(file))

df = pd.concat(all_dfs)
df.sort_values(by='param_cluster__n_clusters', inplace=True)

# df = df.filter(regex=("param_.*|split.*"))
# id_cols = list(df.filter(regex=("param_.*")).columns)
# var_cols = list(df.filter(regex=("split.*")).columns)
#
# df = pd.melt(df,id_vars=id_cols, value_vars=var_cols, value_name='GMRQ')
# df['Data'] = df['variable'].str.extract('(test|train)', expand=True)

parameter = 'n_clusters'
fig, ax = plt.subplots()
ax.errorbar(x=df['param_cluster__{}'.format(parameter)], y=df['mean_test_score'], yerr=df['std_test_score']*2, label='Test score')
ax.errorbar(x=df['param_cluster__{}'.format(parameter)], y=df['mean_train_score'], yerr=df['std_train_score']*2, label='Train score')

ax.set_xlabel('{} value'.format(parameter))
ax.set_ylabel('GMRQ Score')
plt.legend()
plt.savefig('results.png', transparent=True)