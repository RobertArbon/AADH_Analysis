from msmbuilder.io import load_generic
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.pylab import plt
import sys
import seaborn as sns
colors = sns.color_palette("colorblind", 8)
import pandas as pd


df = load_generic('grid_search.pickl')
df = df.filter(regex=("param_.*|split.*"))
id_cols = list(df.filter(regex=("param_.*")).columns)
var_cols = list(df.filter(regex=("split.*")).columns)

df = pd.melt(df,id_vars=id_cols, value_vars=var_cols)


# g = sns.factorplot(x="param_tica__lag_time", y="mean_test_score", hue="param_tica__n_components",
#                    col='param_cluster__n_clusters', data=df)
# parameter = 'n_clusters'
# fig, ax = plt.subplots()
# ax.errorbar(x=df['param_cluster__{}'.format(parameter)], y=df['mean_test_score'], yerr=df['std_test_score']*2)
# ax.set_xlabel('{} value'.format(parameter))
# ax.set_ylabel('GMRQ Score')
# plt.savefig('results.pdf')