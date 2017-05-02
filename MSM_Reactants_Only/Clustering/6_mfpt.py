# The pipeline for constructing an MSM
from msmbuilder.io import load_trajs
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
from matplotlib.pylab import plt
import pandas as pd
import math
sns.set_style("white")


def plot(data, fname):
    plt.clf()
    width = 10
    data = np.array(data)
    bins = np.arange(0,max(data),step=width)
    sns.distplot(data, norm_hist=True, kde=True, bins=bins, label='No ones')
    plt.ylim((0,0.001))
    plt.xlabel('First passage time (ps)')
    plt.text(1000, 0.0005, 'MFPT = {0:4.2f} +/- {1:4.2f} ps'.format(data.mean(), 2*data.std()))
    plt.savefig(fname, transparanet=True)

if __name__ == "__main__":
    meta, ctraj = load_trajs('pcca-2-traj')
    fpt = {}
    fpt[(0, 1)] = []
    fpt[(1, 0)] = []
    for k, v in ctraj.items():
        count = 0
        for i in range(len(v)-1):
            v1,v2 = v[i:(i+2)]
            if math.isnan(v1) or math.isnan(v2):
                count = 0
            else:
                count += 1
                if v1 != v2:
                    fpt[(v1, v2)].append(count)
                    count = 0

    plot(fpt[(1,0)], fname='figures/1-0-mfpt.png')
    plot(fpt[(0,1)], fname='figures/0-1-mfpt.png')