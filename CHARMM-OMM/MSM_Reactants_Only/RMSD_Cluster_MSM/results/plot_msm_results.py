from msmbuilder.io import load_generic
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.pylab import plt
import sys
import seaborn as sns
colors = sns.color_palette("colorblind", 8)
import pandas as pd



def plot_lv(mod):
    lvs = mod.left_eigenvectors_
    nstates = min(5,lvs.shape[1])
    fig, axes = plt.subplots(nrows=nstates, sharey=False, sharex=True)
    for idx, ax in enumerate(axes):
        ax.bar(range(lvs.shape[0]), lvs[:,idx])
    plt.savefig('msm-lvs.png')

def plot_rv(mod):
    lvs = mod.right_eigenvectors_
    nstates = min(5,lvs.shape[1])
    fig, axes = plt.subplots(nrows=nstates, sharey=False, sharex=True)
    for idx, ax in enumerate(axes):
        ax.bar(range(lvs.shape[0]), lvs[:,idx])
    plt.savefig('msm-rvs.png')

msm = load_generic('msm-lag-4000-nclusters-20.pickl')

plot_lv(msm)
plot_rv(msm)
print(msm.timescales_)
print(msm.uncertainty_timescales())
plt.matshow(msm.transmat_)
plt.savefig('msm-transmat.png')


