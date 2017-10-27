# The pipeline for constructing an MSM
from msmbuilder.io import load_generic
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
from matplotlib.pylab import plt
from matplotlib.ticker import FormatStrFormatter

sns.set_style("white")

def print_timescales(timescales):
    pass


if __name__ == "__main__":
    all_msms = load_generic('rmsd_msms.pickl')

    ps_to_ns = 1000

    timescales = []
    lags = []
    n_ts = 10
    sns.set_palette("Blues_r", n_ts*2)
    for msm in all_msms:
        timescales.append(msm.timescales_[:n_ts])
        lags.append(msm.get_params()['lag_time'])
    lags = np.array(lags)
    timescales = np.array(timescales).T/ps_to_ns

    with sns.plotting_context("notebook", font_scale=1.5):
        for i in range(len(timescales)):

            plt.plot(lags, timescales[i])
        plt.xlabel('Lag time $ps$')
        plt.ylabel('Implied timescale $ns$')
        plt.savefig('figures/rmsd_msm_ts_vs_lag.png', transparent=True)
        ymax = (timescales[1].max()//10 + 1)*10
        plt.ylim((0,ymax))
        plt.savefig('figures/rmsd_msm_ts_vs_lag-detail.png', transparent=True)

    msm = all_msms[np.extract(lags == 2000,np.arange(len(lags)))[0]]

    m = 2
    plt.clf()
    with sns.plotting_context("talk", font_scale=1.5):
        vec = msm.left_eigenvectors_
        n_states = vec.shape[0] # may be less than 200 as T may be non-ergodic.
        fig, axes = plt.subplots(nrows=m, sharex=True)
        for i in range(m):
            axes[i].bar(range(n_states), vec[:,i],label='Eigenvector {}'.format(i+1), align='center', width=1)
            axes[i].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            axes[i].set_ylabel('Cluster Projection')
            axes[i].legend()

        plt.xlabel('Cluster')
        plt.savefig('figures/rmsd_msm_left_eigenvectors.png', transparent=True)

    plt.clf()
    with sns.plotting_context("talk", font_scale=1.5):
        vec = msm.right_eigenvectors_
        n_states = vec.shape[0] # may be less than 200 as T may be non-ergodic.
        fig, axes = plt.subplots(nrows=m, sharex=True)
        for i in range(m):
            axes[i].bar(range(n_states), vec[:,i],label='Eigenvector {}'.format(i+1), align='center', width=1)
            axes[i].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            axes[i].set_ylabel('Cluster projection'.format(i+1))
        plt.legend()
        plt.xlabel('Cluster')
        plt.savefig('figures/rmsd_msm_right_eigenvectors.png', transparent=True)




