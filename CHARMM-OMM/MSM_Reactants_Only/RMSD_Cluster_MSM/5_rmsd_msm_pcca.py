# The pipeline for constructing an MSM
from msmbuilder.io import load_generic, load_trajs, save_trajs, save_generic
from msmbuilder.msm import MarkovStateModel
from msmbuilder.lumping import PCCAPlus
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
    meta, ctraj_dict = load_trajs('ctraj-200')
    long_ctrajs = [np.squeeze(traj) for traj in ctraj_dict.values() if traj.shape[0] > 1000]

    ps_to_ns = 1000
    n_ts = 10
    timescales = []
    lags = []
    for msm in all_msms:
        timescales.append(msm.timescales_[:n_ts])
        lags.append(msm.get_params()['lag_time'])
    lags = np.array(lags)
    timescales = np.array(timescales).T/ps_to_ns
    msm = all_msms[np.extract(lags == 2000,np.arange(len(lags)))[0]]

    m = 2
    # 1 timescales means 2 states
    pcca = PCCAPlus.from_msm(msm, n_macrostates=m)
    pcca_mapping = pcca.x
    print(len(pcca_mapping))
    plt.clf()
    sns.color_palette('colorblind', m)
    with sns.plotting_context("notebook", font_scale=1.5):
        vec = msm.right_eigenvectors_
        n_states = vec.shape[0] # may be less than 200 as T may be non-ergodic.
        fig, axes = plt.subplots(nrows=m, sharex=True)
        for i in range(m):
            for j in range(m):
                mask = pcca_mapping==j
                axes[i].bar(np.arange(n_states)[mask], vec[mask,i],label='PCCA State {}'.format(j), align='center')
            axes[i].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            axes[i].legend()
            axes[i].set_ylabel('Cluster projection')

        plt.xlabel('Cluster')
        plt.savefig('figures/rmsd_msm_right_eigenvectors-pcca.png', transparent=True)

    plt.clf()
    sns.color_palette('colorblind', m)
    with sns.plotting_context("notebook", font_scale=1.5):
        vec = msm.left_eigenvectors_
        n_states = vec.shape[0] # may be less than 200 as T may be non-ergodic.
        fig, axes = plt.subplots(nrows=m, sharex=True)
        for i in range(m):
            for j in range(m):
                mask = pcca_mapping==j
                axes[i].bar(np.arange(n_states)[mask], vec[mask,i],label='PCCA State {}'.format(j), align='center')
            axes[i].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            axes[i].legend()
            axes[i].set_ylabel('Cluster projection')

        plt.xlabel('Cluster')
        plt.savefig('figures/rmsd_msm_left_eigenvectors-pcca.png', transparent=True)

    # Transforms:
    msm_traj = {}
    pcca_traj = {}
    for k, v in ctraj_dict.items():
        print(k)
        msm_traj[k] = msm.partial_transform(np.squeeze(v),mode='fill')
        pcca_traj[k] = pcca.partial_transform(np.squeeze(v), mode='fill')

    save_trajs(msm_traj, 'msm-traj-200', meta)
    save_generic(msm, 'msm-200.pickl')
    save_trajs(pcca_traj, 'pcca-2-traj', meta)
    save_generic(pcca, 'pcca-2.pickl')