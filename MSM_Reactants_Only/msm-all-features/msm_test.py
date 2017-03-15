import numpy as np
from msmbuilder.msm import MarkovStateModel
from pyemma.msm import estimate_markov_model
import msmtools.generation as msmgen
import numpy.linalg as LA
import matplotlib
matplotlib.use('Agg')
from matplotlib.pylab import plt

R = 0.998
P = np.array([[R, 1-R],
              [1-R, R]])


eig = LA.eigvals(P)
t = -1/np.log(eig[1:])
fac = 10
ts = []
all_T = np.linspace(10, fac*t, num=5).astype(int)
for T in all_T:
# for i in range(10):
    X = np.zeros((T, 2))
    dtraj = msmgen.generate_traj(P, T)
    try:
        msm = estimate_markov_model([dtraj], lag=1)
        try:
            its = msm.timescales()[0]
        except IndexError:
            its = np.nan
    except RuntimeError:
        its = np.nan

    # ave = np.nanmean(np.array(temp))
    print(T, its)
    ts.append(its)

plt.scatter(all_T, ts, label='Single trajectory')



num_traj = 100
ts = []
all_T = np.linspace(10, fac*t, num=5).astype(int)
for T in all_T:
    trajs = []
    for i in range(num_traj):
        X = np.zeros((T, 2))
        dtraj = msmgen.generate_traj(P, T)
        trajs.append(dtraj)
    try:
        msm = estimate_markov_model(trajs, lag=1)
        try:
            its = msm.timescales()[0]
        except IndexError:
            its = np.nan
    except RuntimeError:
        its = np.nan

    print(T, its)
    ts.append(its)

plt.scatter(all_T, ts, label='{} trajectories'.format(num_traj))
plt.hlines(t, xmin=all_T[0], xmax=all_T[-1])
plt.savefig('figures/msm_traj_test.png')