import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
from glob import glob
import pickle
import numpy as np
import matplotlib.pyplot as plt
import sys

# GLOBALS
indir = '../data'
ttw = np.arange(0,45)
asp = np.arange(45,57)
groups = [ttw, asp]

# LOAD TOPOLOGY

# LOAD FEATURES
# these must be zero indexed.  MDAnalysis and PDBs have 1-indexing.
angles = pickle.load(open('act_site_ang.p', 'rb'))-1
dihedrals = pickle.load(open('act_site_dihe.p', 'rb')) - 1
bonds = pickle.load(open('act_site_bond.p', 'rb')) - 1

# TOPOLOGY
topfile = 'act_site.pdb'
feat = coor.featurizer(topfile)

# Had to comment out the CRYST record from PDB.
feat.add_angles(angles, deg=True, cossin=True, periodic=False)
traj_list = glob('/aligned/*.dcd')
inp = coor.source(traj_list, feat)

# print(inp.get_output()[0])
description = inp.describe()

# fixed: 'ANGLE: COS(TTW 39 N 0 - TTW 39 CA 2 - TTW 39 HA 3)', 'ANGLE: SIN(TTW 39 N 0 - TTW 39 CA 2 - TTW 39 HA 3)',
# not-fixed: 'ANGLE: COS(TTW 39 N 0 - TTW 39 CA 2 - TTW 39 HA 3)', 'ANGLE: SIN(TTW 39 N 0 - TTW 39 CA 2 - TTW 39 HA 3)'
# cossin=False: 'ANGLE: TTW 39 N 0 - TTW 39 CA 2 - TTW 39 HA 3 '
sys.exit()

# feat.add_dihedrals(dihedrals, deg=True, cossin=True, periodic=False)
# feat.add_distances(bonds, periodic=False)
# feat.add_group_mindist(groups, periodic=False)

# LOAD TRAJECTORIES


lags = np.arange(900,1000,50)
dim = 15
times = np.zeros((len(lags), dim))
for idx, lag in enumerate(lags):
    print lag
    tica_obj = coor.tica(inp, lag=lag, var_cutoff=0.9, kinetic_map=True)
    times[idx] = -lag/np.log(tica_obj.eigenvalues[:dim])

plt.plot(lags, times)
plt.show()