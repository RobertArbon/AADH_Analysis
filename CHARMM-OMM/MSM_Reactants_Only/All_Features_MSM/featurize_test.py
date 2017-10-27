#
# FEATURIZES DATA
#

from msmbuilder.featurizer import *
from msmbuilder.feature_selection import *
import numpy as np
import mdtraj as md
from msmbuilder.io import load_meta, preload_tops, save_trajs, save_generic
from multiprocessing import Pool
import matplotlib
matplotlib.use('Agg')
from matplotlib.pylab import plt
from utilities import plot_box, msmb_feat, pyemma_feat
import sys
from itertools import combinations
import pickle
from os.path import join

if __name__ == '__main__':

    # Load
    meta = load_meta()
    tops = preload_tops(meta)
    dihed_feat = DihedralFeaturizer(types=['phi', 'psi', 'omega', 'chi1', 'chi2', 'chi3', 'chi4'])

    ## Featurize logic
    def feat(irow):
        i, row = irow
        traj = md.load(row['traj_fn'], top=tops[row['top_fn']])
        feat_traj = dihed_feat.partial_transform(traj)
        return i, feat_traj


    ## Do it in parallel
    with Pool() as pool:
        dihed_trajs = dict(pool.imap_unordered(feat, meta.iterrows()))

    irow = [(i, row) for i, row in meta.iterrows()]

    row = irow[0][1]
    traj = md.load(row['traj_fn'], top=tops[row['top_fn']])
    for d in dihed_feat.describe_features(traj):
        print(d['atominds'])

    # By looking at these indices and the structure, it's clear that 'these are not the dihedrals you're looking for'

    # ## Save
    # save_trajs(dihed_trajs, 'ftrajs', meta)
    # save_generic(dihed_feat, 'featurizer.pickl')

    FeatureSelector()