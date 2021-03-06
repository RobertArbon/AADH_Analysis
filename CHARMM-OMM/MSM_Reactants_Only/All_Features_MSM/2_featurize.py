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

    # FEATURIZER PARAMETERS
    reference = md.load('topology.pdb')
    atom_pairs = np.array([[x, y] for x, y in
                           combinations(range(reference.topology.n_atoms), 2)])
    idx_dir = '/Users/robert_arbon/Code/AADH/Analysis/msm/active_site'
    bonds = pickle.load(open(join(idx_dir, 'act_site_bond.p'), 'rb')) - 1
    angles = pickle.load(open(join(idx_dir,'act_site_ang.p'), 'rb')) - 1
    dihedrals = pickle.load(open(join(idx_dir,'act_site_dihe.p'), 'rb')) - 1


    # MSMBuilder FEATURIZERS
    RawPos = RawPositionsFeaturizer(ref_traj=reference)
    Pairs = AtomPairsFeaturizer(pair_indices=atom_pairs)
    DRID = DRIDFeaturizer()
    contacts = ContactFeaturizer(contacts=np.array([[0,1]]))
    bonds = AtomPairsFeaturizer(pair_indices=bonds)
    featurizers = [('contacts', contacts),
                   ('raw_positions', RawPos),
                   ('atom_pairs', Pairs),
                   ('drid', DRID),
                   ('bonds', bonds)]

    for name, feat in featurizers:
        print('Featurizing {}'.format(name))
        args = zip(meta.iterrows(), [feat]*meta.shape[0], [tops]*meta.shape[0])

        with Pool() as pool:
            feature_trajs = dict(pool.imap_unordered(msmb_feat, args))

        # # Create save objects
        # featurizer = dict([(x[0], x[2]) for x in feature_trajs])
        # feature_trajs = dict([(x[0], x[1]) for x in feature_trajs])

        selector = VarianceThreshold()
        selector.fit([traj for traj in feature_trajs.values()])
        ftrajs = {}
        for k, v in feature_trajs.items():
            ftrajs[k] = np.squeeze(selector.transform([v]))

        # SAVE
        save_trajs(ftrajs, 'featurized_trajectories/{}-ftraj'.format(name), meta)
        save_generic(feat, 'featurized_trajectories/{}-featurizer.pickl'.format(name))


    # pyEMMA FEATURIZERS
    featurizers = [('angles', 'add_angles', angles),
                   ('dihedrals', 'add_dihedrals', dihedrals)]

    for name, feat, indices in featurizers:
        print('Featurizing {}'.format(name))
        args = zip(meta.iterrows(), [feat] * meta.shape[0], [tops] * meta.shape[0],
                   [indices]*meta.shape[0])

        # Fit features
        with Pool() as pool:
            feature_trajs = dict(pool.imap_unordered(pyemma_feat, args))

        # # Create save objects
        # feat_desc = dict([(x[0], x[2].describe()) for x in feature_trajs])
        # feature_trajs = dict([(x[0], x[1]) for x in feature_trajs])

        # Remove zero variance features
        selector = VarianceThreshold()
        selector.fit([traj for traj in feature_trajs.values()])
        ftrajs = {}
        for k, v in feature_trajs.items():
            ftrajs[k] = np.squeeze(selector.transform([v]))
            # mask = selector.get_support()
            # feat_desc[k] = feat_desc[k][mask]


        save_trajs(ftrajs, 'featurized_trajectories/{}-ftraj'.format(name), meta)
        # save_generic(feat_desc, 'featurized_trajectories/{}-description.pickl'.format(name))



