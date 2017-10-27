from msmbuilder.io import load_trajs, save_generic, load_meta, preload_tops
from msmbuilder.feature_extraction import AtomPairsFeaturizer, FunctionFeaturizer
from msmbuilder.feature_selection import FeatureSelector, VarianceThreshold
import mdtraj as md
from multiprocessing import Pool
from parmed.amber import AmberParm
import numpy as np
import matplotlib.pyplot as plt



# Load data
meta = load_meta()
tops = preload_tops(meta)


def traj_load(irow):
    i, row = irow
    traj = md.load(row['traj_fn'], top=tops[row['top_fn']])
    return i, traj


def get_timings(meta):
    frames_tot = meta['nframes'].sum()
    n_frames = meta['nframes'].unique()
    assert (len(n_frames) == 1, 'Different trajectory lengths')
    n_frames = n_frames[0]
    dt = meta['step_ps'][0]
    to_ns = dt/1000
    t_max = n_frames*to_ns
    return to_ns, t_max, frames_tot


def get_ff_terms(top):
    parm = AmberParm(top)
    bonds = get_bonds(parm)
    angles = get_angles(parm)
    dihedrals = get_dihedrals(parm)
    return bonds, angles, dihedrals


def get_bonds(parm):
    bonds = np.array(list(set([(x.atom1.idx, x.atom2.idx) for x in parm.bonds])))
    return bonds


def get_angles(parm):
    angles = np.array(list(set([(x.atom1.idx, x.atom2.idx, x.atom3.idx) for x in parm.angles])))
    return angles


def get_dihedrals(parm):
    dihedrals = np.array(list(set([(x.atom1.idx, x.atom2.idx, x.atom3.idx, x.atom4.idx) for x in parm.dihedrals])))
    return dihedrals


def compute_angles(traj, indices, is_periodic=False, sincos=True):
    angles = md.compute_angles(traj, angle_indices=indices, periodic=is_periodic)
    if sincos:
        angles = np.concatenate((np.sin(angles), np.cos(angles)), axis=1)
    return angles

def compute_dihedrals(traj, indices, is_periodic=False, sincos=True):
    dihedrals = md.compute_dihedrals(traj, indices=indices, periodic=is_periodic)
    if sincos:
        dihedrals = np.concatenate((np.sin(dihedrals), np.cos(dihedrals)), axis=1)
    return dihedrals


if __name__ == '__main__':

    # Get data
    with Pool() as pool:
        trajs_dct = dict(pool.imap_unordered(traj_load, meta.iterrows()))
    trajs = [traj for traj in trajs_dct.values()]

    # Set up individual features
    bonds, angles, dihedrals = get_ff_terms(meta['top_abs_fn'][0])
    bonds_feat = AtomPairsFeaturizer(pair_indices=bonds, periodic=False)
    angles_feat = FunctionFeaturizer(compute_angles,
                                     func_args={'indices': angles, 'is_periodic': False, 'sincos':True})
    dihedrals_feat = FunctionFeaturizer(compute_dihedrals,  func_args={'indices': dihedrals, 'is_periodic': False,
                                                                       'sincos':True})

    # Build pipeline




