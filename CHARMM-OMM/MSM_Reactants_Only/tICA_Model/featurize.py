#
# FEATURIZES DATA
#

from msmbuilder.featurizer import RawPositionsFeaturizer
import numpy as np
import mdtraj as md
from msmbuilder.io import load_meta, preload_tops, save_trajs, save_generic
from multiprocessing import Pool
import matplotlib
matplotlib.use('Agg')
from matplotlib.pylab import plt
from utilities import plot_box, feat


if __name__=='__main__':

    # Load
    meta = load_meta()
    tops = preload_tops(meta)

    # Select featurizer
    feature_name = 'Positions'
    reference = md.load('topology.pdb')
    featurizer = RawPositionsFeaturizer(ref_traj=reference)

    args = zip(meta.iterrows(), [featurizer]*meta.shape[0], [tops]*meta.shape[0])

    # Do it in parallel
    with Pool() as pool:
        feature_trajs = dict(pool.imap_unordered(feat, args))


    # Plot unscaled features
    ftrajs = np.concatenate([fx[::100] for fx in feature_trajs.values()])
    fig, ax = plt.subplots(figsize=(15, 5))
    plot_box(ax, fxx=ftrajs, feature_name='Unscaled {}'.format(feature_name))
    fig.tight_layout()
    fig.savefig("Unscaled-{}-box.pdf".format(feature_name))

    ## Save
    save_trajs(feature_trajs, 'Unscaled-{}-ftraj'.format(feature_name), meta)
    save_generic(featurizer, 'Unscaled-{}-featurizer.pickl'.format(feature_name))

