#
# SCALES DATA
#
# NB: Can't do this in parallel (easily) as the algorithm requires all trajectory data to do transform
#

from msmbuilder.preprocessing import RobustScaler
import numpy as np
from msmbuilder.io import load_trajs, save_trajs, save_generic
import matplotlib
matplotlib.use('Agg')
from matplotlib.pylab import plt
from utilities import plot_box

if __name__ == '__main__':

    # Load
    feature_name = 'Positions'
    meta, feature_trajs = load_trajs('Unscaled-{}-ftraj'.format(feature_name))

    # Select scaler
    featurizer = RobustScaler()

    # Transform values
    featurizer.fit_transform(feature_trajs.values())
    scaled_trajs = {}
    for k, v in feature_trajs.items():
        scaled_trajs[k] = featurizer.partial_transform(v)

    # Plot unscaled features
    ftrajs = np.concatenate([fx[::100] for fx in scaled_trajs.values()])
    fig, ax = plt.subplots(figsize=(15, 5))
    plot_box(ax, fxx=ftrajs, feature_name='Scaled {}'.format(feature_name))
    fig.tight_layout()
    fig.savefig("Scaled-{}-box.pdf".format(feature_name))

    # Save
    save_trajs(scaled_trajs, 'Scaled-{}-ftraj'.format(feature_name), meta)
    save_generic(featurizer, 'Scaled-{}-featurizer.pickl'.format(feature_name))

