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
    meta, feature_trajs = load_trajs('ftraj')

    # Select scaler
    featurizer = RobustScaler()

    # Transform values
    featurizer.fit_transform(feature_trajs.values())
    scaled_trajs = {}
    for k, v in feature_trajs.items():
        scaled_trajs[k] = featurizer.partial_transform(v)

    # Save
    sample = np.concatenate([fx for fx in scaled_trajs.values()])
    sample = sample[np.random.choice(sample.shape[0], 1000, replace=False), :]
    variance = np.apply_along_axis(np.var, axis=0, arr=sample)
    order = np.argsort(variance)
    ord_var = variance[order]
    labels = [str(x) for x in ord_var[::10]]
    ind = range(variance.shape[0])
    fig, ax = plt.subplots()
    ax.plot(ind, ord_var)
    plt.savefig('ScaledFeatureVariance.png')


    save_trajs(scaled_trajs, 'straj', meta)
    save_generic(featurizer, 'scaler.pickl')

