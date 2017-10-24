from msmbuilder.io import save_trajs, GenericParser, save_meta, gather_metadata, preload_tops, save_generic
from msmbuilder.featurizer import RawPositionsFeaturizer
from msmbuilder.preprocessing import RobustScaler
from msmbuilder.decomposition import tICA
import os
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt


def identity(x):
    return x

# TODO pass a dictionary of parameters
def load_metadata(traj_dir, top):
    """
    Loads metadata of features and saves them.

    :param traj_dir: directory containing trajectories
    :param top: topology file name
    :return: metadata data frame
    """
    re_pattern = '(\w+)-([0-9]{3})k-([0-9])atm-prod([0-9]+\.[0-9]+).*BT([0-9]+)*'
    captured_group_names = ['PDB', 'Temp', 'Pressure', 'Prod_Round', 'Act_Site']
    captured_group_transforms = [identity, float, float, identity, int]
    time_step = 1 # in picoseconds
    file_type = 'dcd'

    parser = GenericParser(re_pattern,
                           group_names=captured_group_names,
                           group_transforms=captured_group_transforms,
                           top_fn=top, step_ps=time_step)
    meta = gather_metadata(os.path.join(traj_dir, "*.{}".format(file_type)), parser)
    save_meta(meta)
    return meta


def featurize(featurizer, meta_data):

    tops = preload_tops(meta)

    def feat(irow):
        i, row = irow
        traj = md.load(row['traj_fn'], top=tops[row['top_fn']])
        feat_traj = featurizer.partial_transform(traj)
        return i, feat_traj

    feature_trajs = dict(map(feat, meta.iterrows()))

    save_trajs(feature_trajs, 'ftrajs', meta)
    save_generic(featurizer, 'featurizer.pickl')

    return feature_trajs


# TODO make this accept meta-data
def combine(trajs):
    try:
        all_trajs = list(trajs.values())
        all_trajs = np.concatenate(all_trajs, axis=0)
        return all_trajs
    except AttributeError:
        print("Trajs must be dictionary")
        raise


def plot_features(feat, name='Features.png', feature_name='Feature', ordered=False):
    """
    plots val vs feat_index
    :param feat: 1 x num_features array o
    :return: None
    """
    fig, ax = plt.subplots()
    ind = np.arange(len(feat))
    labels = ind
    xlab = 'Feature index'
    if ordered:
        labels = np.argsort(feat)
        feat = np.sort(feat)
        xlab += ' (sorted)'
    width = 1
    stride = 10
    ax.bar(ind, feat, width=width, align='center')
    plt.tight_layout()
    ax.set_xlabel(xlab)
    ax.set_ylabel('{0} value'.format(feature_name))
    ax.set_xticks(ind[::stride])
    ax.set_xticklabels(labels[::stride])
    plt.savefig(name)



if __name__ == "__main__":

    trajectory_dir = '/Volumes/REA_Data/AADH/traj_5_rxts'
    topology_file = '/Users/robert_arbon/Code/AADH/Analysis/MSM_Reactants_Only/2agy_rxt.psf'
    reference_file = '/Users/robert_arbon/Code/AADH/Analysis/MSM_Reactants_Only/2agy_rxt.pdb'
    reference_traj = md.load(reference_file)

    # Load the meta data
    meta = load_metadata(traj_dir=trajectory_dir, top=topology_file)

    # Featurize
    feature = RawPositionsFeaturizer(ref_traj=reference_traj)
    ftrajs = featurize(featurizer=feature, meta_data=meta)

    # Summarize
    variance = np.var(combine(ftrajs), axis=0)
    plot_features(variance, name='Variance.png', feature_name='Variance', ordered=False)

    # Normalize
    scaler = RobustScaler()
    strajs = scaler.fit_transform(ftrajs)

    # perform tICA
    tica_obj = tICA(n_components=10, lag_time=10, kinetic_mapping=True)
    tica_traj = tica_obj.fit_transform(strajs)









