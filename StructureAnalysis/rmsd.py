import mdtraj as md
import pandas as pd
import numpy as np


def chart_rmsd(traj_list, row_max):
    """
    :param rmsd_trajs: list of np.arrays corresponding to rmsd of trajectories
    :param chart_name: save to chart_name
    :return: None
    """
    if len(traj_list) > 0:
        fig = plt.figure()
        ntraj = len(traj_list)
        nrow = min(row_max,ntraj)
        ncol = ntraj//row_max
        fig, axes = plt.subplots(nrows=nrow, ncols=ncol, sharex=True, sharey=True)
        for row in range(nrow):
            if ncol > 1:
                for col in range(ncol):
                    idx = min(row + col*nrow, ntraj)
                    print(row, col, idx)
                    axes[row][col].plot(traj_list[idx])
            else:
                axes[row].plot(traj_list[row])
    else:
        raise ValueError("List of trajectories can't be empty")


def rmsd_per_residue_per_frame(target, reference, frame_number):
    """
    Computes the rmsd per residue
    :param target:
    :param reference:
    :param frame_number:
    :return: pandas dataframe of form: RESID RESNAME FRAME RMSD
    """
    num_frames = target.n_frames
    frames = np.arange(num_frames)
    all_dfs = []
    for residue in target.topology.residues:
        resname = np.repeat(residue.name, num_frames)
        resid = np.repeat(residue.index, num_frames)
        indices = [x.index for x in residue.atoms]
        rmsd = md.rmsd(target, reference, atom_indices=indices, frame=frame_number)
        data = {'RESID': resid, 'RESNAME': resname, 'FRAME': frames, 'RMSD': rmsd}
        all_dfs.append(pd.DataFrame(data=data))
    df = pd.concat(all_dfs)
    return df

def rmsd_per_residue(df):
    """
    returns average, mean, max of RMSD per residue
    :param df:
    :return: df of min, mean, max of RMSD per residue
    """
    grouped = df.groupby('RESID').aggregate([np.min, np.mean, np.max])['RMSD']
    return grouped

def plot_rmsd_per_residue(df, ax):
    """

    :param df:
    :return:
    """

