# import MDAnalysis
# from os.path import join, isfile
# from glob import glob
# from copy import copy
#

import os
import shutil
import re
from multiprocessing import Pool, cpu_count


# def write_backbone((dirname, sele, traj)):
#     psf = '../common/2agy_final.psf'
#     bname = traj.split('/')[-1].split('.dcd')[0]
#     outname = join(dirname, bname+'-'+dirname+'.dcd')
#     if not isfile(outname):
#         print 'Writing ', bname
#         u = MDAnalysis.Universe(psf, traj)
#         backbone = u.select_atoms('protein and backbone')
#         natoms = backbone.n_atoms
#         with MDAnalysis.Writer(outname, multiframe=True, n_atoms=backbone.n_atoms) as W:
#             for ts in u.trajectory:
#                 W.write(backbone)
#
# if __name__ == '__main__':
#     traj_list = glob('/Volumes/REA_data/AADH/MD/Trajectories/*.dcd')
#     p = Pool(cpu_count())
#     p.map(write_backbone, traj_list)

# This script imports data from the RDSF storage and places it on local storage for entry into analysis pipeline.
# Turns out it's not so easy to use shutils over a network.  Can't be bothered to figure it out so I'll do a
# Copy-paste of the RDSF files
# 1. Finds all the DCDs in rounds 2 and onwards and imports them into the same directory
# 2. Renames them so that they all follow a sensible naming convention.
# Explanation:
#   prod1 : 10 x initial runs
#   prod1.1 : Further 10 per initial run.  All started from different random seeds.
#   prod1.1.1 : Continuation of prod1.1. with same random seeds.
# Renaming:
#   prod1.1 -> prod1.1-1 i.e. first ns
#   prod1.1.2 -> prod1.1-2 i.e second ns.


def move_dirs(args):
    """
    Moves the data from RDSF to local drive
    :param src: source directory
    :return: None
    """
    src = args[0]
    dst = args[1]
    print("Moving from: {}".format(src))
    print("         to: {}".format(dst))
    shutil.move(src, dst)
    return


def find_dcds(src):
    """
    gets all dcds in a root directory
    :param src: root directory in which dcd's can be found
    :return:
    """

    dcd_paths = []

    for root, dirs, files in os.walk(src):
        for filename in files:
            if filename.endswith(".dcd"):
                dcd_paths.append(os.path.join(root, filename))

    return dcd_paths


def rename_file(file_path, pattern, replacement):
    """
    Renames the DCD files in the round-2 by inserting something before file extension
    :param file_path: the old path of the file
    :param pattern: the ending that should be replaced
    :param replacement: the new ending
    :return:
    """
    old_file_name = os.path.basename(file_path)
    new_file_name = re.sub(pattern, replacement, old_file_name)
    return new_file_name


if __name__ == "__main__":

    file_root = '2agy-310k-1atm-prod'
    dir_root = 'prod-'
    dst_root = '/Volumes/REA_Data/AADH/trajectories'

    src_files = find_dcds('/Volumes/REA_Data/AADH/round_2')
    dst_files = [os.path.join(dst_root, rename_file(x, ".dcd", "-1.dcd")) for x in src_files]
    args = zip(src_files, dst_files)

    p = Pool(cpu_count()-1)
    p.map(move_dirs, args)

    src_files = find_dcds('/Volumes/REA_Data/AADH/round_3')
    dst_files = [os.path.join(dst_root, rename_file(x, r'(?<=\.[\d])\.(?=[\d]{1,2}.dcd)', '-')) for x in src_files]
    args = zip(src_files, dst_files)

    p = Pool(cpu_count()-1)
    p.map(move_dirs, args)




