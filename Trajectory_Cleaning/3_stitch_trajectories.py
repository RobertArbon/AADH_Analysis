import mdtraj as md
from os.path import join
from subprocess import run, PIPE
import os
import datetime

# dir = 'active_site'
# pdb = 'act_site.pdb'
# as_num = 2
# label = 'as{}'.format(as_num)
# pdb_path = join(dir, pdb)
# file_dir = join(dir, 'unaligned')
#
# for i in range(1,11):
#     u = md.load(join(file_dir,'2agy-310k-1atm-prod{0}-{1}.dcd'.format(i, label)), top=pdb_path)
#     for j in range(1,11):
#         v = md.load(join(file_dir, '2agy-310k-1atm-prod{0}.{1}-{2}.dcd'.format(i,j, label)), top=pdb_path)
#         w = u.join(v, discard_overlapping_frames=True)
#         out_file = '2agy-310k-1atm-prod{0}.{1}-{2}-2ns.dcd'.format(i,j, label)
#         print 'Writing {}'.format(out_file)
#         w.save(join(file_dir, out_file))


if __name__ == '__main__':

    pdb = '/Users/robert_arbon/Code/AADH/MD/common/2agy_final.pdb'
    input_dir = '/Volumes/REA_data/AADH/traj_2_no_water'
    output_dir = '/Volumes/REA_data/AADH/traj_3_combined'

    # If more trajectories from second round are needed then
    # then initial 1-ns of trajectory will possibly need to be removed due to same starting point.
    os.chdir(input_dir)

    start_r1, end_r1 = 1, 10    # The first round of trajectory sampling
    start_r2, end_r2 = 1, 1     # second round (new random seeds)
    start_r3, end_r3 = 1, 5     # indices to combine
    with open(os.path.join(output_dir,'stich.log'), 'w') as fh:
        for i in range(start_r1, end_r1+1):
            for j in range(start_r2, start_r2+1):
                inp_string = []
                for k in range(start_r3, end_r3+1):  # These will concatenated
                    inp_string.append('2agy-310k-1atm-prod{0}.{1}-{2}.nowat.dcd'.format(i, j, k))

                out_string = os.path.join(output_dir, '2agy-310k-1atm-prod{0}.{1}-nowat-comb.dcd'.format(i,j))
                output = run(['catdcd', '-o', out_string] + inp_string, stdout=PIPE, universal_newlines=True)
                fh.write(output.stdout)
                if output.stderr is not None:
                    fh.write(output.stderr)
                fh.write('-'*100+'\n')





