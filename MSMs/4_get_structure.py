from msmbuilder.io import save_generic
import parmed as pmd
from parmed.amber import AmberParm


struct = AmberParm('proc_traj/2agy-as1.prmtop')
angles = [(x.atom1.idx, x.atom2.idx, x.atom3.idx) for x in struct.angles]
bonds = [(x.atom1.idx, x.atom2.idx) for x in struct.bonds]
dihedrals = [(x.atom1.idx, x.atom2.idx, x.atom3.idx, x.atom4.idx) for x in struct.dihedrals]
print('Number of angles {}'.format(len(angles)))
print('Number of bonds {}'.format(len(bonds)))
print('Number of dihedrals {}'.format(len(dihedrals)))

# save_generic(bonds, 'proc_traj/2agy_as-1_bonds.pickl')
# save_generic(angles, 'proc_traj/2agy_as-1_angles.pickl')
# save_generic(dihedrals, 'proc_traj/2agy_as-1_dihedrals.pickl')

struct2 = AmberParm('proc_traj/2agy-as2.prmtop')
angles2 = [(x.atom1.idx, x.atom2.idx, x.atom3.idx) for x in struct2.angles]
bonds2 = [(x.atom1.idx, x.atom2.idx) for x in struct2.bonds]
dihedrals2 = [(x.atom1.idx, x.atom2.idx, x.atom3.idx, x.atom4.idx) for x in struct2.dihedrals]
tuples = [bonds, angles, dihedrals]
tuples2 = [bonds2, angles2, dihedrals2]

for k in range(len(tuples)):
    tuple1 = tuples[k]
    tuple2 = tuples2[k]
    for i in range(len(tuple1)):
        print('Testing {0} and {1}'.format(tuple1[i], tuple2[i]))
        for j in range(len(tuple1[i])):
            if tuple1[i][j] != tuple2[i][j]:
                print(tuple1[i], tuple2[i])
                break
            else:
                print('\t OK')
