## Gets the bonds, angles and dihedrals in the forcefield.

from parmed.amber import AmberParm
from parmed.tools import strip, outparm

parm = AmberParm('2agy_final_min.prmtop')
strip(parm, "!(:400,419)").execute()
outparm(parm, '2agy-as1.prmtop').execute()

for x in parm.bonds:
    print(x.atom1.idx, x.atom2.idx)

