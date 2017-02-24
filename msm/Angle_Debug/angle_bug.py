import pyemma.coordinates as coor
import numpy as np
import pickle

feat = coor.featurizer('act_site.pdb')
feat.add_angles(np.array([[0,1,2], [1,2,3]]), deg=False, cossin=False, periodic=False)
inp = coor.source('prod0-as1-aligned.dcd', feat)

# Save the feature description for comparison later
pickle.dump(inp.describe(), open('nocos_desc.p', 'wb'))


# Comparisons
if True:
    fixed = pickle.load(open('fixed_desc.p', 'rb'))     # rad = rad.reshape(rad.shape[0], rad.shape[1]*rad.shape[2])
    broken = pickle.load(open('broken_desc.p', 'rb'))   # rad = rad.reshape(functools.reduce(lambda x, y: x * y, rad.shape),)
    nocos = pickle.load(open('nocos_desc.p', 'rb'))     # rad = rad.reshape(functools.reduce(lambda x, y: x * y, rad.shape),)
                                                        # & cossin=False

    fixed = fixed[1:]
    broken = broken[1:]
    nocos = nocos[1:]

    # Check consistency with cossin = False
    for idx in range(len(fixed)):
        if not nocos[idx//2][7:].strip() == fixed[idx][11:-1].strip():
            print idx, nocos[idx//2][7:].strip(), fixed[idx][11:-1].strip()
        else:
            print idx, ' - OK'

    # Check consistency between fixed and broken descriptions
    for idx in range(len(fixed)):
        if not broken[idx] == fixed[idx]:
            print idx, broken[idx], fixed[idx]
        else:
            print idx, ' - OK'