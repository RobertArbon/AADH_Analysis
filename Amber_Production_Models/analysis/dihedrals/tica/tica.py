"""Reduce dimensionality with tICA

msmbuilder autogenerated template version 2
created 2017-05-23T16:38:49.125259
please cite msmbuilder in any publications

"""

from msmbuilder.io import load_trajs, save_trajs, save_generic
from msmbuilder.decomposition import tICA

## Load
tica = tICA(n_components=5, lag_time=10, kinetic_mapping=True)
meta, ftrajs = load_trajs("ftrajs")

## Fit
tica.fit(ftrajs.values())

## Transform
ttrajs = {}
for k, v in ftrajs.items():
    ttrajs[k] = tica.partial_transform(v)

## Save
save_trajs(ttrajs, 'ttrajs', meta)
save_generic(tica, 'tica.pickl')