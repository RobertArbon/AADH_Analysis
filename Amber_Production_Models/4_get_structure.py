from msmbuilder.io import load_meta, save_generic
from simtk.openmm import app
from simtk import unit
import parmed as pmd
from os.path import join

# # Get topology
# meta = load_meta()
# top_fn = join('proc_traj',meta['top_fn'][0])
#
# # Set up OpenMM system
# ffdir = '/home/robert/Research/AADH/MD/common'
# fffiles = ['top_all36_prot.rtf', 'par_all36_prot_mod.prm',
#             'aadh_parameters.prm', 'aadh_topologies.rtf']
# ffpaths = [join(ffdir, x) for x in fffiles]
# print(ffpaths)
# prmtop = app.AmberPrmtopFile(top_fn)
# # forcefield = app.ForceField(*ffpaths)

