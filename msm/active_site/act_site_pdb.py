import MDAnalysis
import numpy.linalg
import numpy as np
import pickle
from MDAnalysis.analysis.align import *
import os.path

pdb = '../common/2agy_final.pdb'
psf = '../common/2agy_final.psf'


trj = MDAnalysis.Universe(psf, pdb)
act_site = trj.select_atoms('segid BT1 and (resid 39 or resid 58)')

act_site.write('act_site.pdb')



