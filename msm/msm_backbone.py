import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
from glob import glob
import pickle
import numpy as np
import matplotlib.pyplot as plt
import sys
from os.path import join, isfile
import pickle
import multiprocessing as mp
import copy
from pyemma.coordinates import source, load


if __name__ == '__main__':





    # LOAD TRAJECTORIES
    print('Loading Trajectories')
    traj_list = ['backbone/2agy-310k-1atm-prod1.2-backbone-2ns.dcd'] #glob('./backbone/.dcd')

    inp = coor.source(traj_list, feat)
    feature_data = inp.get_output()

    # TICA
    # Use lag time of 1ns
    tica_obj = coor.tica(inp, lag=1000, kinetic_map=True, var_cutoff=.9)
    print('TICA dimension ', tica_obj.dimension())
    # Sort out eigenvalue sign
    for i in range(2):
        if tica_obj.eigenvectors[0, i] > 0:
            tica_obj.eigenvectors[:, i] *= -1

    Y = tica_obj.get_output() # get tica coordinates
    print('number of trajectories = ', np.shape(Y)[0])
    print('number of frames = ', np.shape(Y)[1])
    print('number of dimensions = ',np.shape(Y)[2])



    # # INVESTIGATE LAG TIMES
    # print('Producing lag-time chart')
    # # inp = coor.source(copy.deepcopy(traj_list), copy.deepcopy(feat))
    # #
    # # data = [(100, copy.deepcopy(inp)), (200, copy.deepcopy(inp))] #[(lag, coor.source(traj_list, feat)) for lag in lags]
    # # dim = 10

    # data = [100,200]
    # p = mp.Pool(2)
    # print(p.map(mp_worker, data))



    # plt.plot(lags, times)
    # plt.show()