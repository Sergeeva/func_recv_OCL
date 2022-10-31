#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec
from matplotlib import cm

from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

import csv

#########################################################################################################################
#########################################################################################################################

if (__name__=='__main__'):

#    def read_csv(filename, shape):
#        source_data = np.zeros(shape)
#        with open(filename, 'rt') as f:
#            rows = csv.reader(f)
#            for rno, row in enumerate(rows):
#                # print(row)
#                source_data[rno, :] = row
#                # print(perf_fbf[rno, :])
#        return source_data
#
#    shape = (128, 128)
##
    
    filename = 'func_recv_2d__result.dat'
    frecv_2d_res = np.loadtxt(filename)

    #print(frecv_2d_res)
    
    MX, MY = np.meshgrid(range(16*4), range(16*4))
    fig1 = plt.figure(1)
#    ax1 = fig1.gca(projection='3d')
    ax1 = fig1.add_subplot(1, 1, 1, projection='3d')
#    surf1 = ax1.plot_surface(MX, MY, frecv_2d_res , cmap = cm.jet, rstride = 1, cstride = 1, linewidth = 0, antialiased = False)   
    ax1.scatter(MX, MY, frecv_2d_res, c='r', marker='*')
#    surf1 = ax1.plot_wireframe(MX, MY, frecv_2d_res , rstride = 1, cstride = 1) #, rcount = 128, ccount = 128)   

    filename = 'func_recv_2d__input.dat'
    frecv_2d_in = np.loadtxt(filename)

    #print(frecv_2d_in)
    
    KX, KY = np.meshgrid(range(0,16*4,4), range(0,16*4,4))
#    fig2 = plt.figure(2)
#    ax2 = fig2.gca(projection='3d')
#    ax2 = fig1.add_subplot(1, 2, 2, projection='3d')
    #surf2 = ax2.plot_surface(KX, KY, frecv_2d_in , cmap = cm.jet, rstride = 1, cstride = 1, linewidth = 0, antialiased = False)   
    ax1.scatter(KX, KY, frecv_2d_in, c='b', marker='*')
    #surf2 = ax2.plot_wireframe(KX, KY, frecv_2d_in , rstride = 1, cstride = 1) #, rcount = 32, ccount = 32)
    

   


    plt.show()
else:
    print(__name__, 'is imported')
