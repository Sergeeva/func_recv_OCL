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
    
    filename = 'kern4.dat'
    kern4_dat = np.loadtxt(filename)

##    print(kern4_dat)
    
    MX, MY = np.meshgrid(range(-8, 8+1), range(-8, 8+1))
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111, projection='3d')
    #surf1 = ax1.plot_wireframe(MX, MY, kern4_dat, rstride = 1, cstride = 1) #, rcount = 128, ccount = 128)
    surf1 = ax1.plot_surface(MX, MY, kern4_dat, cmap=cm.coolwarm, antialiased=False)

    # Plot a sin curve using the x and y axes.
#    x = np.linspace(-8, 1, 8+1)
#    y = [0, 0, 0, 0, 0, 0, 0, 0, 0]
#    ax1.plot(x, y, zs=0, zdir='z')



    plt.show()
else:
    print(__name__, 'is imported')
