#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 19:12:09 2023

@author: liang
"""
#importing the required modules and defining constants
import pyshtools

import matplotlib.pyplot as plt
#from pyshtools import constant, can't do for Mars+Venus
import pdb 
import numpy as np
import scipy.io as io
import sys
import time
import copy
from scipy.interpolate import UnivariateSpline


a = 1738000
gm = 4902.80011526323e9

a = 6378136.3
gm = 398600441500000.0

lmax_calc = 240 
#extracting the necessary data

clm0 = pyshtools.SHGravCoeffs.from_file('GO_CONS_GCF_2_TIM_R6.gfc', format='icgem',lmax=lmax_calc)
shape = pyshtools.SHCoeffs.from_file('Earth2014.RET2014.degree2160.bshc',format='bshc',lmax=lmax_calc)

shape_ret = copy.deepcopy(shape)
shape_coeff = shape_ret.to_array()

shape_grid0 = shape_ret.expand(grid='DH2',lmax=lmax_calc)
interface0 = shape_grid0.to_array()


interface0 = interface0 + a

interface0 = np.delete(interface0,0,1)
interface0 = np.delete(interface0,0,0)

# interface0 = nd.gaussian_filter(interface0, 2)

incilm = pyshtools.expand.SHExpandDH(interface0,lmax_calc=lmax_calc,sampling=2)
shape_ret = pyshtools.SHCoeffs.from_array(incilm,lmax=lmax_calc)#,a=a)

#basic data cleaning and pre-processing. Actual data cleaning/pre-processing takes 10x more lines (see mare_fill.py)

interface0 = shape.to_array()
bc = pyshtools.SHGravCoeffs.from_shape(shape_ret, rho=2600, gm=gm, lmax=lmax_calc,nmax=9)#3100
bc2 = bc.change_ref(r0=a,gm=gm,lmax=lmax_calc)
clm0 = clm0.change_ref(r0=a,gm=gm,lmax=lmax_calc)

clm = clm0 - bc2

cilm = clm.to_array()
cilm[0][0][0] = 1
cilm[:,1:5,:] = 0

L = np.arange(0,lmax_calc+1)
L1 = round(0.9*lmax_calc)#-25
L2 = lmax_calc

taper0 = (1-np.cos((L2-L)/(L2-L1)*np.pi))/2
taper = np.tile(taper0,(L2+1,1)).transpose()

taperf = np.zeros((2,L2+1,L2+1))
taperf[0,:,:] = taper
taperf[1,:,:] = taper

cilm[:,L1:,:] = taperf[:,L1:,:]*cilm[:,L1:,:]

#calculating gravity

r,t,p,earth_grav,pot = pyshtools.gravmag.MakeGravGridDH(cilm,gm,a,lmax=lmax_calc,a=a,f=0,omega=0,normal_gravity=1)


exit()
#plotting gravity

plt.figure(figsize=(20,10))
ax = plt.axes([0, 0.05, 0.9, 0.9 ])
ax.set_title('Earth Gravity Map (mGal)',fontsize=50)
fig = plt.imshow(earth_grav*1e5,cmap='jet',aspect='equal',vmin=-300,vmax=300)

cax = plt.axes([0.92, 0.05, 0.03,0.9 ])
cb = plt.colorbar(fig,cax=cax)
cb.ax.tick_params(labelsize=30)


ax.axis('off')

np.save('Earth_grav_demo.npy',earth_grav)



