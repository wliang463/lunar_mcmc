#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 22:56:27 2022

@author: liang
"""



import pyshtools

import matplotlib.pyplot as plt
#from pyshtools import constant, can't do for Mars+Venus
import pdb 
import numpy as np
import scipy.io as io
import scipy.ndimage as nd
import sys
import time
import copy
from scipy.interpolate import UnivariateSpline

def random_coeff(power):
    
    l = power.shape[0]
    
    coeff = np.zeros((2,l,l))
    
    for x in range(l):
        
        nums = np.random.randn(2*x)
        coeff_t = nums/(sum(nums**2)**0.5)*power[x]**0.5#**0.5/(2*x+1)**0.5
        
        
        #coeff[0,x,:x+1] = coeff_t[:x+1]
        #coeff[1,x,:x] = coeff_t[x+1:]
        
        
        coeff[0,x,:x] = coeff_t[:x]
        coeff[1,x,:x] = coeff_t[x:]
    
    power2 = pyshtools.spectralanalysis.spectrum(coeff)
        
    return coeff,power2

def cuntopo(topo,lmax_calc):
    
    siz = topo.shape   
    
    topof = np.zeros((2003,4005))
    
    topof[235:1770] = topo
    topof[0:235] = 1738e3
    topof[1770:] = 1738e3
    
    L1 = 0
    L2_1 = int(siz[0]*0.25)
    L2_2 = int(siz[0]*0.75)
    L2 = lmax_calc
    
    L = np.arange(0,L2_1)
    
    # taper1 = -(np.cos((L2-L)/(L2-L1)*np.pi))
    # taper0 = np.flip(taper1)
     
    # taper00 = np.tile(taper0,(siz[1],1)).transpose()
    # taper11 = np.tile(taper1,(siz[1],1)).transpose()
    
    # taperf = np.ones(siz)
    # taperf[0:L2_1,:] = np.flipud(taper00)
    # taperf[L2_2+1:,:] = taper11
    
    topo[0:L2_1,:] = 1737e3
    topo[L2_2+1:,:] = 1737e3
    
    #topo = taperf*topo
    
    return topof
    
    
now = 'earth'#replaced with mars/venus if necesssary
fill = 'fill'
#Jeff: the for loop for the Bouguer elevation part starts at line 161
start_time = time.time()
# your code

mode = 'spectra'

if now in 'moon':
    rho_c = 2550
    rho_d = -450
    a = 1738000
    gm = 4902.80011526323e9
    l = 1200#500
    nmax = 14
    shape = pyshtools.SHCoeffs.from_file('MoonTopo2600p.shape.sh',lmax=l)
 #   shape = pyshtools.SHCoeffs.from_file('lro_ltm05_2050_sha.tab.txt',lmax=l)
    omega = 0
    f=0
    near = 1735000
    far = 1742500
    dh = far-near
    mass = 7.3457892e22
elif now in 'mars':
    rho_c = 2580
    rho_d = -520
    mass = 6.417e23
    l = 100
    shape = pyshtools.SHCoeffs.from_file('MarsTopo2600.shape.sh',lmax=l)
    gm = 4.2828376e13
    nmax = 7
    a = 3395428
    f = 0
    omega = 0
    dh = 0
elif now in 'venus':
    rho_c = 2580
    rho_d = -520
    mass = constant.Venus.mass_venus.value
    shape = pyshtools.SHCoeffs.from_file('VenusTopo719.shape')
    l = 90
    nmax = 7
    a = constant.Venus.r_venus.value
    f = 0
    omega = 0
elif now in 'earth':

    a = 6378136.3
    gm = 398600441500000.0
    l = 500#500
    nmax = 7
    shape = pyshtools.SHCoeffs.from_file('Earth2014.RET2014.degree2160.bshc',format='bshc',lmax=l)
 #   shape = pyshtools.SHCoeffs.from_file('lro_ltm05_2050_sha.tab.txt',lmax=l)
    omega = 0
    f=0


#original
lmax_calc = 500#450#500#l
  
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


temp = shape.to_array()
temp[:,31:,:] = 0
interface2 = pyshtools.expand.MakeGridDH(temp,lmax=l,sampling=2)

interface2 = interface2 + a
# shape2 = pyshtools.SHCoeffs.from_array(temp,lmax=lmax_calc)
# shape_grid2 = shape2.expand(grid='DH2',lmax=lmax_calc)
# interface2 = shape_grid2.to_array()


#rho = rho_d+np.zeros(interface2.shape)

size = interface2.shape
#interface2 = np.concatenate((interface[:,:size[1]//4],interface[:,size[1]//4:3*size[1]//4],interface[:,3*size[1]//4:]),axis=1)


b_step = 100


pass_num =  5#20#5

shape_sur = pyshtools.SHCoeffs.from_file('Earth2014.SUR2014.degree2160.bshc',format='bshc',lmax=l)

shape_coeff2 = shape_sur.to_array()

shape_grid2 = shape_sur.expand(grid='DH2',lmax=lmax_calc)
interface02 = shape_grid2.to_array()


interface02 = interface02 + a

interface02 = np.delete(interface02,0,1)
interface02 = np.delete(interface02,0,0)

# interface0 = nd.gaussian_filter(interface0, 2)

incilm2 = pyshtools.expand.SHExpandDH(interface02,lmax_calc=lmax_calc,sampling=2)
shape_sur = pyshtools.SHCoeffs.from_array(incilm2,lmax=lmax_calc)#,a=a)

added = 0#need 1000 for diffusion

mmin = np.min(interface2) 
mmax = np.max(interface2)


clm = pyshtools.SHGravCoeffs.from_file('GO_CONS_GCF_2_TIM_R6.gfc', format='icgem')

mod_clm0 = clm.change_ref(r0=a,gm=gm,lmax=250)

mod_clm1 = mod_clm0.to_array()
temp = np.zeros(shape_coeff.shape)
temp[:,:251,:251] = mod_clm1

mod_clm0 = pyshtools.SHGravCoeffs.from_array(temp,r0=a,gm=gm)

bc = pyshtools.SHGravCoeffs.from_shape(shape_ret, rho=2600., gm=gm, lmax=lmax_calc)
bc = bc.change_ref(r0=a)

bc2 = pyshtools.SHGravCoeffs.from_shape(shape_sur, rho=500., gm=gm, lmax=lmax_calc)
bc2 = bc.change_ref(r0=a)

mod_clm0 = mod_clm0 - bc - bc2
mod_clm0 = mod_clm0.to_array()
mod_clm0[0][0][0] = 1

mod_clm0[:,1:pass_num,:] = 0

temp_clm = mod_clm0
mmin = np.min(interface2)
mmax = np.max(interface2)
potential = np.zeros([interface2.shape[0],interface2.shape[1],np.arange(mmin,mmax,b_step).shape[0]])
ind = 0
new_potential3 = np.zeros(interface2.shape)
bouguer3 = np.zeros(interface2.shape)
bg4 = pyshtools.SHGravCoeffs.from_array(temp_clm,gm,clm.r0,omega=0,lmax=l)

for i in np.arange(mmin,mmax,b_step):
    bg3 = bg4.change_ref(r0=i,gm=gm,lmax=l)#bg4
    temp = bg3.to_array()

    potential2 = pyshtools.expand.MakeGridDH(temp,lmax=l,sampling=2)
    new_potential3[(interface2>=i)*(interface2<=i+b_step)] = potential2[(interface2>=i)*(interface2<=i+b_step)]
    r,t,p,gravb,pot = pyshtools.gravmag.MakeGravGridDH(temp_clm,gm,a,lmax=l,a=i,f=0,omega=0,normal_gravity=1)
    bouguer3[(interface2>=i)*(interface2<=i+b_step)] = gravb[(interface2>=i)*(interface2<=i+b_step)]
    if ind == 0:
        new_potential3[interface2<=i] = potential2[interface2<=i]
        bouguer3[interface2<=i] = gravb[interface2<=i]
    elif ind == np.arange(mmin,mmax,b_step).shape[0]-1:
        new_potential3[interface2>=i] = potential2[interface2>=i]
        bouguer3[interface2>=i] = gravb[interface2>=i]
    ind+=1

cilm = pyshtools.expand.SHExpandDH(new_potential3,lmax_calc=lmax_calc,sampling=2)

cilm[0][0][0] = 1
#cilm0 = cilm[:]


#low pass
ratio = 6378/1737
L1 = int(25*ratio)
L2 = int(35*ratio)

L = np.arange(0,L2+1)


taper0 = (1-np.cos((L2-L)/(L2-L1)*np.pi))/2
taper = np.tile(taper0,(L2+1,1)).transpose()

taperf = np.zeros((2,L2+1,L2+1))
taperf[0,:,:] = taper
taperf[1,:,:] = taper

cilm[:,L1:L2,:L2+1] = taperf[:,L1:L2,:]*cilm[:,L1:L2,:L2+1]
cilm[:,L2+1:,:] = 0

#high pass

# ratio = 6378/1737
# L1 = int(5*ratio)
# L2 = int(15*ratio)

# L = np.arange(0,L2+1)


# taper0 = (np.cos((L2-L)/(L2-L1)*np.pi))/2
# taper = np.tile(taper0,(L2+1,1)).transpose()

# taperf = np.zeros((2,L2+1,L2+1))
# taperf[0,:,:] = taper
# taperf[1,:,:] = taper

# cilm[:,L1:L2,:L2+1] = taperf[:,L1:L2,:]*cilm[:,L1:L2,:L2+1]
# cilm[:,:L1-1,:] = 0

# #cosine taper at end of coeffs
# L1 = lmax_calc-100
# L2 = lmax_calc

# L = np.arange(0,L2+1)


# taper0 = (1-np.cos((L2-L)/(L2-L1)*np.pi))/2
# taper = np.tile(taper0,(L2+1,1)).transpose()

# taperf = np.zeros((2,L2+1,L2+1))
# taperf[0,:,:] = taper
# taperf[1,:,:] = taper

# cilm[:,L1:L2,:L2+1] = taperf[:,L1:L2,:]*cilm[:,L1:L2,:L2+1]

ea = np.array((3,0.5,0))#nnw

# ea = np.array((0,-0.7,0))#FSP_NW_R
# ea = np.array((0,-0.85,0))#FSP_NW_L
# ea = np.array((-2.5,-0.66,0))#NSW
# ea = np.array((0,-0.64,0))#FN_W
# ea = np.array((0,-0.2,0))#inter

# ea = np.array((-3.1,-2.21,0))#far sp, center

cilm[0][0][0] = 1

temp_clm2 = pyshtools.rotate.SHRotateRealCoef(cilm,ea,pyshtools.rotate.djpi2(lmax_calc))

temp_clm2[0][0][0] = 1

# clm2 = pyshtools.SHGravCoeffs.from_array(temp_clm2,gm,a,omega=0,lmax=lmax_calc)
r,t,p,grav2,pot = pyshtools.gravmag.MakeGravGridDH(temp_clm2,gm,a,lmax=lmax_calc,a=a,f=0,omega=0,normal_gravity=1)



r,t,p,grav,pot = pyshtools.gravmag.MakeGravGridDH(cilm,gm,a,lmax=l,a=a+0,f=0,omega=0,normal_gravity=1)

clm = pyshtools.SHGravCoeffs.from_array(cilm,gm,a+0,omega=0,lmax=lmax_calc)

tensorr = clm.tensor(a=a+0,f=0,lmax=lmax_calc)
tensorr.compute_eigh()
eig = tensorr.eighh.data

    
eig22 = np.concatenate((eig[:,size[1]//4:3*size[1]//4],eig[:,3*size[1]//4:],eig[:,:size[1]//4]),axis=1)
grav22 = np.concatenate((grav[:,size[1]//4:3*size[1]//4],grav[:,3*size[1]//4:],grav[:,:size[1]//4]),axis=1)

    
elapsed_time = time.time() - start_time
print(elapsed_time)
sys.exit()

#cosine taper for crustal thickness

# R1=1700; R2=1670;
# L=[1:110];
# Rscl = (R1/R2).^(L+1);
# L1=30; L2=110;
# f = (1+cos((L-L1)/(L2-L1)*pi))/2;
# f(1:L1)=1; f(L2:110)=0;
# Rscl2=(Rscl.*f + 1*(1-f));
# plot(L,Rscl2);

#spectra time
#temp_clm1 = mod_clm0
    
theta0 = np.pi/20#og 20
lmax = 19#og 19
#n_dh = 2*(lmax*4+1)
#OG maps

maskn0 = np.zeros(interface2.shape)
maskf0 = np.zeros(interface2.shape)

#maskn0[400:500,1800:1900] = 1
#maskn0[300:400,1660:1760] = 1
maskn0[475:575,1900:2000] = 1
maskf0[300:400,1200:1300] = 1
maskf0[430:530,600:700] = 1

mask_far = np.zeros(interface2.shape)
mask_far[253:747,506:1494] = 1

#eig is >0.98,>0.99

tapersn_og, eigenvalues_nog = pyshtools.spectralanalysis.SHReturnTapersMap(maskn0, 40, sampling=2,ntapers=3)#4,3
mtsen_og, sdnog = pyshtools.spectralanalysis.SHMultiTaperMaskSE(cilm, tapersn_og)#temp_clm1

tapersf_og, eigenvalues_fog = pyshtools.spectralanalysis.SHReturnTapersMap(maskf0, 80, sampling=2,ntapers=2)#23,22
mtsef_og, sdfog = pyshtools.spectralanalysis.SHMultiTaperMaskSE(cilm, tapersf_og)

mtsef_og2, sdf = pyshtools.spectralanalysis.SHMultiTaperMaskSE(cilm+cilm_topo, tapersf_og)


#MATLAB ROI maps

from_matlab = {}
from_matlab = io.loadmat('masks.mat')
maskn = from_matlab['maskn']
maskf = from_matlab['maskf']

all_mask = maskn[0,0]+maskn[0,1]+maskn[0,2]+maskn[0,3]+maskf[0,0]+maskf[0,1]+maskn0+maskf0+maskf[0,2]

#theta0 = 100/interface.shape[1]*np.pi
#lmax = int(np.floor(1*np.pi/theta0-1))
theta0 = np.pi/20
lmax = 19



#clm1 = pyshtools.SHGravCoeffs.from_array(cilm,gm,a,omega=0,lmax=l)
#clm2 = clm1.change_ref(r0=a+0,gm=gm,lmax=l)
#cilm1 = clm2.to_array()

nla = (501-450)/501*90;nlon = 1850/2004*360;#og
fla = (501-350)/501*90;flon = 1250/2004*360#og


#fla = 36;flon = 207;#"smooth" patch

ea = np.array((2.5,1.91,0))#far north
ea = np.array((-3.1,-1.91,0))#far sp
ea = np.array((3,0.5,0))#nnw

ea = np.array((0,-0.7,0))#FSP_NW_R
ea = np.array((0,-0.85,0))#FSP_NW_L
ea = np.array((-2.5,-0.66,0))#NSW
ea = np.array((0,-0.64,0))#FN_W
ea = np.array((0,-0.2,0))#inter

ea = np.array((-3.1,-2.21,0))#far sp, center

plt.figure()
plt.plot(np.arange(mtsen.shape[0]),np.log10(mtsen),'k',linewidth=3)
plt.plot(np.arange(mtsen.shape[0]),np.log10(mtsef_1km),'tab:gray',linewidth=3)
plt.plot(np.arange(mtsen.shape[0]),np.log10(mtsef_400),'k',linewidth=1.5,linestyle='dashed')
plt.plot(np.arange(mtsef_tp1.shape[0]),np.log10(mtsef_tp2),'tab:gray',linewidth=1.5,linestyle='dashed')
plt.plot(np.arange(mtsen.shape[0]),np.log10(mtsef_diff),'k',linewidth=1.5,linestyle='dotted')
plt.plot(np.arange(mtsef_tp1.shape[0]),np.log10(mtsef_final),'tab:gray',linewidth=1.5,linestyle='dotted')
plt.legend(['nearside mare region', 'flooded farside model','400 kg/m^3 density contrast','2-km pre-mare mare','runaway diffusion','400 contrast + 2 km pre-mare'],fontsize=12)
plt.ylim((-17.5,-12))
plt.xlim((30,420))
plt.xlabel('Spherical harmonic degree',fontsize=16)
plt.ylabel('log(Power)',fontsize=16)

temppp = pyshtools.rotate.SHRotateRealCoef(cilm,ea,pyshtools.rotate.djpi2(lmax_calc))
if now in 'moon':
    tapersnn0, eigenvalues_n0 = pyshtools.spectralanalysis.SHReturnTapersMap(maskn[0,0], 40, sampling=2,ntapers=4)#4,2
    mtsenn0, sdn0 = pyshtools.spectralanalysis.SHMultiTaperMaskSE(temp_clm1, tapersnn0)
    
    tapersnn1, eigenvalues_n1 = pyshtools.spectralanalysis.SHReturnTapersMap(maskn[0,1], 40, sampling=2,ntapers=2)#2,1
    mtsenn1, sdn1 = pyshtools.spectralanalysis.SHMultiTaperMaskSE(temp_clm1, tapersnn1)
    
    tapersnn2, eigenvalues_n2 = pyshtools.spectralanalysis.SHReturnTapersMap(maskn[0,2], 40, sampling=2,ntapers=2)#2,1
    mtsenn2, sdn2 = pyshtools.spectralanalysis.SHMultiTaperMaskSE(temp_clm1, tapersnn2)
    
    tapersnn3, eigenvalues_n3 = pyshtools.spectralanalysis.SHReturnTapersMap(maskn[0,3], 40, sampling=2,ntapers=3)#3,2
    mtsenn3, sdn3 = pyshtools.spectralanalysis.SHMultiTaperMaskSE(temp_clm1, tapersnn3)
    
    Rn = np.mean(interface0[maskn[0,1]==1])
    mlist = (a/Rn)**np.arange(0,l+1)
    marray = np.tile(mlist,(l+1,1)).transpose()
    
    multn = np.zeros((2,l+1,l+1))
    multn[0,:,:] = marray
    multn[1,:,:] = marray
    
    temp_clm2 = mod_clm0*multn
#    temp_clm_2 = temp_clm_*multn
    
    
    
    mtsenn1_2, sdn = pyshtools.spectralanalysis.SHMultiTaperMaskSE(temp_clm2, tapersnn1)
    mtsenn1_0, sdn = pyshtools.spectralanalysis.SHMultiTaperMaskSE(bg_clm_grail, tapersnn1)
    
#    mtsenn1_2_2, sdn = pyshtools.spectralanalysis.SHMultiTaperMaskSE(temp_clm_2, tapersnn1_2)
    
    #far side
    
    tapersff0, eigenvalues_f0 = pyshtools.spectralanalysis.SHReturnTapersMap(maskf[0,0], 80, sampling=2,ntapers=7)#8,7
    mtseff0, sdf0 = pyshtools.spectralanalysis.SHMultiTaperMaskSE(cilm, tapersff0)
    #west 5,4
    #east 4,3
    
    tapersff1, eigenvalues_f1 = pyshtools.spectralanalysis.SHReturnTapersMap(maskf[0,1], 80, sampling=2,ntapers=4)#5,4
    mtseff1, sdf1 = pyshtools.spectralanalysis.SHMultiTaperMaskSE(cilm, tapersff1)
    
    tapersff2, eigenvalues_f2 = pyshtools.spectralanalysis.SHReturnTapersMap(maskf[0,2], 80, sampling=2,ntapers=16)#17,16
    mtseff2, sdf2 = pyshtools.spectralanalysis.SHMultiTaperMaskSE(cilm, tapersff2)
    
    #26,22 central rim
    
    eigenvalues = np.concatenate((eigenvalues_n0,eigenvalues_n1,eigenvalues_n2,eigenvalues_n3,eigenvalues_f0,eigenvalues_f1,eigenvalues_f2))

    
    
#    cilm2 = (cilm0 + bg_clm_grail)*multf
#    cilm2 = bg_clm_grail*multf
    
    lf = lmax_calc-lmax
    
    from_matlab = {}
    from_matlab = io.loadmat('grav_p.mat')
    grav_p = from_matlab['grav_p']
    tempp = pyshtools.expand.SHExpandDH(grav_p,lmax_calc=lf,sampling=2)
    tempp[:,1:pass_num,:] = 0
    
    Rf = np.mean(interface0[maskf[0,0]==1])
    mlist = (a/Rf)**np.arange(0,lf+1)
    marray = np.tile(mlist,(lf+1,1)).transpose()
    
    multf = np.zeros((2,lf+1,lf+1))
    multf[0,:,:] = marray
    multf[1,:,:] = marray
    
    temp_clm_p = tempp*multf
    
#    ind = 0
#    new_potential4 = np.zeros(interface2.shape)
#    bouguer5 = np.zeros(interface2.shape)
#    bg6 = pyshtools.SHGravCoeffs.from_array(tempp,gm,bg.r0,omega=0,lmax=l)
#
#    for i in np.arange(mmin,mmax,b_step):
#        bg5 = bg6.change_ref(r0=i,gm=gm,lmax=l)#bg4
#        temp = bg5.to_array()
#
#        potential3 = pyshtools.expand.MakeGridDH(temp,lmax=l,sampling=2)
#        new_potential4[(interface2>=i)*(interface2<=i+b_step)] = potential3[(interface2>=i)*(interface2<=i+b_step)]
#        r,t,p,gravb,pot = pyshtools.gravmag.MakeGravGridDH(temp_clm,gm,a,lmax=l,a=i,f=0,omega=0,normal_gravity=1)
#        bouguer5[(interface2>=i)*(interface2<=i+b_step)] = gravb[(interface2>=i)*(interface2<=i+b_step)]
#        if ind == 0:
#            new_potential4[interface2<=i] = potential3[interface2<=i]
#            bouguer5[interface2<=i] = gravb[interface2<=i]
#        elif ind == np.arange(mmin,mmax,b_step).shape[0]-1:
#            new_potential4[interface2>=i] = potential3[interface2>=i]
#            bouguer5[interface2>=i] = gravb[interface2>=i]
#        ind+=1
#    
#    temp_clm_p = pyshtools.expand.SHExpandDH(new_potential4,lmax_calc=lf,sampling=2)
    temp_clm_p[0][0][0] = 1
    
    far_bg = random_coeff(mtsef_og)
    testt = pyshtools.spectralanalysis.spectrum(temp_clm_p)
    f0_syn = pyshtools.spectralanalysis.spectrum(temp_clm_p + far_bg)
    
    tempp2 = tempp*multf
    
    def forward_grav(td,th,xs):
        G = 6.67408e-11;
        drho = 300;
        grav = np.zeros((1002,2004));
            
        topo0 = 1738000*np.ones((1002,2004));
        
        lat = np.transpose(np.tile(np.linspace(0,180-180/1002,1002),(2004,1)))*np.pi/180;
        lon = np.tile(np.linspace(0,360-360/2004,2004),(1002,1))*np.pi/180;
        xp = topo0*np.sin(lat)*np.cos(lon);
        yp = topo0*np.sin(lat)*np.sin(lon);
        zp = topo0*np.cos(lat);
        
        latc = lat[120-1,940-1];
        lonc = lon[120-1,940-1];
        z_trans = np.array(((np.cos(lonc),np.sin(lonc),0),
            (-np.sin(lonc), np.cos(lonc), 0),
            (0, 0, 1)))
        y_trans = np.array(((np.cos(latc), 0, -np.sin(latc)),
            (0, 1, 0),
            (np.sin(latc), 0, np.cos(latc))))
        trans = np.matmul(y_trans,z_trans);
                
#        td = 43e3;
#        th = 9e3;
#        xs = 36e3;
        ys = 130e3;
               
        x = np.array((-xs,xs));
        y = np.array((-ys,ys));
        z = np.array((1738e3-td,1738e3-td-th))
    
        xr = np.reshape(xp,(1,2004*1002)); yr = np.reshape(yp,(1,2004*1002)); zr = np.reshape(zp,(1,2004*1002));
        cr = np.matmul(trans,np.squeeze(np.stack((xr,yr,zr))));
        xp = np.reshape(cr[0,:],(1002,2004));
        yp = np.reshape(cr[1,:],(1002,2004));
        zp = np.reshape(cr[2,:],(1002,2004));
        
        for i in np.arange(2):
            dxi = x[i]-xp;
            for j in np.arange(2):
                dyj = y[j]-yp;
                for k in np.arange(2):
                    mu = (-1)**i*(-1)**j*(-1)**k;
                    dzk = z[k]-zp;
                    R = (dxi**2 + dyj**2 + dzk**2)**0.5;
                    grav = grav + mu*(dzk*(np.arctan(dxi*dyj/(dzk*R)))
                        - dxi*np.log(R+dyj)
                        - dyj*np.log(R+dxi));
        
                    
        
        grav = G*drho*grav;
        return grav
    
#    th = 10e3
    far_bg = random_coeff(mtsef_og)
    bar = 1e10
    xs_grid = []
    th_grid = []
    td_grid = []
    bar_grid = []
    
    
    for bd in np.arange(20,50,2)*1e3:
        for xs in np.arange(4,34,2)*1e3:
            for td in np.arange(2,36,2)*1e3:
                if bd <= td:
                    continue
                grav_p = forward_grav(td,bd-td,xs/2)
                tempp = pyshtools.expand.SHExpandDH(grav_p,lmax_calc=lf,sampling=2)
                tempp[0][0][0] = 1
                tempp[:,1:pass_num,:] = 0
                f0_syn = pyshtools.spectralanalysis.spectrum(tempp + far_bg)
                testt = (sum((np.log(f0_syn[100:250])-np.log(mtseff0[100:250]))**2))**0.5#original lower limit 50
                if testt < bar:
                    bar = copy.deepcopy(testt)
                    xss = copy.deepcopy(xs)
                    tdd = copy.deepcopy(td)
                    thh = copy.deepcopy(bd-td)
                    f0_synn = copy.deepcopy(f0_syn)
                if testt < 6.7962598763087305:
                    xs_grid.append(xs)
                    th_grid.append(bd-td)
                    td_grid.append(td)
                    bar_grid.append(testt)
                    
#    for th in np.arange(0,55,5)*1e3:
#        for xs in np.arange(0,55,5)*1e3:
#            for td in np.arange(10,60,10)*1e3:
#                grav_p = forward_grav(td,th,xs/2)
#                tempp = pyshtools.expand.SHExpandDH(grav_p,lmax_calc=lf,sampling=2)
#                tempp[0][0][0] = 1
#                tempp[:,1:pass_num,:] = 0
#                f0_syn = pyshtools.spectralanalysis.spectrum(tempp + far_bg)
#                testt = (sum((np.log(f0_syn[50:250])-np.log(mtseff0[50:250]))**2))**0.5
#                if testt < bar:
#                    bar = copy.deepcopy(testt)
#                    xss = copy.deepcopy(xs)
#                    tdd = copy.deepcopy(td)
#                    thh = copy.deepcopy(th)
#                    f0_synn = copy.deepcopy(f0_syn)
                

   
        
#plt.plot(np.arange(l),mtsef)
#
###return grav,grav2
    
    plt.figure()
    fig = plt.plot(np.arange(mtsenn0.shape[0]),np.log10(mtsen_og))
    plt.hold(True)
    plt.plot(np.arange(mtsenn0.shape[0]),np.log10(mtsef_og2))
    plt.legend(['near og','far og'])
    plt.xlabel('Spherical harmonic degree')
    plt.ylim((-20,-2))
    plt.xlim((0,500))
    plt.savefig('og comp.png',dpi=1000)
    
    plt.figure()
    fig = plt.plot(np.arange(mtsenn0.shape[0]),np.log10(mtsen_og))
    plt.hold(True)
    plt.plot(np.arange(mtsenn0.shape[0]),np.log10(mtsenn0))
    plt.legend(['near og','near nw'])
    plt.xlabel('Spherical harmonic degree')
    plt.ylim((-20,-2))
    plt.xlim((0,500))
    plt.savefig('near nw.png',dpi=1000)
    
    plt.figure()
    fig = plt.plot(np.arange(mtsenn0.shape[0]),np.log10(mtsen_og))
    plt.hold(True)
    plt.plot(np.arange(mtsenn0.shape[0]),np.log10(mtsenn1))
    plt.legend(['near og','near ne'])
    plt.xlabel('Spherical harmonic degree')
    plt.ylim((-20,-2))
    plt.xlim((0,500))
    plt.savefig('near ne.png',dpi=1000)
    
    plt.figure()
    fig = plt.plot(np.arange(mtsenn0.shape[0]),np.log10(mtsen_og))
    plt.hold(True)
    plt.plot(np.arange(mtsenn0.shape[0]),np.log10(mtsenn2))
    plt.legend(['near og','near se'])
    plt.xlabel('Spherical harmonic degree')
    plt.ylim((-20,-2))
    plt.xlim((0,500))
    plt.savefig('near se.png',dpi=1000)
    
    plt.figure()
    fig = plt.plot(np.arange(mtsenn0.shape[0]),np.log10(mtsen_og))
    plt.hold(True)
    plt.plot(np.arange(mtsenn0.shape[0]),np.log10(mtsenn3))
    plt.legend(['near og','near s'])
    plt.xlabel('Spherical harmonic degree')
    plt.ylim((-20,-2))
    plt.xlim((0,500))
    plt.savefig('near s.png',dpi=1000)
    
    plt.figure()
    fig = plt.plot(np.arange(mtsenn0.shape[0]),np.log10(mtsef_og))
    plt.hold(True)
    plt.plot(np.arange(mtsenn0.shape[0]),np.log10(mtseff0))
    plt.legend(['far og','far north'])
    plt.xlabel('Spherical harmonic degree')
    plt.ylim((-20,-2))
    plt.xlim((0,500))
    plt.savefig('far north.png',dpi=1000)
    
    plt.figure()
    fig = plt.plot(np.arange(mtsenn0.shape[0]),np.log10(mtsef_og))
    plt.hold(True)
    plt.plot(np.arange(mtsenn0.shape[0]),np.log10(mtseff1))
    plt.legend(['far og','far sp'])
    plt.xlabel('Spherical harmonic degree')
    plt.ylim((-20,-2))
    plt.xlim((0,500))
    plt.savefig('far sp.png',dpi=1000)
    
    plt.figure()
    fig = plt.plot(np.arange(mtsenn0.shape[0]),np.log10(mtsef_og))
    plt.hold(True)
    plt.plot(np.arange(mtsenn0.shape[0]),np.log10(mtseff2))
    plt.legend(['far og','far sp rim'])
    plt.xlabel('Spherical harmonic degree')
    plt.ylim((-20,-2))
    plt.xlim((0,500))
    plt.savefig('far sp rim.png',dpi=1000)
    
    
    
    plt.figure()
    fig = plt.imshow(grav2*1e5,cmap='RdYlBu_r',aspect='equal',vmin=-200,vmax=200)
    cb = plt.colorbar(orientation="horizontal")
    cb.set_label('40 high pass, mGal')
    plt.savefig('bouguer 40 high.png',dpi=1000)
    
    plt.figure()
    fig = plt.imshow(eig2,cmap='RdYlBu_r',aspect='equal',vmin=-30,vmax=30)
    cb = plt.colorbar(orientation="horizontal")
    cb.set_label('40 high pass')
    plt.savefig('eig 40 high.png',dpi=1000)
#def func(x, a, b, c):
#    return a * np.exp(-b * x) + c
    
    plt.figure()
    x = np.arange(20,mtsen_og.shape[0])
    y = mtsen_og[20:]
    
    y_spl = UnivariateSpline(x,y)
    plt.semilogy(x,y,'ro',label = 'data')
    x_range = np.linspace(x[0],x[-1],1000)
    plt.semilogy(x_range,y_spl(x_range))
    
    plt.figure()    
    y_spl_2d = y_spl.derivative(n=2)

    plt.plot(x_range,y_spl_2d(x_range))
    plt.savefig('near og inflection.png',dpi=1000)
    
    plt.figure()
    x = np.arange(20,mtsen_og.shape[0])
    y = mtsef_og[20:]
    
    y_spl = UnivariateSpline(x,y,s=0,k=4)
    
    
    x_range = np.linspace(x[0],x[-1],1000)
    
    y_spl_2d = y_spl.derivative(n=2)

    plt.plot(x_range,y_spl_2d(x_range))
    plt.savefig('far og inflection.png',dpi=1000)
    
    plt.figure()
    x = np.arange(20,mtsen_og.shape[0])
    y = mtsenn0[20:]
    
    y_spl = UnivariateSpline(x,y,s=0,k=4)
    
    
    x_range = np.linspace(x[0],x[-1],1000)
    y_spl_2d = y_spl.derivative(n=2)

    plt.plot(x_range,y_spl_2d(x_range))
    plt.savefig('near nw inflection.png',dpi=1000)
    
    plt.figure()
    x = np.arange(20,mtsen_og.shape[0])
    y = mtsenn1[20:]
    
    y_spl = UnivariateSpline(x,y,s=0,k=4)
    
    
    x_range = np.linspace(x[0],x[-1],1000)
    y_spl_2d = y_spl.derivative(n=2)

    plt.plot(x_range,y_spl_2d(x_range))
    plt.savefig('near ne inflection.png',dpi=1000)
    
    plt.figure()
    x = np.arange(20,mtsen_og.shape[0])
    y = mtsenn2[20:]
    
    y_spl = UnivariateSpline(x,y,s=0,k=4)
    
    
    x_range = np.linspace(x[0],x[-1],1000)
    y_spl_2d = y_spl.derivative(n=2)

    plt.plot(x_range,y_spl_2d(x_range))
    plt.savefig('near se inflection.png',dpi=1000)
    
    plt.figure()
    x = np.arange(20,mtsen_og.shape[0])
    y = mtsenn3[20:]
    
    y_spl = UnivariateSpline(x,y,s=0,k=4)
    
    
    x_range = np.linspace(x[0],x[-1],1000)
    y_spl_2d = y_spl.derivative(n=2)

    plt.plot(x_range,y_spl_2d(x_range))
    plt.savefig('near s inflection.png',dpi=1000)
    
    plt.figure()
    x = np.arange(20,mtsen_og.shape[0])
    y = mtseff0[20:]
    
    y_spl = UnivariateSpline(x,y,s=0,k=4)
    
    
    x_range = np.linspace(x[0],x[-1],1000)
    y_spl_2d = y_spl.derivative(n=2)

    plt.plot(x_range,y_spl_2d(x_range))
    plt.savefig('far n inflection.png',dpi=1000)
    
    plt.figure()
    x = np.arange(20,mtsen_og.shape[0])
    y = mtseff1[20:]
    
    y_spl = UnivariateSpline(x,y,s=0,k=4)
    
    
    x_range = np.linspace(x[0],x[-1],1000)
    y_spl_2d = y_spl.derivative(n=2)

    plt.plot(x_range,y_spl_2d(x_range))
    plt.savefig('far sp inflection.png',dpi=1000)
    
    
def func(x, a, b, c):

    return a * x**-b + c



#df: a*(-b)*x**(-b-1)
#d2f: a*(-b)*(-b-1)*x**(-b-2)

#y2 = np.log10(mtsen_og[19:-1])

from scipy.optimize import curve_fit
poptn_og, pcovn_og = curve_fit(func, np.arange(20,mtsen_og.shape[0]), np.log10(mtsen_og[20:]),p0 = np.asarray([22,0.4,-18]))
poptnn0, pcovnn0 = curve_fit(func, np.arange(20,mtsen_og.shape[0]), np.log10(mtsenn0[20:]),p0 = np.asarray([22,0.4,-18]))
poptnn1, pcovnn1 = curve_fit(func, np.arange(20,mtsen_og.shape[0]), np.log10(mtsenn1[20:]),p0 = np.asarray([22,0.4,-18]))
poptnn2, pcovnn2 = curve_fit(func, np.arange(20,mtsen_og.shape[0]), np.log10(mtsenn2[20:]),p0 = np.asarray([22,0.4,-18]))
poptnn3, pcovnn3 = curve_fit(func, np.arange(20,mtsen_og.shape[0]), np.log10(mtsenn3[20:]),p0 = np.asarray([22,0.4,-18]))

poptf_og, pcovf_og = curve_fit(func, np.arange(20,mtsen_og.shape[0]), np.log10(mtsef_og[20:]),p0 = np.asarray([22,0.4,-18]))
poptff0, pcovff0 = curve_fit(func, np.arange(20,mtsen_og.shape[0]), np.log10(mtseff0[20:]),p0 = np.asarray([22,0.4,-18]))
poptff1, pcovff1 = curve_fit(func, np.arange(20,mtsen_og.shape[0]), np.log10(mtseff1[20:]),p0 = np.asarray([22,0.4,-18]))

popt_all = np.column_stack((poptn_og,poptnn0,poptnn1,poptnn2,poptnn3,poptf_og,poptff0,poptff1))
pcov_all = np.column_stack((np.sqrt(np.diag(pcovn_og)),np.sqrt(np.diag(pcovnn0)),np.sqrt(np.diag(pcovnn1)),np.sqrt(np.diag(pcovnn2)),np.sqrt(np.diag(pcovnn3)),np.sqrt(np.diag(pcovf_og)),np.sqrt(np.diag(pcovff0)),np.sqrt(np.diag(pcovff1))))







