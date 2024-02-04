#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 23:14:25 2022

@author: liang
"""


#tts: bottom depth, width, top depth, ellipse number
#center row+1, center column+1, anomaly length, grav col1, 2, grav row1, 2,
#anomaly number

import pyshtools


import numpy as np
import scipy.io as io
import time
import copy
import sys
start_time = time.time()

import pandas as pd







def mcmc(iterations):

    #known dimensions of anomaly


    
    grav_col1 = 120
    grav_col2 = 200
    grav_row1 = 150
    grav_row2 = 240
    y_cen = 196
    x_cen = 160
    ys = 256e3
    
    
    from_matlab = {}
    from_matlab = io.loadmat('nne.mat')

        
    gravv0 = np.squeeze(from_matlab['grav2'])
    #get min to 0
    gravv = np.mean(gravv0[:,grav_col1:grav_col2]*1e5,axis=1)
    gravv = gravv - min(gravv[grav_row1:grav_row2])
    
    # if np.size(sys.argv) < 6:
    #     #gravv = np.mean(gravv0[:,900:950]*1e5,axis=1)+14
    # else:
    #     #gravv = np.mean(gravv0[:,900:950]*1e5,axis=1)+int(sys.argv[5])
    
    #36,37,34
    G = 6.67408e-11;
    drho = 400;#Hess, Parmentier 1995
 
    y = np.array((-ys/2,ys/2)); 

    
    topo0 = 1738000*np.ones((1002,2004));
    
    lat = np.transpose(np.tile(np.linspace(0,180-180/1002,1002),(2004,1)))*np.pi/180;
    lon = np.tile(np.linspace(0,360-360/2004,2004),(1002,1))*np.pi/180;
    xp = topo0*np.sin(lat)*np.cos(lon);
    yp = topo0*np.sin(lat)*np.sin(lon);
    zp = topo0*np.cos(lat);
    
    latc = lat[y_cen-1,x_cen-1];
    lonc = lon[y_cen-1,x_cen-1];
    z_trans = np.array(((np.cos(lonc),np.sin(lonc),0),
        (-np.sin(lonc), np.cos(lonc), 0),
        (0, 0, 1)))
    y_trans = np.array(((np.cos(latc), 0, -np.sin(latc)),
        (0, 1, 0),
        (np.sin(latc), 0, np.cos(latc))))
    trans = np.matmul(y_trans,z_trans);
            
    xr = np.reshape(xp,(1,2004*1002)); yr = np.reshape(yp,(1,2004*1002)); zr = np.reshape(zp,(1,2004*1002));
    cr = np.matmul(trans,np.squeeze(np.stack((xr,yr,zr))));
    xp0 = np.reshape(cr[0,:],(1002,2004));
    yp0 = np.reshape(cr[1,:],(1002,2004));
    zp0 = np.reshape(cr[2,:],(1002,2004));
    
    xp = np.reshape(xp0[:,x_cen-1],(1002,1))
    yp = np.reshape(yp0[:,x_cen-1],(1002,1))
    zp = np.reshape(zp0[:,x_cen-1],(1002,1))
    
    
    def forward_grav(td,th,xs,ys,div):
        grav = np.zeros((1002,1));
        grav_f= np.zeros((1002,2004))
    
        
        x_s = np.flip((np.arange(div)+1)*xs/div/2)
        z_s = td+(np.arange(div+1))*th/div
    
        
        for n in np.arange(div):
            x = np.array((-x_s[n],x_s[n]));    
            z = np.array((1738e3-z_s[n],1738e3-z_s[n+1]))
        
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
        
                    
        
        grav_f = G*drho*np.repeat(grav,2004,axis=1);
        return grav_f
    
 
    def likelihood(td,bd,xs,ys,gravv):#no bd
    
        grav_p0 = forward_grav(td,bd-td,xs,ys,100)
        grav_p = np.mean(grav_p0[:,grav_col1:grav_col2]*1e5,axis=1)#900:950
        #gravv is assumed already mean-ed
         
        #summ = np.mean(abs(gravv[grav_row1:grav_row2]+grav_p[grav_row1:grav_row2]))#484:523
        summ = (np.mean((gravv[grav_row1:grav_row2]+grav_p[grav_row1:grav_row2])**2))**0.5
        
        llh = np.exp(np.float128(-summ**2/297))
        #np.sqrt(sum((p1-p2)**2)/501)*1e5*2; p1 and p2 are half of nnw's :,1620
        
        #modified for poor fits
        
        #llh = summ**-2/2
     
        return llh, summ
        
    bar = 1e10

    lf = 500
    start_time = time.time()
        
    g0 = 1.6230976292682058
    
    L = np.arange(0,lf+1)
    taper0 = 1/((L+1)*g0)
    taper = np.tile(taper0,(lf+1,1)).transpose()
    
    taperf = np.zeros((2,lf+1,lf+1))
    taperf[0,:,:] = taper
    taperf[1,:,:] = taper
    
    # a = 500; b = -350
    
                
    # import pdb 
    # pdb.set_trace()
    xs = np.random.uniform(50e3,200e3)
    td = np.random.uniform(20e3,50e3)
    bd = np.random.uniform(td,100e3)
    
    bd_c = np.array([]);
    td_c = np.array([]);
    xs_c = np.array([]);
    bars = np.array([]);

    
    total = iterations
    
    acc = 0
    rej = 0
    
    bdd = 0
    tdd = 0
    xss = 0
    
    #MCMC algorithm
    
    for step in np.arange(total):
        # if step%1000 == 0 and step//1000>0:
        #     elapsed_time = time.time() - start_time
        #     print(str(step) + ' interations, it has been ' + str(elapsed_time//3600)
        #            + ' hr')
        # llho, summ1 = likelihood(td, xs, far_bg, mtseff0)
        llho, summ1 = likelihood(td, bd, xs, ys, gravv)
        
        td1 = td+np.random.normal(0,1.5e3)
        bd1 = bd+np.random.normal(0,7e3)
        xs1 = xs+np.random.normal(0,7e3)
        llhn, summ2 = likelihood(td1, bd1, xs1, ys, gravv)
        # llhn, summ2 = likelihood(td1, xs1, far_bg, mtseff0)
        
        judge = np.random.uniform(0, 1)
        
        if xs1 <= 0 or td1 <= 0 or td1 >= bd1 or bd1 >= 150e3:
        # if td1 < 0 or xs1 < 0:
            llhn = 0
            
        test = llhn/llho
        
        if test < judge:
            rej+=1
        else:
            acc+=1
            bd = copy.deepcopy(bd1)
            xs = copy.deepcopy(xs1)
            td = copy.deepcopy(td1)
            
            bd_c = np.append(bd_c,bd)
            xs_c = np.append(xs_c,xs)
            td_c = np.append(td_c,td)
            bars = np.append(bars,summ2)
            
            if summ2 < bar:
                bar = copy.deepcopy(summ2)
                xss = copy.deepcopy(xs1)
                bdd = copy.deepcopy(bd1)
                tdd = copy.deepcopy(td1)
                
                if xss < 0 or bdd < 0:
                    print('Error. Please notify me via email.')
                    print(llhn)
                    print(llho)
                    print(judge)
            
    
        
       
    # filename = 'dike_results_mcmc' + '_' + str(int(td)) + '_' + str(int(xs)) + '_' + str(int(bd)) + '_' + str(ano_num) + '.mat'
    # filename2 = 'dike_results_mcmc' + '_' + str(int(td)) + '_' + str(int(xs)) + '_' + str(int(bd)) + '_' + str(ano_num) + '_bestfits.mat'
    
    # io.savemat(filename,{'bd_c':bd_c,'xs_c':xs_c,'td_c':td_c,'acc':acc,'rej':rej,'bars':bars})
    # io.savemat(filename,{'xs_c':xs_c,'td_c':td_c,'acc':acc,'rej':rej})
    
    from math import log10, floor
    def r(x, sig=3):
        if x == 0:
            return 0
        else:
            return round(x, sig-int(floor(log10(abs(x))))-1)
    
    def lb(chain):
        siz = chain.shape[0]
        c_sort = np.sort(chain)
        return c_sort[round((siz-1)*0.16)]
    
    def rb(chain):
        siz = chain.shape[0]
        c_sort = np.sort(chain)
        return c_sort[round((siz-1)*0.84)]
    
            
    
    thh1 = bdd-tdd
    th_c1 = bd_c-td_c
    

    rr1 = copy.deepcopy(thh1/xss)
    r_c1 = th_c1/xs_c
    
    # io.savemat(filename2,{'tdd':tdd,'xss':xss,'bdd':bdd,'thh':thh1,'rr':rr1})
    
    
    tda = np.mean(td_c)
    bda = np.mean(bd_c)
    xsa = np.mean(xs_c)
    th1a = np.mean(th_c1)
    rr1a = np.mean(r_c1)
    
    td_state = str(r(tdd)/1e3) + ' (' + str(r(tda/1e3)) + '+' + str(r(rb(td_c)-tda)/1e3) + ',-' + str(r(tda-lb(td_c))/1e3) + ')'
    bd_state = str(r(bdd)/1e3) + ' (' + str(r(bda/1e3)) + '+' + str(r(rb(bd_c)-bda)/1e3) + ',-' + str(r(bda-lb(bd_c))/1e3) + ')'
    xs_state = str(r(xss)/1e3) + ' (' + str(r(xsa/1e3)) + '+' + str(r(rb(xs_c)-xsa)/1e3) + ',-' + str(r(xsa-lb(xs_c))/1e3) + ')'
    
    
    print('Top Depth: ' + td_state + ' km')
    print('Width: ' + xs_state + ' km')
    
    h_state = str(r(thh1)/1e3) + ' (' + str(r(th1a/1e3)) + '+' + str(r(rb(th_c1)-th1a)/1e3) + ',-' + str(r(th1a-lb(th_c1))/1e3) + ')'
    
    print('Thickness = ' + h_state + ' km')    
    
    print('RMS: ' + str(r(bar)) + ' mGal')
    
         
    elapsed_time = time.time() - start_time
    print('Total runtime: ' + str(round(elapsed_time/60,3))
           + ' min')


    # Read existing excel file
    df = pd.read_excel('tri.xlsx')

    num = 0
    
    # Adding data
    df.iloc[num, 1] = td_state
    df.iloc[num, 2] = xs_state
    df.iloc[num, 3] = h_state
    df.iloc[num, 4] = r(bar)
    
    
    # Save the updated DataFrame to excel
    df.to_excel('tri.xlsx', index=False)
    
    
    # calculating the MCMC errorbar
    inds = np.argsort(bars)

    
    td_c = td_c[inds]
    bd_c = bd_c[inds]
    xs_c = xs_c[inds]
    
    
    num = round(inds.shape[0]/2)
    
    profs = np.zeros((num,1002))
    
    
    for jj in np.arange(num):

    
        td = td_c[jj]
        bd = bd_c[jj]
        xs = xs_c[jj]
        
            
        grav_p0 = forward_grav(td,bd-td,xs,ys,100)
        grav_p01 = np.mean(grav_p0[:,grav_col1:grav_col2]*1e5,axis=1)#900:950
        # std_p01 = np.std(grav_p0[:,grav_col1:grav_col2]*1e5,axis=1)#900:950
        #gravv is assumed already mean-ed
         
        # spec2 = grav_p[grav_row1:grav_row2];
        
        profs[jj,:] = grav_p01
    
    np.save('mcmc_error.npy',profs)
        
