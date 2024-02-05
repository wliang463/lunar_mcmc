#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Code to perform the Markov chain Monte Carlo method using the 
Metropolis-Hastings algorithm on GRAIL gravity data to constrain
the dimensions of the IBC bodies beneath the lunar surface

@author: Weigang Liang
"""


import pyshtools


import numpy as np
import scipy.io as io
import time
import copy
import sys
start_time = time.time()

def likelihood(td,bd,xs,ys,gravv):
    
    grav_p0 = forward_grav(td,bd-td,xs,ys,100)
    grav_p = np.mean(grav_p0[:,grav_col1:grav_col2]*1e5,axis=1)
     
    summ = (np.mean((gravv[grav_row1:grav_row2]+grav_p[grav_row1:grav_row2])**2))**0.5
    
    llh = np.exp(np.float128(-summ**2/297))    
 
    return llh, summ

#Information for the various gravity anomalies of the Moon. 
#The NW, NE, and S anomalies are in the third, fourth, and fifth indices, respectively.


grav_col1_0 = np.array((900,1185,700,120,1872))
grav_col2_0 = np.array((945,1251,850,200,1940))
grav_row1_0 = np.array((484,575,291,150,627))
grav_row2_0 = np.array((523,629,360,240,682))
y_cen_0 = np.array((501,595,333,196,641))
x_cen_0 = np.array((921,1221,776,160,1906))
ys_0 = np.array((350e3,344e3,702e3,256e3,336e3))


for ii in np.array([2,3,4]):

    numm = ii
    
    grav_col1 = grav_col1_0[numm];
    grav_col2 = grav_col2_0[numm];
    grav_row1 = grav_row1_0[numm];
    grav_row2 = grav_row2_0[numm];
    
    
    from_matlab = {}
    
    if ii==2:
        from_matlab = io.loadmat('nnw.mat')
    elif ii==3:
        from_matlab = io.loadmat('nne.mat')
    elif ii==4:
        from_matlab = io.loadmat('ns.mat')
        
    gravv0 = np.squeeze(from_matlab['grav2'])
    #normalize minimum gravity of region of interest to 0
    gravv = np.mean(gravv0[:,grav_col1:grav_col2]*1e5,axis=1)
    gravv = gravv - min(gravv[grav_row1:grav_row2])
    
    G = 6.67408e-11;
    drho = 400;#density contrast between mantle and IBC
    
    y_cen = y_cen_0[numm]
    x_cen = x_cen_0[numm]
    
    ys = ys_0[numm]
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
    
    for jj in np.arange(5):
    
        ano_num = jj;
            
        
        bar = 1e10
        xs_grid = []
        th_grid = []
        td_grid = []
        bar_grid = []
        lf = 500
        start_time = time.time()
        
        pass_num = 6
        
        g0 = 1.6230976292682058
        
        L = np.arange(0,lf+1)
        taper0 = 1/((L+1)*g0)
        taper = np.tile(taper0,(lf+1,1)).transpose()
        
        taperf = np.zeros((2,lf+1,lf+1))
        taperf[0,:,:] = taper
        taperf[1,:,:] = taper
        
        xs = np.random.uniform(50e3,200e3)
        td = np.random.uniform(20e3,30e3)
        bd = np.random.uniform(td,500e3)
        
        bd_c = np.array([]);
        td_c = np.array([]);
        xs_c = np.array([]);
        bars = np.array([]);
    
        
        total = 10000 #number of iterations
        
        acc = 0
        rej = 0
        
        bdd = 0
        tdd = 0
        xss = 0
        
        
        
        for step in np.arange(total):
            if step%1000 == 0:
                elapsed_time = time.time() - start_time
                print(elapsed_time)

            llho, summ1 = likelihood(td, bd, xs, ys, gravv)
            
            td1 = td+np.random.normal(0,1.5e3)
            bd1 = bd+np.random.normal(0,7e3)
            xs1 = xs+np.random.normal(0,7e3)
            llhn, summ2 = likelihood(td1, bd1, xs1, ys, gravv)
            
            judge = np.random.uniform(0, 1)
            
            if xs1 <= 0 or td1 >= 0 or td1 >= bd1:
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
                        print('wtf')
                        print(llhn)
                        print(llho)
                        print(judge)
                
        
            
           
        filename = 'dike_results_mcmc' + '_' + str(int(td)) + '_' + str(int(xs)) + '_' + str(int(bd)) + '_' + str(ano_num) + '.mat'
        filename2 = 'dike_results_mcmc' + '_' + str(int(td)) + '_' + str(int(xs)) + '_' + str(int(bd)) + '_' + str(ano_num) + '_bestfits.mat'
        
        io.savemat(filename,{'bd_c':bd_c,'xs_c':xs_c,'td_c':td_c,'acc':acc,'rej':rej,'bars':bars})
        
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
        
        io.savemat(filename2,{'tdd':tdd,'xss':xss,'bdd':bdd,'thh':thh1,'rr':rr1})
        
        
        tda = np.mean(td_c)
        bda = np.mean(bd_c)
        xsa = np.mean(xs_c)
        th1a = np.mean(th_c1)
        rr1a = np.mean(r_c1)
        
        td_state = str(r(tdd)/1e3) + ' (' + str(r(tda/1e3)) + '+' + str(r(rb(td_c)-tda)/1e3) + ',-' + str(r(tda-lb(td_c))/1e3) + ')'
        bd_state = str(r(bdd)/1e3) + ' (' + str(r(bda/1e3)) + '+' + str(r(rb(bd_c)-bda)/1e3) + ',-' + str(r(bda-lb(bd_c))/1e3) + ')'
        xs_state = str(r(xss)/1e3) + ' (' + str(r(xsa/1e3)) + '+' + str(r(rb(xs_c)-xsa)/1e3) + ',-' + str(r(xsa-lb(xs_c))/1e3) + ')'
        
        
        print(td_state)
        print(bd_state)
        print(xs_state)
        
        h_state = str(r(thh1)/1e3) + ' (' + str(r(th1a/1e3)) + '+' + str(r(rb(th_c1)-th1a)/1e3) + ',-' + str(r(th1a-lb(th_c1))/1e3) + ')'
        r_state = str(r(rr1)) + ' (' + str(r(rr1a)) + '+' + str(r(rb(r_c1)-rr1a)) + ',-' + str(r(rr1a-lb(r_c1))) + ')'
        
        print('height 1 = ' + h_state)
        print('ratio 1 = ' + r_state)
        
        
        
        print(str(r(bar)))
        
        
        from xlrd import open_workbook
        from xlutils.copy import copy as cpy
        
        rb = open_workbook("tri.xls")
        wb = cpy(rb)
        
        num = jj+1
        
        print('SHEET NUM AT ' + str(ano_num))
        s = wb.get_sheet(int(ii-2))
        s.write(num,1,td_state)
        s.write(num,2,bd_state)
        s.write(num,3,xs_state)
        s.write(num,4,h_state)
        s.write(num,5,r_state)
        s.write(num,6,r(bar))
        wb.save('tri.xls')
