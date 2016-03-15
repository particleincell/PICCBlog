# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 19:01:45 2016

@author: lbrieda
"""

import numpy as np
import pylab as pl
from random import random as rand

NP = 20000;		#number of samples
NUM_BINS = 500;	#number of uniques we would like

mass = 16*1.66e-27  #16amu
Kb = 1.3806e-23     #boltzmann constant
T = 300           #temperature
vth = np.sqrt(2*Kb*T/mass)		#thermal velocity

bins = np.zeros((NUM_BINS,1));

bin_min = 0;
bin_max = max(2,6*vth)
bin_max = 5000
dbin = (bin_max-bin_min)/NUM_BINS

def maxw1d(vth):
    M = 12			#parameter for Birdsall's method
    #sample Maxwellian using Birdsall's method
    sum = 0
    for i in range(M):
        sum = sum+rand()	
    return np.sqrt(0.5)*vth*(sum-M/2)/np.sqrt(M/12);

# sampling by taking three 1D distributiona
def maxw(vth):
    v1 = maxw1d(1)
    v2 = maxw1d(1)
    v3 = maxw1d(1)
    return vth* np.sqrt(v1*v1+v2*v2+v3*v3)

# sampling directly from distribution function	
def maxw_direct(vth):
    #pick random velocity in this bin
    #v = rand()*(dbin)+(v_min+(i-1)*dbin);
		
    #evaluate normalized distribution function
    fm = vel**2*np.exp(-vel**2/vth**2);
    return fm	

 
#storage for particles
vels = np.zeros((NP,3))
      
for s in range (NP):	
    
    mag = maxw(vth)
	
    #bin result, add to nearest
    bin = np.floor((mag-bin_min)/dbin+0.5)
    if (bin<NUM_BINS):
        bins[bin] = bins[bin] + 1
    
bins = bins/np.max(bins);	#%normalize bins

vel = np.arange(bin_min,bin_max,dbin)	#velocities for plotting
#fm = np.exp(-vel**2/vth**2)		#theoretical maxwellian
gm = vel**2*np.exp(-vel**2/vth**2)
print(np.max(gm))
gm = gm/np.max(gm)
#vel = vel/vth;			#normalize by thermal velocity

pl.figure(1)
pl.plot(vel,bins,'x');
pl.hold(True)
pl.plot(vel,gm,'r-');		#analytical 
pl.grid()
pl.hold(False)

pl.figure(2);
pl.semilogy(vel,bins,'x');
pl.hold(True)
pl.semilogy(vel,gm,'r-');
pl.ylim([1e-4,1]);
pl.grid()
pl.hold(False)
