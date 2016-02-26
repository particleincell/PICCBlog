# -*- coding: utf-8 -*-
"""
demo of loading maxwellian VDF

@author: lbrieda
"""
import numpy as np
import pylab as pl
from random import random as rand

NP = 50000;		#number of samples
NUM_BINS = 500;	#number of uniques we would like


mass = 16*1.66e-27  #16amu
Kb = 1.3806e-23     #boltzmann constant
T = 1000           #temperature
vdrift = 100       #drift velocity
vth = np.sqrt(2*Kb*T/mass)		#thermal velocity

bins = np.zeros((NUM_BINS,1));

bin_min = 0;
bin_max = max(2*vdrift,10*vth)
dbin = (bin_max-bin_min)/NUM_BINS

def maxw1d(vth):
    M = 12			#parameter for Birdsall's method
    #sample Maxwellian using Birdsall's method
    sum = 0
    for i in range(M):
        sum = sum+rand()	
    return np.sqrt(0.5)*vth*(sum-M/2)/np.sqrt(M/12);
 
#storage for particles
vels = np.zeros((NP,3))
      
for s in range (NP):	
    
    v1 = maxw1d(vth)
    v2 = maxw1d(vth)
    v3 = maxw1d(vth)+vdrift
    
    mag = np.sqrt(v1*v1+v2*v2+v3*v3)
	
    #bin result, add to nearest
    bin = np.floor((mag-bin_min)/dbin+0.5)
    if (bin<NUM_BINS):
        bins[bin] = bins[bin] + 1
    
    #save for later
    vels[s] = [v1,v2,v3]

#compute temperature
#part 1), get average velocity
vel_ave = [np.sum(vels[:,0]),
           np.sum(vels[:,1]),
           np.sum(vels[:,2])]
vel_ave = np.divide(vel_ave,NP)
print("vel_ave = [%.3g %.3g %.3g]"%(vel_ave[0],vel_ave[1],vel_ave[2]))    

#part 2, compute mv^2
v2_sum = 0
for s in range (NP):
    dv = vels[s] - vel_ave
    v2 = dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2]
    v2_sum = v2_sum + v2  

v2_ave = np.divide(v2_sum,NP)
temp_computed = v2_ave * mass/(3.*Kb)
print("Computed temperature %.3g"%temp_computed)


bins = bins/np.max(bins);	#%normalize bins

vel = np.arange(bin_min,bin_max,dbin)	#velocities for plotting
#fm = np.exp(-vel**2/vth**2)		#theoretical maxwellian
gm = (vel)**2*np.exp(-(vel-vdrift)**2/vth**2)
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
