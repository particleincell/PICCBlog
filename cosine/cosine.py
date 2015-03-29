# -*- coding: utf-8 -*-
"""
sampling from a cosine distribution
https://www.particleincell.com/2015/cosine-distribution/

@author: lubos brieda
"""
from random import random
import math
import matplotlib.pyplot as plt

plt.close("all")
u = []
v = []
w = []

count = []
bin_theta = []
nbins = 91
delta_bin = 90/(nbins-1)

#surface properties
tang1 = [1,0,0]
tang2 = [0,1,0]
norm =  [0,0,1]

#initialize data
for bin in range(0,nbins):
    count.append(0);
    bin_theta.append(bin*delta_bin)

# sample random points
for it in range(1,100000):
    sin_theta = random()    
    cos_theta = math.sqrt(1-sin_theta*sin_theta)

    #random in plane angle
    psi = random()*2*math.pi;
    
    #three vector components
    a = sin_theta*math.cos(psi)
    b = sin_theta*math.sin(psi)
    c = cos_theta
    
    #multiply by corresponding directions
    v1 = [a*n for n in tang1]
    v2 = [b*n for n in tang2]
    v3 = [c*n for n in norm]
    
    #add up to get velocity, vel=v1+v2+v3
    vel = []
    for i in range(0,3):
        vel.append(v1[i]+v2[i]+v3[i])
    
    #update histogram
    theta = math.acos(cos_theta)*180/math.pi;
    bin = int( (theta)/delta_bin)
    count[bin] += 1
    
    #add every 100th particle to the visualization list
    if (it%100==0):
        u.append(vel[0])
        v.append(vel[1])
        w.append(vel[2])
    
#plot results
fig = plt.figure()
ax = fig.add_subplot(211, projection='3d')
ax.scatter(u, v, w, c='r', marker='.')
ax.set_xlabel('u')
ax.set_ylabel('v')
ax.set_zlabel('w')

#normalize data
c0 = count[0]
count[:] = [float(c) / c0 for c in count]
cos = [math.cos(th*math.pi/180) for th in bin_theta]

#add xy plot
ax2 = fig.add_subplot(212)
ax2.plot(bin_theta,count,c='r');
ax2.plot(bin_theta,cos,c='k')
ax2.set_xlabel('angle')
ax2.set_ylabel('normalized count')
plt.show()
    