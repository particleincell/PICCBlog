# axisymmetric (RZ) particle in cell code example
#
# see https://www.particleincell.com/2015/rz-pic/ for more info
# simulates a simplistic ion source in which ions are
# produced from a volumetric source with constant electron and 
# neutral density (i.e. not taking into account avalanche ionization
# or source depletion)
#
# code illustrates velocity and position rotation in RZ
# step 1: domain setup
#
# requires numpy, scipy, pylab, and mathplotlib

import numpy
import pylab as pl
import math
from random import (seed,random)

def XtoL(pos):
    lc = [pos[0]/dz, pos[1]/dr]
    return lc

def Pos(lc):
    pos = [lc[0]*dz,lc[1]*dr]
    return pos
 
def R(j):
    return j*dr
 
    
def plot(ax,data,scatter=False):
    pl.sca(ax)
    pl.cla()
    cf = pl.contourf(pos_z, pos_r, numpy.transpose(data),8,alpha=.75,linewidth=1,cmap='jet')
    #cf = pl.pcolormesh(pos_z, pos_r, numpy.transpose(data))
    if (scatter):
        ax.hold(True);
        (ZZ,RR)=pl.meshgrid(pos_z,pos_r)
        ax.scatter(ZZ,RR,c=numpy.transpose(cell_type),cmap='winter')
    ax.set_yticks(pos_r)
    ax.set_xticks(pos_z)
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    pl.xlim(min(pos_z),max(pos_z))
    pl.ylim(min(pos_r),max(pos_r))
    ax.grid(b=True,which='both',color='k',linestyle='-')
    ax.set_aspect('equal', adjustable='box')
  #  pl.colorbar(cf,ax=pl.gca(),orientation='horizontal',shrink=0.75, pad=0.01)

    
#---------- INITIALIZATION ----------------------------------------

pl.close('all')
seed()
       
# allocate memory space
nz = 35
nr = 12
dz = 1e-3
dr = 1e-3    
dt = 5e-9

QE = 1.602e-19
AMU =  1.661e-27
EPS0 = 8.854e-12

charge = QE
m = 40*AMU  #argon ions  
qm = charge/m                   
spwt = 50

#solver parameters
n0 = 1e12
phi0 = 100
phi1 = 0
kTe = 5

phi = numpy.zeros([nz,nr])
efz = numpy.zeros([nz,nr])
efr = numpy.zeros([nz,nr])
rho_i = numpy.zeros([nz,nr])
den = numpy.zeros([nz,nr])

# ---- sugarcube domain --------------------
cell_type = numpy.zeros([nz,nr]);
tube1_radius = 6*dr;
tube1_length = 0.01;
tube1_aperture_rad = 4*dr;
tube2_radius = tube1_radius+dr;
tube2_length = tube1_length+2*dz;
tube2_aperture_rad = 3*dr;
[tube_i_max, tube_j_max] = map(int, XtoL([4*dz, tube1_radius]))

for i in range(0,nz):
    for j in range(0,nr):
        pos = Pos([i,j])  # node position
        
        #inner tube
        if ((i==0 and pos[1]<tube1_radius) or 
        (pos[0]<=tube1_length and pos[1]>=tube1_radius and pos[1]<tube1_radius+0.5*dr) or
        (pos[0]>=tube1_length and pos[0]<tube1_length+0.5*dz and 
            pos[1]>=tube1_aperture_rad and pos[1]<tube1_radius) ) :
            cell_type[i][j]=1
            phi[i][j] = phi0
        
        if ((pos[0]<=tube2_length and pos[1]>=tube2_radius and pos[1]<tube2_radius+0.5*dr) or
        (pos[0]>=tube2_length and pos[0]<=tube2_length+1.5*dz and 
            pos[1]>=tube2_aperture_rad and pos[1]<=tube2_radius) ) :
            cell_type[i][j]=2
            phi[i][j] = phi1
                       
#----------- COMPUTE NODE VOLUMES ------------------------
node_volume = numpy.zeros([nz,nr])
for i in range(0,nz):
    for j in range(0,nr):
        j_min = j-0.5
        j_max = j+0.5
        if (j_min<0): j_min=0
        if (j_max>nr-1): j_max=nr-1
        a = 0.5 if (i==0 or i==nz-1) else 1.0
        #note, this is r*dr for non-boundary nodes
        node_volume[i][j] = a*dz*(R(j_max)**2-R(j_min)**2)  

#create an array of particles
particles = []
    
#counter for fractional particles    
mpf_rem = numpy.zeros([nz,nr])
rho_i = numpy.zeros([nz,nr])

lambda_d = math.sqrt(EPS0*kTe/(n0*QE))
print ("Debye length is %.4g, which is %.2g*dz"%(lambda_d,lambda_d/dz))
print ("Expected ion speed is %.2f m/s"%math.sqrt(2*phi0*qm))

#positions for plotting
pos_r = numpy.linspace(0,(nr-1)*dr,nr)
pos_z = numpy.linspace(0,(nz-1)*dz,nz)
fig1 = pl.figure(num=None, figsize=(10, 10), dpi=80, facecolor='w', edgecolor='k')
sub = (pl.subplot(111))
    
#----------- END OF MAIN LOOP ------------------------
plot(sub,cell_type,scatter=True)               
pl.draw()
      

#this will block execution until figure is closed
pl.show()

