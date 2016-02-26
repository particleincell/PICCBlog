# axisymmetric (RZ) particle in cell code example
#
# see https://www.particleincell.com/2015/rz-pic/ for more info
# simulates a simplistic ion source in which ions are
# produced from a volumetric source with constant electron and 
# neutral density (i.e. not taking into account avalanche ionization
# or source depletion)
#
# code illustrates velocity and position rotation in RZ
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

#simple Jacobian solver, does not do any convergence checking
def solvePotential(phi,max_it=100):
    
    #make copy of dirichlet nodes
    P = numpy.copy(phi)  
    
    g = numpy.zeros_like(phi)
    dz2 = dz*dz
    dr2 = dr*dr

    rho_e = numpy.zeros_like(phi)
    
    #set radia
    r = numpy.zeros_like(phi)
    for i in range(nz):
        for j in range(nr):
            r[i][j] = R(j)
    
    for it in range (max_it):
        #compute RHS
        #rho_e = QE*n0*numpy.exp(numpy.subtract(phi,phi0)/kTe)
        
#        for i in range(1,nz-1):
#            for j in range(1,nr-1):
                
#                if (cell_type[i,j]>0):
#                    continue
                
#                rho_e=QE*n0*math.exp((phi[i,j]-phi0)/kTe)
#                b = (rho_i[i,j]-rho_e)/EPS0;
#                g[i,j] = (b + 
#                         (phi[i,j-1]+phi[i,j+1])/dr2 +
#                         (phi[i,j-1]-phi[i,j+1])/(2*dr*r[i,j]) +
#                         (phi[i-1,j] + phi[i+1,j])/dz2) / (2/dr2 + 2/dz2)
#                         
 #               phi[i,j]=g[i,j]

        #compute electron term                                         
        rho_e = QE*n0*numpy.exp(numpy.subtract(P,phi0)/kTe)
        b = numpy.where(cell_type<=0,(rho_i - rho_e)/EPS0,0)

        #regular form inside        
        g[1:-1,1:-1] = (b[1:-1,1:-1] + 
                       (phi[1:-1,2:]+phi[1:-1,:-2])/dr2 +
                       (phi[1:-1,0:-2]-phi[1:-1,2:])/(2*dr*r[1:-1,1:-1]) +
                       (phi[2:,1:-1] + phi[:-2,1:-1])/dz2) / (2/dr2 + 2/dz2)
        
        #neumann boundaries
        g[0] = g[1]       #left
        g[-1] = g[-2]     #right
        g[:,-1] = g[:,-2] #top
        g[:,0] = g[:,1]
        
        #dirichlet nodes
        phi = numpy.where(cell_type>0,P,g)
        
    return phi

#computes electric field                    
def computeEF(phi,efz,efr):
    
    #central difference, not right on walls
    efz[1:-1] = (phi[0:nz-2]-phi[2:nz+1])/(2*dz)
    efr[:,1:-1] = (phi[:,0:nr-2]-phi[:,2:nr+1])/(2*dr)
    
    #one sided difference on boundaries
    efz[0,:] = (phi[0,:]-phi[1,:])/dz
    efz[-1,:] = (phi[-2,:]-phi[-1,:])/dz
    efr[:,0] = (phi[:,0]-phi[:,1])/dr
    efr[:,-1] = (phi[:,-2]-phi[:,-1])/dr
    
def plot(ax,data,scatter=False):
    pl.sca(ax)
    pl.cla()
    cf = pl.contourf(pos_z, pos_r, numpy.transpose(data),8,alpha=.75,linewidth=1,cmap='jet')
    #cf = pl.pcolormesh(pos_z, pos_r, numpy.transpose(data))
    if (scatter):
        ax.hold(True);
        (ZZ,RR)=pl.meshgrid(pos_z,pos_r)
        ax.scatter(ZZ,RR,c=numpy.transpose(cell_type),cmap='jet')
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
sub = (pl.subplot(211),pl.subplot(212))

#solve potential
phi = solvePotential(phi,1000)
computeEF(phi,efz,efr);

plot(sub[0],efz)
plot(sub[1],phi,scatter=True)               
#Q = pl.quiver(pos_z, pos_r, numpy.transpose(efz), numpy.transpose(efr),units='xy')
pl.draw()
      

#this will block execution until figure is closed
pl.show()

