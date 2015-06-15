# axisymmetric (RZ) particle in cell code example
#
# simulates a simplistic ion source in which ions are
# produced from a volumetric source with constant electron and 
# neutral density (i.e. not taking into account avalanche ionization
# or source depletion)
#
# code illustrates velocity and position rotation in RZ

# part 2: adds potential solver
#
# requires numpy, scipy, pylab, and mathplotlib

import numpy
import pylab as pl
import math
from random import (seed,random)
from matplotlib.colors import LogNorm

def XtoL(pos):
    lc = [pos[0]/dz, pos[1]/dr]
    return lc

def Pos(lc):
    pos = [lc[0]*dz,lc[1]*dr]
    return pos
 
def R(j):
    return j*dr
 
def gather(data,lc):
    i = numpy.trunc(lc[0])
    j = numpy.trunc(lc[1])
    di = lc[0] - i
    dj = lc[1] - j
    return  (data[i][j]*(1-di)*(1-dj) +
          data[i+1][j]*(di)*(1-dj) + 
          data[i][j+1]*(1-di)*(dj) + 
          data[i+1][j+1]*(di)*(dj)) 
    
def scatter(data,lc,value):
    i = numpy.trunc(lc[0])
    j = numpy.trunc(lc[1])
    di = lc[0] - i
    dj = lc[1] - j
            
    data[i][j] += (1-di)*(1-dj)*value
    data[i+1][j] += (di)*(1-dj)*value
    data[i][j+1] += (1-di)*(dj)*value
    data[i+1][j+1] += (di)*(dj)*value
       
#particle definition
class Particle:
    def __init__(self,pos,vel):
        self.pos=[pos[0],pos[1],0]    #xyz position
        self.vel=[vel[0],vel[1],vel[2]]

# --- helper functions ----
def sampleIsotropicVel(vth):
    #pick a random angle
    theta = 2*math.pi*random()
    
    #pick a random direction for n[2]
    R = -1.0+2*random()
    a = math.sqrt(1-R*R)
    n = (math.cos(theta)*a, math.sin(theta)*a, R)
    
    #pick maxwellian velocities
    vm = numpy.zeros(3)
    vm[0:3] = vth*(2*(random()+random()+random()-1.5))
    
    vel = (n[0]*vm[0], n[1]*vm[1], n[2]*vm[2]) 
    return vel

#simple Jacobian solver, does not do any convergence checking
def solvePotential(phi, rho_i):
    
    #reset potential, some instability is driving the solution
    #towards divergence if this is not used
    #phi[:][:] = numpy.where(cell_type>0,P,0)
  
    #make copy of dirichlet nodes
    P = numpy.copy(phi)  
    
    t1 = numpy.zeros_like(phi)
    t2 = numpy.zeros_like(phi)
    dz2 = dz*dz
    dr2 = dr*dr

    #set radia
    r = numpy.zeros_like(phi)
    for i in range(nz):
        for j in range(nr):
            r[i][j] = R(j)
    
    for it in range (100):
        t1[1:-1] = (phi[2:]  + phi[:-2]) / dz2
        t2[:,1:-1] = (phi[:,2:]*0.5*(r[:,2:]+r[:,1:-1]) + 
                      phi[:,:-2]*0.5*(r[:,:-2]+r[:,1:-1]))/(dr2*r[:,1:-1])
               
        Ax = t1+t2
        rho_e = QE*n0*numpy.exp(numpy.subtract(P,phi0)/kTe)
        b = -(rho_i - rho_e)/EPS0
   
        g = (b-Ax)
        g[:,1:-1] /= (-2/dz2 - 
                      (r[:,1:-1]+0.5*r[:,0:-2]+0.5*r[:,2:])/(dr2*r[:,1:-1]))
         
        #neumann boundaries
        g[0] = phi[1]
        g[-1] = phi[-2]
        g[:,0] = phi[:,1]
        g[:,-1] = phi[:,-2]
         
        #dirichlet nodes
        phi = numpy.where(cell_type>0,P,g)
    return phi

#computes electric field                    
def computeEF(phi,efz,efr):
    #central difference
    efz[1:-1] = (phi[0:nz-2]-phi[2:nz+1])/(2*dz)
    efr[:,1:-1] = (phi[:,0:nr-2]-phi[:,2:nr+1])/(2*dr)
    
    #one sided difference on boundaries
    efz[0,:] = (phi[0,:]-phi[1,:])/dz
    efz[-1,:] = (phi[-2,:]-phi[-1,:])/dz
    efr[:,0] = (phi[:,0]-phi[:,1])/dr
    efr[:,-1] = (phi[:,-2]-phi[:,-1])/dr
    
    
#---------- INITIALIZATION ----------------------------------------

pl.close('all')
seed()
       
# allocate memory space
nz = 35
nr = 12
dz = 1e-3
dr = 1e-3    
dt = 5e-8

QE = 1.602e-19
AMU =  1.661e-27
EPS0 = 8.854e-12

charge = QE
m = 40*AMU  #argon ions  
qm = charge/m                   
spwt = 2e2

#solver parameters
n0 = 1e13
phi1 = 250
phi0 = 0
kTe = 2

phi = numpy.zeros([nz,nr])
efz = numpy.zeros([nz,nr])
efr = numpy.zeros([nz,nr])
rho_i = numpy.zeros([nz,nr])

# ---- sugarcube domain --------------------
cell_type = numpy.zeros([nz,nr]);
tube1_radius = 0.005;
tube1_length = 0.01;
tube1_aperture_rad = 3*dr;
tube2_radius = tube1_radius+dr;
tube2_length = tube1_length+2*dz;
tube2_aperture_rad = 2*dr;
[tube_i_max, tube_j_max] = map(int, XtoL([tube1_length, tube1_radius]))

for i in range(0,nz):
    for j in range(0,nr):
        pos = Pos([i,j])  # node position
        
        #inner tube
        if ((i==0 and pos[1]<tube1_radius) or 
        (pos[0]<=tube1_length and pos[1]>=tube1_radius and pos[1]<tube1_radius+0.5*dr) or
        (pos[0]>=tube1_length and pos[0]<tube1_length+0.5*dz and 
            pos[1]>=tube1_aperture_rad and pos[1]<tube1_radius) ) :
            cell_type[i][j]=1
            phi[i][j] = phi1
        
        if ((pos[0]<=tube2_length and pos[1]>=tube2_radius and pos[1]<tube2_radius+0.5*dr) or
        (pos[0]>=tube2_length and pos[0]<=tube2_length+0.5*dz and 
            pos[1]>=tube2_aperture_rad and pos[1]<tube2_radius) ) :
            cell_type[i][j]=2
            phi[i][j] = phi0
                       
#create an array of particles
particles = []

lambda_d = math.sqrt(EPS0*kTe/(n0*QE))
print ("Debye length is %.4g, which is %.2g*dz"%(lambda_d,lambda_d/dz))


#solve potential
phi=solvePotential(phi,rho_i)

#positions for plotting
pos_r = numpy.linspace(0,(nr-1)*dr,nr)
pos_z = numpy.linspace(0,(nz-1)*dz,nz)
fig1 = pl.figure(num=None, figsize=(10, 10), dpi=80, facecolor='w', edgecolor='k')
sub = (fig1.add_subplot(111))
        
        
ax=sub
ax.hold(False)
cf =ax.contourf(pos_z, pos_r, numpy.transpose(phi),8,alpha=.75,linewidth=1,cmap='jet')
ax.hold(True);
(ZZ,RR)=pl.meshgrid(pos_z,pos_r)
ax.scatter(ZZ,RR,c=numpy.transpose(cell_type),cmap='jet')
ax.set_yticks(pos_r)
ax.set_xticks(pos_z)
pl.grid(b=True,which='both',color='k',linestyle='-')
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels([])
pl.xlim(min(pos_z),max(pos_z))
pl.ylim(min(pos_r),max(pos_r))
pl.gca().set_aspect('equal', adjustable='box')
pl.colorbar(cf,orientation='horizontal',shrink=0.75, pad=0.01)
pl.draw()
      

#this will block execution until figure is closed
pl.show()

