#Demo 1D sheath simulation

import numpy
import pylab as pl
import math
from random import (seed,random)


# constants
EPS_0 = 8.85418782e-12      # F/m, vacuum permittivity
K =    1.38065e-23        # J/K, Boltzmann constant
ME = 9.10938215e-31        # kg, electron mass
QE = 1.602176565e-19    # C, elementary charge
AMU  = 1.660538921e-27    # kg, atomic mass unit
EV_TO_K =    11604.52        # 1eV in Kelvin, QE/K

#simulation parameters, these could come from an input file
PLASMA_DEN = 1e16        # plasma density to load
NUM_IONS = 100000        # number of ions
NUM_ELECTRONS = 100000    # number of electrons
dx = 1e-4                # cell spacing
ni = 101                # number of nodes
NUM_TS =    1000            # number of time steps
DT = 1e-11            # time step size
ELECTRON_TEMP = 3.0        # electron temperature in eV
ION_TEMP = 1.0            # ion temperature in eV

def XtoL(x):
    lc = pos/dx
    return lc

def Pos(lc):
    pos = x0 + lc*dx
    return pos
  
def gather(data,lc):
    i = math.trunc(lc)
    di = lc - i
    return  (data[i]*(1-di) +
             data[i+1]*(di))
    
def scatter(data,lc,value):
    i = numpy.trunc(lc)
    di = lc - i
            
    data[i] += (1-di)*value
    data[i+1] += (di)*value
       
#particle definition
class Particle:
    def __init__(self,pos,vel):
        self.x=pos 
        self.v=vel

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
    vm[0:3] = math.sqrt(2)*vth*(2*(random()+random()+random()-1.5))
    
    vel = (n[0]*vm[0], n[1]*vm[1], n[2]*vm[2]) 
    return vel

#Thomas algorithm for a tri-diagonal matrix*/
def solvePotential(x,max_it=100):
    #set coefficients, this should be pre-computed
    dx2 = dx*dx
    a = np.zeros(ni)
    b = np.zeros(ni)
    c = np.zeros(ni)
    
    #central difference on internal nodes
    a[1:-1] = 1
    b[1:-1] = -2
    c[1:-1] = 1
         
    #dirichlet b.c. on boundaries
    a[0]=0
    b[0]=1
    c[0]=0
    a[ni-1]=0
    b[ni-1]=1
    c[ni-1]=0
         
    #multiply RHS
    x = -rho*dx2/EPS_0
    x[0] = 0
    x[ni-1] = 0
    
    #Modify the coefficients, b[0] implies singular matrix/
    c[0] /= b[0]
    x[0] /= b[0]                
    
    for i in range(1,ni):
         d = (b[i] - c[i-1] * a[i])    
         c[i] /= d                            
         x[i] = (x[i] - x[i-1] * a[i])/d
    
    #Now back substitute.
    for i in range (ni-2,0,-1):
        x[i] = x[i] - c[i] * x[i + 1]

#computes electric field                    
def computeEF(phi,efx):
    
    #central difference, not right on walls
    efx[1:-1] = (phi[0:]-phi[2:])/(2*dx)
    
    #one sided difference on boundaries
    efz[0] = (phi[0]-phi[1])/dx
    efz[-1] = (phi[-2]-phi[-1])/dx
    
def plot(ax,data,scatter=False):
    pl.sca(ax)
    pl.cla()
    ax.set_xticks(pos_x)
    ax.xaxis.set_ticklabels([])
    pl.xlim(min(pos_x),max(pos_x))
    ax.grid(b=True,which='both',color='k',linestyle='-')
    ax.set_aspect('equal', adjustable='box')
 #   pl.colorbar(cf,ax=pl.gca(),orientation='horizontal',shrink=0.75, pad=0.01)

   
#samples random velocity from Maxwellian distribution using Birdsall's method
def SampleVel(v_th):
    M = 12;    
    sum = 0
    for i in range(M):
        sum+=random()

    return sqrt(0.5)*v_th*(sum-M/2.0)/sqrt(M/12.0)


#scatter particles of species to the mesh
def ScatterSpecies(particles, den):
    #initialize densities to zero
    den[:] = 0
    
    #scatter particles to the mesh
    for p in range(p):
        lc = XtoL(part[p,0]);
        scatter(lc, spwt, den);
    
    #divide by cell volume
    den[:]/=dx;
        
    #only half cell at boundaries
    den[0] *=2.0;
    den[ni-1] *= 2.0;

#adds new particle to the species, returns pointer to the newly added data
def AddParticle(species, x, v):
    #store position and velocity of this particle
    species.part[species.np,0] = x
    species.part[species.np,1] = v
    
    #increment particle counter
    species.np+=1

#computes charge density by adding ion and electron data
def ComputeRho(ni, ne):
    rho[:] = QE*(ni[:]-ne[:])

#moves particles of a single species, returns wall charge
def PushParticles(species, ef):
    #precompute q/m
    qm = species.charge / species.mass
    
    #loop over particles
    for p in range(species.np):
        #compute particle node position
        lc = XtoL(particles[p,0])
        
        #gather electric field onto particle position
        part_ef = gather(lc,ef)
        
        #advance velocity/
        particles[p,1] += DT*qm*part_ef    
        
        #advance position
        particles[p,0] += DT*particles[p,1]
            
        #remove particles leaving the domain
        if (particles[p,0] < x0 or particles[p,0] >= xmax):
            #replace this particle slot with last one*/
            particles[p,:] = particles[p,species.np-1,:]
            species.np -=1  #reduce number of particles
            p-=1    #decrement to reprocess this slot
        
# rewinds particle velocities by -0.5DT*/
def rewindParticle(particles, ef):
    #precompute q/m
    qm = species.charge / species.mass
    
    for p in range(len(particles)):

        #compute particle node position
        lc = XtoL(part[p,0]);
        
        #gather electric field onto particle position
        part_ef = gather(lc,ef)
        
        #advance velocity
        part[p,1] -= 0.5*DT*qm*part_ef            
    
#--------- main code
    

#---1) initialize domain---*/
    
phi =  np.zeos(ni)
    
    
#set material data*/
ions.mass = 16*AMU        
ions.charge = QE
ions.spwt = PLASMA_DEN*domain.xl/NUM_IONS
ions.np = 0
ions.np_alloc = NUM_IONS
ions.part = new Particle[NUM_IONS]
    
electrons.mass = ME;    // electrons
electrons.charge = -QE;
electrons.spwt = PLASMA_DEN*domain.xl/NUM_ELECTRONS;
electrons.np = 0;
electrons.np_alloc = NUM_ELECTRONS;
electrons.part = new Particle[NUM_ELECTRONS];

ions = np.zeros((2,NUM_IONS))
electrons = np.zeros((2,NUM_ELECTRONS))
    
#*randomize RNG*/
seed()      #reset random number generator
    
#/*load uniformly spaced ions and electrons*/
delta_ions = domain.xl/NUM_IONS
v_thi = sqrt(2*K*ION_TEMP*EV_TO_K/ions.mass)    

for p in range(NUM_IONS):
    x = x0 + p*delta_ions
    v = SampleVel(v_thi)
    AddParticle(ions,x,v)
        
#now do the same for electrons*/
delta_electrons = xl/NUM_ELECTRONS
v_the = sqrt(2*K*ELECTRON_TEMP*EV_TO_K/electrons.mass)    
for p in range(NUM_ELECTRONS):
    x = x0 + p*delta_electrons;
    v = SampleVel(v_the);
    AddParticle(&electrons,x,v);
            
    #/*compute number density*/
    ScatterSpecies(&ions,ndi);
    ScatterSpecies(&electrons,nde);
        
    #/*compute charge density and solve potential*/
    ComputeRho(&ions, &electrons);
    SolvePotential(phi,rho);
    ComputeEF(phi,ef);
    
     #
    RewindSpecies(&ions,ef);
    RewindSpecies(&electrons,ef);
    
    #/*OUTPUT*/
        
    clock_t start = clock();    // grab starting clock time

    #/* MAIN LOOP*/
    for (ts = 1; ts<=NUM_TS; ts++)
    {
        /*compute number density*/
        ScatterSpecies(&ions,ndi);
        ScatterSpecies(&electrons,nde);
        
        ComputeRho(&ions, &electrons);
        SolvePotential(phi,rho);
        ComputeEF(phi,ef);
        
        /*move particles*/
        PushSpecies(&electrons,ef);
        PushSpecies(&ions,ef);
        
        /*write diagnostics*/
        if (ts%25==0)
        {
            /*max phi*/
            double max_phi = phi[0];
            for (int i=0;i<domain.ni;i++)
                if (phi[i]>max_phi) max_phi=phi[i];
            
            printf("TS:%i\tnp_i:%d\tnp_e:%d\tdphi:%.3g\n", 
                ts,ions.np,electrons.np,max_phi-phi[0]);            
        }
        
        /*save data*/
        if (ts%1000==0) WriteResults(ts);
    }
    
    clock_t end=clock();
    
    fclose(file_res);
    
    
    printf("Time per time step: %.3g ms\n",1000*(end-start)/(double)(CLOCKS_PER_SEC*NUM_TS));

    return 0;

#***** HELPER FUNCTIONS *********************************************************/
