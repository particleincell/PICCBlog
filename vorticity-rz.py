'''flow through a cavity, RZ version
See below URL for more info
https://www.particleincell.com/2016/vorticity-streamfunction-cylindrical
'''

#enable this if running on shell system without graphics
#import matplotlib
#matplotlib.use('Agg')

import numpy as np
import pylab as pl

#flags nodes as follows:
OPEN = 0
WALL = 1
INLET = -1
OUTLET= -2    
def makeGeometry():    
    node_type = np.zeros((ni,nj))
        
    #left wall
    node_type[0,0:ni] = WALL
    node_type[0,0:inlet_nn] = INLET
    
    #top wall
    node_type[:ni-outlet_dz,nj-1] = WALL
        
    #right wall
    node_type[ni-3-outlet_dz:ni-outlet_dz, outlet_nn:] = WALL
   
    return node_type

#sets boundary conditions on psi
def initPsi():
    psi = np.zeros((ni,nj))
    
    psi[:0] = 0     #streamline along axis of rotations
    
    #volumetric flow = 2*pi*psi
    psi_wall = 0.5*u0*(pos_r[inlet_nn]**2)
    print("Inlet flow rate: %g"%(psi_wall))
    
    psi[node_type==WALL] = psi_wall
            
    return psi
    
#solves d^2psi/dz^2 + d^2psi/dr^2  = -w
def computePsi(w, psi, u, v):              
    psi2 = np.copy(psi)
    
    idz2 = 1/(dz*dz)
    idr2 = 1/(dr*dr)
    for it in range(1001):   
        psi2[1:ni-1,1:nj-1] = (w[1:ni-1,1:nj-1]*pos_r[1:nj-1] + 
                             idz2*(psi[0:ni-2,1:nj-1]+psi[2:ni,1:nj-1])+
                             idr2*(psi[1:ni-1,0:nj-2]+psi[1:ni-1,2:nj])-
                             1/(2*dr*pos_r[1:nj-1])*(psi[1:ni-1,2:nj]-psi[1:ni-1,0:nj-2])
                             )/(2*(idz2+idr2))        
                            
        #replace values on boundary nodes with previous values
        psi2[node_type>0] = psi[node_type>0]
                    
        #inlet, dpsi/dz = -v = 0           
        for j in range(nj):
            if (node_type[0,j]==INLET):
                psi2[0,j] = psi2[1,j] 
        
        #dpsi/dz = -v on zmax 
        for j in range(nj):
            if (node_type[ni-1,j]!=WALL):
                psi2[ni-1,j] = psi2[ni-2,j] - dz*(v[ni-1,j])*pos_r[j]
        
        #dpsi/dr = -u = 0 on rmax
        psi2[ni-outlet_dz:,nj-1] = psi2[ni-outlet_dz:,nj-2]
 
        #copy down solution        
        psi = np.copy(psi2) 
            
        #check for convergence                
        if (it%25==0):
            R = np.zeros_like(psi)
            R[1:ni-1,1:nj-1] = (w[1:ni-1,1:nj-1]*pos_r[1:nj-1] + 
                    idz2*(psi[0:ni-2,1:nj-1]-2*psi[1:ni-1,1:nj-1]+psi[2:ni,1:nj-1])+
                    idr2*(psi[1:ni-1,0:nj-2]-2*psi[1:ni-1,1:nj-1]+psi[1:ni-1,2:nj])-
                    1/(2*dr*pos_r[1:nj-1])*(psi[1:ni-1,2:nj]-psi[1:ni-1,0:nj-2]))
            R[node_type>0] = 0
            norm = np.linalg.norm(R)
            if (norm<1e-8): 
                return psi2
                
    print("computePsi failed to converge, norm = %g"%norm)
    return psi

#---------------- VORTICITY -------------------------------
def applyVorticityBoundaries(w,psi,u,v):
    dz2 = dz*dz
    dr2 = dr*dr    
    #apply boundaries
    for i in range(ni):
        for j in range(nj):
            count = 0
            ww = 0
            
            #left wall
            if (i<ni-1 and node_type[i,j]==WALL and node_type[i+1,j]==OPEN):
                ww += 2*(psi[i,j]-psi[i+1,j])/(pos_r[j]*dz2) - 2*v[i,j]/dz
                count += 1
                
            #right wall
            if (i>0 and node_type[i,j]==WALL and node_type[i-1,j]==OPEN):
                ww += 2*(psi[i,j]-psi[i-1,j])/(pos_r[j]*dz2) + 2*v[i,j]/dz
                count += 1
            
            #top wall
            if (j>0 and node_type[i,j]==WALL and node_type[i,j-1]==OPEN):
                ww += 2*(psi[i,j]-psi[i,j-1])/(pos_r[j]*dr2) - 2*u[i,j]/dr + u[i,j]/pos_r[j]
                count +=1 
            
            #set values
            if (count>0):
                w[i,j] = ww/count
                
    #outlet on right side, dw/dz = 0
    w[ni-1,:] = w[ni-2,:]
    
    #outlet on rmax, dw/dr = 0
    for i in range(ni-outlet_dz,ni):
        w[i,nj-1]=w[i,nj-2]
    
    #inlet on left side
    for j in range(1,inlet_nn):
        if (j<nj-1):
            du_dr = (u[0,j+1]-u[0,j-1])/(2*dr)
        else:
            du_dr = (u[0,j]-u[0,j-1])/dr                    
        #w[i,j] = 2*(psi[i,j]-psi[i+1,j])/dz2 - 2*v[i,j]/dz + du_dr                 
        w[0,j] = -du_dr
    
    #this should already be set, w=0 on axis
    w[:,0] = 0

#computes RHS for vorticity equation
def R(w):
    dz2 = dz*dz
    dr2 = dr*dr    
    #make copy so we use consistent data
    r = np.zeros_like(w)
    for i in range(1,ni-1):
        for j in range(1,nj-1):
            if (node_type[i,j]>0): continue
                
            #viscous term, d^2w/dz^2+d^2w/dr^2+(1/r)dw/dr
            A = nu*(
                (w[i-1][j]-2*w[i][j]+w[i+1][j])/dz2 + 
                (w[i][j-1]-2*w[i][j]+w[i][j+1])/dr2 +
                (w[i][j+1]-w[i][j-1])/(2*dr*pos_r[j]))
            
            #convective term u*dw/dz    
            B = u[i][j]*(w[i+1][j]-w[i-1][j])/(2*dz)
            
            #convective term v*dw/dr
            C = v[i][j]*(w[i][j+1]-w[i][j-1])/(2*dr)
                      
            
            r[i][j] = A - B - C
    return r
            
#advances vorticity equation using RK4
def advanceRK4(w,psi,u,v):
    applyVorticityBoundaries(w,psi,u,v)
        
    #compute the four terms of RK4
    Rk = R(w)
    w1 = w + 0.5*dt*Rk
    
    R1 = R(w1)
    w2 = w + 0.5*dt*R1
    
    R2 = R(w2)
    w3 = w + dt*R2
    
    R3 = R(w3)
    w_new = w + (dt/6.0)*(Rk + 2*R1 + 2*R2 +R3) 
              
    #return new value
    return w_new

#differentiates psi to get velocity
def computeVel(psi):
    u = np.zeros_like(psi)
    v = np.zeros_like(psi)
    for i in range (ni):
        for j in range (1,nj-1):
            
            #skip over walls, otherwise differencing on neighbors will be off
            if (node_type[i,j]==WALL): 
                continue
            #v = -dpsi/dz
            if (i==0):
                v[i,j] = -(psi[i+1,j]-psi[i,j])/(dz*pos_r[j])
            elif (i==ni-1):
                v[i,j] = -(psi[i,j]-psi[i-1,j])/(dz*pos_r[j])
            else:
                v[i,j] = -(psi[i+1,j]-psi[i-1,j])/(2*dr*pos_r[j])
            
            #u = dpsi/dr            
            u[i,j] = (psi[i,j+1] - psi[i,j-1])/(2*dr*pos_r[j])
       
    #u velocity on the axis from q=2*pi*psi
    u[:,0] = 2*psi[:,1]/((pos_r[1])**2)
    
    #similar approach to get u velocity on nj-1
    #first compute u[:,0.5]
    u[:,nj-1] = (psi[:,nj-1] - psi[:,nj-2])/dr
    
    u[node_type==WALL] = 0
    #v=0 on axis
    v[:,0] = 0

    return (u,v)
                        
# animation function
def make_plot(it): 
    pl.figure(fig1.number)
    fig1.clear()
    img=pl.contourf(pos_z,pos_r,np.transpose(w))
    pl.colorbar(img)
    pl.axis('equal')
    pl.xlim((pos_z[0],pos_z[ni-1]))
    pl.ylim((pos_r[0],pos_r[nj-1]))
    pl.title("w, time = %d"%it)    

    pl.figure(fig2.number)
    fig2.clear()
    img=pl.contourf(pos_z,pos_r,
                np.transpose(psi)
                )
    pl.colorbar(img)
    (Z,R)=pl.meshgrid(pos_z,pos_r)
    pl.scatter(Z,R,c=np.transpose(node_type),cmap='jet')
    pl.axis('equal')
    pl.xlim((pos_z[0],pos_z[ni-1]))
    pl.ylim((pos_r[0],pos_r[nj-1]))
    pl.title("psi, time = %d"%it)    

    pl.figure(fig3.number)
    fig3.clear()
    fig3.hold(False)
    img=pl.contourf(pos_z,pos_r,
                np.transpose(np.sqrt(u*u+v*v))
#                np.transpose(u)
#                ,levels=np.linspace(0,0.30,num=10)
                )
    pl.colorbar(img)
    pl.hold(True)
    pl.quiver(pos_z,pos_r, np.transpose(u), np.transpose(v))  
    pl.axis('equal')
    pl.xlim((pos_z[0],pos_z[ni-1]))
    pl.ylim((pos_r[0],pos_r[nj-1]))
    pl.title("|vel|, time = %d"%it)    
    pl.savefig("plots/vel%05d"%it)
    
    #also compute average u velocity
    u_ave = np.zeros(ni)
    flux = np.zeros(ni)
    for i in range (ni):
        A = 0       #cross-sectional area
        for j in range (nj):
            if (node_type[i,j]==WALL):
                continue
            
            if (j==0):
                r1 = pos_r[j]                
                r2 = r1+0.5*dr
            elif (j==nj-1):
                r2 = pos_r[j]
                r1 = r2-0.5*dr
            else:
                r1 = pos_r[j]-0.5*dr
                #if (node_type[i,j+1]==WALL):
                #    r2 = pos_r[j+1]
                #else:
                r2 = pos_r[j]+0.5*dr
                
            dA = np.pi*(r2**2-r1**2)
            A = A + dA
            
            flux[i] = flux[i] + dA*u[i,j]
        u_ave[i] = flux[i]/A
        
    pl.figure(fig4.number)
    fig4.clear()
    fig4.hold(False)                    
    pl.plot(pos_z,flux,'k',linewidth=2)
    pl.title("flux, time = %d"%it)
    print("It %d"%it)
    pl.pause(0.0001)
    return flux,u_ave
          
#main program     
ni = 31
nj = 15

dt=2.0e-3
dz=0.0025
dr=0.0025

#parameters
u0 = 0.1   #inlet velocity
nu0 = 1.568e-5   #air kinematic viscosity at 300K
nu_k = 1        #artifical visocisity factor
nu = nu_k*nu0
inlet_nn = 3  #number of nodes making up inlet
outlet_nn = 3 #number of nodes in the outlet
outlet_dz = 6    #number of nodes outside the cavity
pos_z = np.arange(0,ni)*dz
pos_r = np.arange(0,nj)*dr

#generate geometry
node_type = makeGeometry()
    
#screen output
d=nu*dt*(1/(dz*dz)+1/(dr*dr))
print("d=%g, should be <= 0.5"%d)
u_max=u0*(inlet_nn/outlet_nn) #mass conservation through outlet
Re_x = u_max*dz/nu
cx=u0*dt/dz
Re_main = u0*pos_r[nj-1]*2/nu
Re_outlet = u_max*pos_r[outlet_nn-1]*2/nu
print("%g <= %g <= %g"%(2*cx,Re_x,2/cx))
print("u0 = %g, u_max = %g"%(u0,u_max))
print("Re_main = %g, Re_outlet = %g"%(Re_main,Re_outlet))
#set streamfunction boundary conditions
psi = initPsi()

#set initial values to zeros
w = np.zeros_like(psi)
u = np.zeros_like(psi)
v = np.zeros_like(psi)

#update u and v
it = 0
#pl.ioff()
fig1=pl.figure(1,figsize=(10,4))
fig2=pl.figure(2,figsize=(10,4))
fig3=pl.figure(3,figsize=(10,4))
fig4=pl.figure(4,figsize=(8,4))
pl.ion()

#iterate
print ("Starting main loop")
for it in range (1001):
       
    #solve psi
    psi = computePsi(w,psi,u,v)
  
    #update u and v
    u,v = computeVel(psi)
    
    #advance w
    w = advanceRK4(w,psi,u,v)
        
    if (it%100==0):
        flux,u_ave = make_plot(it)

print ("Done!")