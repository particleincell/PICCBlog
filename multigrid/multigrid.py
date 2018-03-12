# -*- coding: utf-8 -*-
"""
Example 1D Multigrid solver
Written by Lubos Brieda for 
https://www.particleincell.com/2018/multigrid-solver/

Fri Mar  9 21:47:15 2018
"""

import numpy as np
import pylab as pl

#standard GS+SOR solver
def GSsolve(phi,b):
    print("**** GS SOLVE ****")

    #plot for solution vector
    fig_sol = pl.figure()
    #plot analytical solution
    pl.plot(x,phi_true,LineWidth=4,color='red',LineStyle="dashed")
    pl.title('GS Solver: phi vs. x')
    
    #plot for error
    fig_err = pl.figure()
    pl.title('GS Solver: (phi-phi_true) vs. x')
     
    R_freq = 100 #how often should residue be computed
    w = 1.4 #SOR acceleration factor
    
    for it in range(100001):
        phi[0] = phi[1] #neumann BC on x=0
        for i in range (1,ni-1):    #Dirichlet BC at last node so skipping
            g = 0.5*(phi[i-1] + phi[i+1] - dx*dx*b[i])            
            phi[i] = phi[i] + w*(g-phi[i])
        
        #compute residue to check for convergence
        if (it%R_freq==0):    
            r_i = phi[1]-phi[0]
            r_sum = r_i*r_i
            for i in range (1, ni-1):
                r_i = (phi[i-1]-2*phi[i]+phi[i+1])/(dx*dx) - b[i]
                r_sum += r_i*r_i
        
            norm = np.sqrt(r_sum)/ni
            if (norm<1e-4): 
                print("Converged after %d iterations"%it)
                op_count = it*ni*5 + (it/R_freq)*5*ni
                print("Operation count: %d"%op_count)
                pl.figure(fig_sol.number)
                pl.plot(x,phi)
                return op_count
                break
    
        if (it % 250 ==0):
            pl.figure(fig_sol.number)
            pl.plot(x,phi)
            pl.figure(fig_err.number)
            pl.plot(x,phi-phi_true)
            #print("it: %d, norm: %g"%(it, norm))
    return -1

#multigrid solver
def MGsolve(phi,b):
    global eps_f, eps_c, R_f, R_c
    print("**** MULTIGRID SOLVE ****")
    ni = len(phi)
    ni_c = ni>>1    #divide by 2
    dx_c = 2*dx
        
    R_f = np.zeros(ni)
    R_c = np.zeros(ni_c)
    eps_c = np.zeros(ni_c)
    eps_f = np.zeros(ni)
 
    pl.figure()
    #plot analytical solution
    pl.plot(x,phi_true,LineWidth=4,color='red',LineStyle="dashed")
    pl.title('Multigrid Solver: phi vs. x')
   
    for it in range(10001):
        
        #number of steps to iterate at the finest level
        inner_its = 1

        #number of steps to iterate at the coarse level        
        inner2_its = 50
        w = 1.4
        
        #1) perform one or more iterations on fine mesh
        for n in range(inner_its):
            phi[0] = phi[1]            
            for i in range (1,ni-1):
                g = 0.5*(phi[i-1] + phi[i+1] - dx*dx*b[i])
                phi[i] = phi[i] + w*(g-phi[i])
                            
        #2) compute residue on the fine mesh, R = A*phi - b
        for i in range(1,ni-1):
            R_f[i] = (phi[i-1]-2*phi[i]+phi[i+1])/(dx*dx) - b[i]
        R_f[0] = (phi[0] - phi[1])/dx    #neumann boundary
        R_f[-1] = phi[-1] - 0   #dirichlet boundary
                    
        #2b) check for termination
        r_sum = 0
        for i in range(ni):
            r_sum += R_f[i]*R_f[i]
        norm = np.sqrt(r_sum)/ni
        if (norm<1e-4): 
            print("Converged after %d iterations"%it)
            print("This is %d solution iterations on the fine mesh"%(it*inner_its))
            op_count_single = (inner_its*ni*5 + ni*5 + (ni>>1) + 
                               inner2_its*(ni>>1)*5 + ni + ni)
            print("Operation count: %d"%(op_count_single*it))
            pl.plot(x,phi)
            return op_count_single*it
            break
        
        #3) restrict residue to the coarse mesh
        for i in range(2,ni-1,2):
            R_c[i>>1] = 0.25*(R_f[i-1] + 2*R_f[i] + R_f[i+1])            
        R_c[0] = R_f[0]
            
        #4) perform few iteration of the correction vector on the coarse mesh
        eps_c[:] = 0
        for n in range(inner2_its):
            eps_c[0] = eps_c[1] + dx_c*R_c[0]
            for i in range(1,ni_c-1):
                g = 0.5*(eps_c[i-1] + eps_c[i+1] - dx_c*dx_c*R_c[i])
                eps_c[i] = eps_c[i] + w*(g-eps_c[i])
            
        #5) interpolate eps to fine mesh
        for i in range(1,ni-1):
            if (i%2==0):    #even nodes, overlapping coarse mesh
                eps_f[i] = eps_c[i>>1]
            else:
                eps_f[i] = 0.5*(eps_c[i>>1]+eps_c[(i>>1)+1])
            eps_f[0] = eps_c[0]
        
        #6) update solution on the fine mesh
        for i in range(0,ni-1):
            phi[i] = phi[i] - eps_f[i]
              
        if (it % int(25/inner_its) == 0):
            pl.plot(x,phi)
#            print("it: %d, norm: %g"%(it, norm))
            
    return -1

ni = 128

#domain length and mesh spacing
L = 1
dx = L/(ni-1)

x = np.arange(ni)*dx

phi = np.ones(ni)*0
phi[0] = 0
phi[ni-1] = 0
phi_bk = np.zeros_like(phi)
phi_bk[:] = phi[:]


#set RHS
A = 10
k = 4
b = A*np.sin(k*2*np.pi*x/L) 

#analytical solution
C1 = A/(k*2*np.pi/L)
C2 = - C1*L
phi_true = (-A*np.sin(k*2*np.pi*x/L)*(1/(k*2*np.pi/L))**2 + C1*x + C2)

pl.figure(1)
pl.title('Right Hand Side')
pl.plot(x,b)

op_gs = GSsolve(phi,b)

#reset solution
phi = phi_bk
op_mg = MGsolve(phi,b)

print("Operation count ratio (GS/MG): %g"%(op_gs/op_mg))


        