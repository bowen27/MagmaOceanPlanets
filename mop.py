## Goal: simulate a 1D boundary layer of a magma ocean planet

# Pending issues: 
# (1) what are the boundary conditions associated with each variable?
# (2) what are the initial conditions associated with each variable?
# (3) validation and verification of each function
# (4) safeguards that prevent negative quantities, etc.

# Additional issues:
# (1) Since we use upwind scheme, evaluating u[t+1,j=2] in get_u() relies on u[t,j=1],
# which is zero because of the boundary condition. Therefore, the gradient could be large.
# (2) Because of upwind scheme, we have do decide what to do at j=0 for get_e. 
# (3) Because of upwind scheme, we have do decide what to do at j=0 for get_xidot. 

# Import python libraries
import numpy as np
import matplotlib.pyplot as plt 
import math

# Create parameter class

class parameters:
    # Spatial grid
    rp   = 6314e3 
    ndeg = 180    
    ddeg = 1
    dx   = rp*np.pi/ndeg
    x    = np.arange(0,ndeg+ddeg,ddeg)*dx
    xlen = len(x)
    # Spatial indices of the interior
    jmin   = 1
    jmax   = xlen-2
    jmaxp1 = jmax+1
    # Spatial indices of left boundary and right boundary
    jlb   = 0
    jrb   = xlen-1
    # Temporal grid
    dt   = 0.2*dx
    tmin = 0
    tmax = dt*5
    t    = np.arange(tmin,tmax+dt,dt)
    tlen = len(t)
    # CFL criterion must be met
    if not dt/dx < 0.4:
        print('CFL criterion is not met.')
        quit()
    # Thermodynamics **
    M  = 18e-3    # kg/mol
    R  = 8.314/M  # J/kg/K
    cp = 1864     # J/kg/K
    cv = cp-R     # J/kg/K
    L  = 2260e3   # J/kg
    # Clausius-Clapeyron **
    pref = 3533 # Pa
    Tref = 300  # K
    # Free Atmosphere
    p0  = 2e4    # Pa
    g   = 10     # m/s^2
    tau = dt*100 # damping timescale (s)
    # Surface Energy Imbalance (+- 5 W/m^2)
    Fnet = np.cos(np.pi*np.arange(0,ndeg+ddeg,ddeg)/ndeg)*5
    
    def __init__(self):
        self.Ts , self.ps   , self.rhodel, self.xidotdel,\
        self.u  , self.e    , self.T     , self.p       ,\
        self.rho, self.delta, self.xidot , self.Fnet    ,\
        self.Cm , self.Cu   , self.Ce    , self.ubar     \
        = (np.zeros([parameters.tlen,parameters.xlen],dtype='float') for i in range(16))
                
# Saturation Vapor Pressure
def get_esat(x):
    # Input x could be T or Ts
    esat = np.zeros([par.xlen],dtype='float')
    esat = par.pref*np.exp(-(par.L/par.R)*(1/x - 1/par.Tref))
    return esat

def get_Cm():
    # Diagnostic at current time step, i
    Cm = np.zeros([par.xlen],dtype='float')
    Cm = par.xidot[i,:]*par.delta[i,:]
    return Cm

def get_Cu():
    # Diagnostic at current time step, i
    Cu = np.zeros([par.xlen],dtype='float')
    for j in range(par.xlen):
        if par.xidot[i,j]>=0:
            Cu[j] = -par.rho[i,j]*par.delta[i,j]*par.u[i,j]/par.tau
            
        else:
            Cu[j] = -par.rho[i,j]  *par.delta[i,j]*par.u[i,j]/par.tau + \
                     par.xidot[i,j]*par.delta[i,j]*par.u[i,j]
    return Cu       
    
def get_Ce():
    # Diagnostic at current time step, i
    Ce = np.zeros([par.xlen],dtype='float')
    for j in range(par.xlen):
        if par.xidot[i,j]>=0:
            Ce[j] = par.xidot[i,j]*par.delta[i,j]*par.cv*par.Ts[i,j]
        else:
            Ce[j] = par.xidot[i,j]*par.delta[i,j]*par.e[i,j]     
    return Ce  

def get_Ts(): 
    # Prognostic
    Ts = np.zeros([par.xlen],dtype='float')
    Ts = par.Ts[i,:] + par.dt/par.cp*\
            (\
             par.Fnet[i,:]-par.L*par.xidot[i,:]*par.delta[i,:]\
            )     
    return Ts

def get_ps(): 
    # Diagnostic
    ps = np.zeros([par.xlen],dtype='float')
    ps = get_esat(par.Ts[i+1,:])
    return ps

def get_rhodel(): 
    # Diagnostic
    rhodel = np.zeros([par.xlen],dtype='float')
    # interior
    rhodel = (par.ps[i+1,par.jmin:par.jmaxp1]-par.p0)/par.g
    # left boundary
    rhodel[par.jlb] = rhodel[par.jbl+1]
    # right boundary
    rhodel[par.jrb] = rhodel[par.jrb-1]
    return rhodel

def get_xidotdel(): # ** are we solving for Cm(t+1) or Cm(t)? Does it make sense to solve for Cm(t+1)?
    # Prognostic 
    xidotdel = np.zeros([par.xlen],dtype='float')
    for j in range(par.xlen):
        if j!=0:
            xidotdel[j] = (par.rhodel[i+1,:]-par.rhodel[i,:])/par.dt + \
                          (par.u[i,j]*par.rhodel[i,j]-par.u[i,j-1]*par.rhodel[i,j-1])/par.dx
        else: 
            # boundary condition at left wall
            xidotdel[j] = 0
    return xidotdel

def get_u(): 
    # Prognostic
    rhodelu = np.zeros([par.xlen],dtype='float') # (t+1)
    u = np.zeros([par.xlen],dtype='float')       # (t+1)
    for j in range(par.xlen):
        if j!=0:
            rhodelu[j] = par.rho[i,j]*par.delta[i,j]*par.u[i,j] +                                \
                         par.dt*(par.Cu[i,j] -                                                   
                                (                                                                
                                par.delta[i,j]  *(par.rho[i,j]  *par.u[i,j]**2   + par.p[i,j])-   
                                par.delta[i,j-1]*(par.rho[i,j-1]*par.u[i,j-1]**2 + par.p[i,j-1]) 
                                )/                                                               
                                (par.x[j]-par.x[j-1])                                            
                                )
        else: 
            # boundary condition at left wall
            rhodelu[j]=0 # **
    u = rhodelu[:]/par.rhodel[i+1,:]
    return u

def get_e(): # **
    # Prognostic
    rhodele = np.zeros([par.xlen],dtype='float') # (t+1)
    e = np.zeros([par.xlen],dtype='float')       # (t+1)
    for j in range(par.xlen):
        if j!=0:
            rhodele[j] = par.rho[i,j]*par.delta[i,j]*par.e[i,j] +                                                  \
                         par.dt*(par.Ce[i,j] -                                                                     
                                    (                                                                             
                                    par.delta[i,j]  *par.u[i,j]  *(par.rho[i,j]  *par.e[i,j]   + par.p[i,j])-     
                                    par.delta[i,j-1]*par.u[i,j-1]*(par.rho[i,j-1]*par.e[i,j-1] + par.p[i,j-1])    
                                    )/                                                                            
                                    (par.x[j]-par.x[j-1])                                                         
                                )
        else: 
            # boundary condition at left wall
            rhodele[j]=0 #**
    e = rhodele[:]/par.rhodel[i+1,:]
    return e

def get_T():
    # Diagnostic
    T = np.zeros([par.xlen],dtype='float') # (t+1)
    # interior
    T[par.jmin:par.jmaxp1] = (1/par.cv)*(par.e[i+1,par.jmin:par.jmaxp1]-0.5*par.u[i+1,par.jmin:par.jmaxp1]**2)
    # left boundary
    T[par.jlb] = T[par.jlb+1]
    # right boundary
    T[par.jrb] = T[par.jrb-1]
    return T

def get_p():
    # Diagnostic
    p = np.zeros([par.xlen],dtype='float') # (t+1)
    p = get_esat(par.T[i+1,:])
    return p

def get_rho():
    # Diagnostic
    rho = np.zeros([par.xlen],dtype='float') # (t+1)
    rho = par.p[i+1,:]/(par.R*par.T[i+1,:])
    return rho

def get_delta():
    # Diagnostic
    delta = np.zeros([par.xlen],dtype='float') # (t+1)
    delta = par.rhodel[i+1,:]/par.rho[i+1,:]
    return delta
    
def get_xidot(): # **
    # Diagnostic
    xidot = np.zeros([par.xlen],dtype='float') # (t+1)
    for j in range(par.xlen):
        xidot[j] =  (
                    par.Fnet[i+1,j]*(par.ps[i+1,j]*par.L)/(par.R*par.Ts[i+1,j]**2*par.cp)  
                    + par.g*(                                                                 
                            par.rho[i+1,j]  *par.delta[i+1,j]  *par.u[i+1,j]    -              
                            par.rho[i+1,j-1]*par.delta[i+1,j-1]*par.u[i+1,j-1]        
                            )/(par.x[j]-par.x[j-1])
                    )/\
                    (                                                                       
                   par.delta[i+1,j]*(par.g + 
                                    (par.ps[i+1,j]*par.L**2)/(par.R*par.Ts[i+1,j]**2*par.cp)
                                    )
                   )    
    return xidot

def get_initial_conditions():
    # Scenario: intially saturated atmosphere with uniform Ts,T,p,delta,etc.
    par.Ts[0,:]       = np.ones([par.xlen],dtype='float')*350 # uniform initial temp
    par.ps[0,:]       = np.ones([par.xlen],dtype='float')*4e4 # uniform initial pressure
    par.rhodel[0,:]   = (par.ps[0,:]-par.p0)/par.g            # uniform atmospheric mass
    par.u[0,:]        = np.zeros([par.xlen],dtype='float')    # no wind
    par.T[0,:]        = par.Ts[0,:]                           # T = Ts
    par.e[0,:]        = par.cv*par.T[0,:] + 0.5*par.u[0,:]**2 
    par.p[0,:]        = par.ps[0,:]                           # p = ps
    par.rho[0,:]      = par.p[0,:]/(par.R*par.T[0,:])         
    par.delta[0,:]    = par.rhodel[0,:]/par.rho[0,:]          
    par.xidot[i+1,:]  = np.zeros([par.xlen],dtype='float')    # zero surface mass-flux
    
# Time Integration
par = parameters()
get_initial_conditions()

# Plot forward integration
plt.clf()
fig, ax = plt.subplots(7)    

for i in range(par.tlen):
    # Get variables at current time step
    par.Cm[i,:]         = get_Cm()
    par.Cu[i,:]         = get_Cu()
    par.Ce[i,:]         = get_Ce()
    # Get variables at the next time step
    par.Ts[i+1,:]       = get_Ts()       
    par.ps[i+1,:]       = get_ps()       
    par.rhodel[i+1,:]   = get_rhodel()   
    par.u[i+1,:]        = get_u()        
    par.e[i+1,:]        = get_e()        
    par.T[i+1,:]        = get_T()        
    par.p[i+1,:]        = get_p()        
    par.rho[i+1,:]      = get_rho()      
    par.delta[i+1,:]    = get_delta()   
    par.xidot[i+1,:]    = get_xidot()  
    
    ax[0] = ax.plot(par.x/np.rp,par.T[i,:]) 
    ax[1] = ax.plot(par.x/np.rp,par.Ts[i,:]) 
    ax[2] = ax.plot(par.x/np.rp,par.p[i,:]) 
    ax[3] = ax.plot(par.x/np.rp,par.ps[i,:]) 
    ax[4] = ax.plot(par.x/np.rp,par.u[i,:])
    ax[5] = ax.plot(par.x/np.rp,par.delta[i,:])
    ax[6] = ax.plot(par.x/np.rp,par.xidot[i,:])
    
ax[0].set_xlabel('T')
ax[1].set_xlabel('Ts')
ax[2].set_xlabel('p')
ax[3].set_xlabel('ps')
ax[4].set_xlabel('u')
ax[5].set_xlabel('delta')
ax[6].set_xlabel('xidot')

plt.show()
    
    
