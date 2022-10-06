## Goal: simulate a 1D boundary layer of a magma ocean planet

# Pending issues: 
# (1) what are the boundary conditions associated with each variable?
# (2) what are the initial conditions associated with each variable?
# (3) validation and verification of each function
# (4) safeguards that prevent negative quantities, etc.
# (5) Resolve the noise originating at the boundaries

# Import python libraries
import numpy as np
import matplotlib.pyplot as plt 

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
    jlrb = [jlb,jlb+1,jrb-1,jrb]
    # Temporal grid
    dt   = 180 # seconds
    tmin = 0
    tmax = dt*1
    t    = np.arange(tmin,tmax+dt,dt)
    tlen = len(t)
    # CFL criterion must be met
    if not dt/dx < 0.4:
        print('CFL criterion is not met.')
        quit()
    # Thermodynamics of Atmosphere**
    M  = 18e-3    # kg/mol
    R  = 8.314/M  # J/kg/K
    cp = 1864     # J/kg/K
    cv = cp-R     # J/kg/K
    L  = 2260e3   # J/kg
    # Thermodynamics of Ocean
    cpo  = 4184 # J/kg/K
    rhoo = 1000 # kg/m^3
    ho   = 4000 # m
    Co  = cpo*rhoo*ho # J/m^2/K
    # Clausius-Clapeyron **
    pref = 3533 # Pa
    Tref = 300  # K
    # Free Atmosphere
    p0  = 2e4    # Pa
    g   = 10     # m/s^2
    tau = dt*100 # damping timescale (s)

    #plt.clf()
    #fig,ax=plt.subplots()
    #ax.plot(x/rp,Fnet)
    #ax.set_ylabel('Fnet')
    #plt.show()
    
    def __init__(self):
        self.Ts , self.ps   , self.rhodel, self.xidotdel,\
        self.u  , self.e    , self.T     , self.p       ,\
        self.rho, self.delta, self.xidot , self.Fnet    ,\
        self.Cm , self.Cu   , self.Ce    , self.ubar     \
        = (np.zeros([parameters.tlen,parameters.xlen],dtype='float') for i in range(16))

# Surface Energy Imbalance
def get_Fnet():
    # Specified at current time step, i
    Fnet = np.zeros([par.xlen],dtype='float')
    Fnet = np.cos(np.pi*np.arange(0,par.ndeg+par.ddeg,par.ddeg)/par.ndeg)*5 
    return Fnet

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
            Cu[j] = -par.rho[i,j]  *par.delta[i,j]*par.u[i,j]/par.tau \
                    +par.xidot[i,j]*par.delta[i,j]*par.u[i,j]
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
    Ts = par.Ts[i,:] + par.dt/par.Co*\
            (
             par.Fnet[i,:]-par.L*par.xidot[i,:]*par.delta[i,:]
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
    rhodel[par.jmin:par.jmaxp1] = (par.ps[i+1,par.jmin:par.jmaxp1]-par.p0)/par.g
    # left boundary
    rhodel[par.jlb] = rhodel[par.jlb+1]
    # right boundary
    rhodel[par.jrb] = rhodel[par.jrb-1]
    return rhodel

def get_u(): 
    # Prognostic
    rhodelu = np.zeros([par.xlen],dtype='float') # (t+1)
    u = np.zeros([par.xlen],dtype='float')       # (t+1)
    # interior
    for j in range(par.jmin+1,par.jmax):
        if par.u[i,j]>=0:
            rhodelu[j] = par.rho[i,j]*par.delta[i,j]*par.u[i,j] +                                \
                         par.dt*(par.Cu[i,j] -                                                   
                                    (                                                                
                                    par.delta[i,j]  *(par.rho[i,j]  *par.u[i,j]**2   + par.p[i,j])-   
                                    par.delta[i,j-1]*(par.rho[i,j-1]*par.u[i,j-1]**2 + par.p[i,j-1]) 
                                    )/                                                               
                                    (par.x[j]-par.x[j-1])                                            
                                    )
        elif par.u[i,j]<0:
            rhodelu[j] = par.rho[i,j]*par.delta[i,j]*par.u[i,j] +                                \
                         par.dt*(par.Cu[i,j] -                                                   
                                (                                                                
                                par.delta[i,j+1]*(par.rho[i,j+1]*par.u[i,j+1]**2   + par.p[i,j+1])-   
                                par.delta[i,j]  *(par.rho[i,j]  *par.u[i,j]**2     + par.p[i,j]) 
                                )/                                                               
                                (par.x[j+1]-par.x[j])                                            
                                )
    # left and right boundaries are zero by construction
    u = rhodelu[:]/par.rhodel[i+1,:]
    return u

def get_e(): # **
    # Prognostic
    rhodele = np.zeros([par.xlen],dtype='float') # (t+1)
    e = np.zeros([par.xlen],dtype='float')       # (t+1)
    # interior
    for j in range(par.jmin,par.jmaxp1):
        # insert criterion for direction of mean wind
        if par.u[i,j]>=0:
            rhodele[j] = par.rho[i,j]*par.delta[i,j]*par.e[i,j] +                                                  \
                         par.dt*(par.Ce[i,j] -                                                                     
                                    (                                                                             
                                    par.delta[i,j]  *par.u[i,j]  *(par.rho[i,j]  *par.e[i,j]   + par.p[i,j])-     
                                    par.delta[i,j-1]*par.u[i,j-1]*(par.rho[i,j-1]*par.e[i,j-1] + par.p[i,j-1])    
                                    )/                                                                            
                                    (par.x[j]-par.x[j-1])                                                         
                                )
        elif par.u[i,j]<0:
            rhodele[j] = par.rho[i,j]*par.delta[i,j]*par.e[i,j] +                                                  \
                         par.dt*(par.Ce[i,j] -                                                                     
                                    (                                                                             
                                    par.delta[i,j+1]*par.u[i,j+1]*(par.rho[i,j+1]*par.e[i,j+1] + par.p[i,j+1])-     
                                    par.delta[i,j]  *par.u[i,j]  *(par.rho[i,j]  *par.e[i,j]   + par.p[i,j])    
                                    )/                                                                            
                                    (par.x[j+1]-par.x[j])                                                         
                                )
    # left boundary
    rhodele[par.jlb] = par.rho[i,par.jlb]*par.delta[i,par.jlb]*par.e[i,par.jlb] + par.dt*par.Ce[i,par.jlb]    
    # right boundary
    rhodele[par.jrb] = par.rho[i,par.jrb]*par.delta[i,par.jrb]*par.e[i,par.jrb] + par.dt*par.Ce[i,par.jrb]
    e = rhodele[:]/par.rhodel[i+1,:]
    return e

def get_T():
    # Diagnostic
    T = np.zeros([par.xlen],dtype='float') # (t+1)
    # interior
    T[par.jmin:par.jmaxp1] = (1/par.cv)*(par.e[i+1,par.jmin:par.jmaxp1]-0.5*par.u[i+1,par.jmin:par.jmaxp1]**2)
    # left boundary (no flux)
    T[par.jlb] = T[par.jlb+1]
    # right boundary (no flux)
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
    for j in range(par.jmin,par.jmaxp1):
        if par.u[i,j]>=0:
            xidot[j] =  (
                        par.Fnet[i+1,j]*(par.ps[i+1,j]*par.L)/(par.R*par.Ts[i+1,j]**2*par.Co)  
                        + par.g*(                                                                 
                                par.rho[i+1,j]  *par.delta[i+1,j]  *par.u[i+1,j]              
                               -par.rho[i+1,j-1]*par.delta[i+1,j-1]*par.u[i+1,j-1]        
                                )/(par.x[j]-par.x[j-1])
                        )/\
                        (                                                                       
                        par.delta[i+1,j]*(par.g + 
                                         (par.ps[i+1,j]*par.L**2)/(par.R*par.Ts[i+1,j]**2*par.Co)
                                         )
                        )
        elif par.u[i,j]<0:
            xidot[j] =  (
                        par.Fnet[i+1,j]*(par.ps[i+1,j]*par.L)/(par.R*par.Ts[i+1,j]**2*par.Co)  
                        + par.g*(                                                                 
                                par.rho[i+1,j+1]*par.delta[i+1,j+1]*par.u[i+1,j+1]           
                               -par.rho[i+1,j]  *par.delta[i+1,j]  *par.u[i+1,j]        
                                )/(par.x[j+1]-par.x[j])
                        )/\
                        (                                                                       
                        par.delta[i+1,j]*(par.g + 
                                         (par.ps[i+1,j]*par.L**2)/(par.R*par.Ts[i+1,j]**2*par.Co)
                                         )
                        )    
    # left boundary
    j = par.jlb
    xidot[j] =  (
                par.Fnet[i+1,j]*(par.ps[i+1,j]*par.L)/(par.R*par.Ts[i+1,j]**2*par.Co)  
               +par.g*(                                                                 
                        par.rho[i+1,j+1]*par.delta[i+1,j+1]*par.u[i+1,j+1]           
                        -par.rho[i+1,j]  *par.delta[i+1,j]  *par.u[i+1,j]        
                    )/(par.x[j+1]-par.x[j])
                )/\
                (                                                                       
                par.delta[i+1,j]*(par.g + (par.ps[i+1,j]*par.L**2)/(par.R*par.Ts[i+1,j]**2*par.Co))
                )    
    # right boundary
    j = par.jrb
    xidot[j] =  (
                 par.Fnet[i+1,j]*(par.ps[i+1,j]*par.L)/(par.R*par.Ts[i+1,j]**2*par.Co)  
                +par.g*(                                                                 
                        par.rho[i+1,j]  *par.delta[i+1,j]  *par.u[i+1,j]              
                        -par.rho[i+1,j-1]*par.delta[i+1,j-1]*par.u[i+1,j-1]        
                        )/(par.x[j]-par.x[j-1])
                )/\
                (                                                                       
                par.delta[i+1,j]*(par.g + (par.ps[i+1,j]*par.L**2)/(par.R*par.Ts[i+1,j]**2*par.Co))
                )

    print((par.ps[i+1,:]*par.L**2)/(par.R*par.Ts[i+1,:]**2*par.Co))
    return xidot

def get_initial_conditions():
    # Scenario: intially saturated atmosphere with uniform Ts,T,p,delta,etc.
    par.Ts[0,:]       = np.ones([par.xlen],dtype='float')*350 # uniform initial temp
    par.ps[0,:]       = np.ones([par.xlen],dtype='float')*4e4 # uniform initial pressure
    par.rhodel[0,:]   = (par.ps[0,:]-par.p0)/par.g            # uniform atmospheric mass
    par.u[0,:]        = np.ones([par.xlen],dtype='float')     # 1 m/s initial wind
    par.u[0,par.jlb:par.jlb+2]   = 0 
    par.u[0,par.jrb-1:par.jrb+1] = 0 
    par.T[0,:]        = par.Ts[0,:]                           # T = Ts
    par.e[0,:]        = par.cv*par.T[0,:] + 0.5*par.u[0,:]**2 
    par.p[0,:]        = par.ps[0,:]                           # p = ps
    par.rho[0,:]      = par.p[0,:]/(par.R*par.T[0,:])         
    par.delta[0,:]    = par.rhodel[0,:]/par.rho[0,:]          
    par.xidot[0,:]    = np.zeros([par.xlen],dtype='float')    # zero surface mass-flux
    #print(par.rhodel[0,:])
    #print(par.e[0,:])
    #print(par.p[0,:])
    #print(par.rho[0,:])
    #print(par.delta[0,:])
    
# Time Integration
par = parameters()
get_initial_conditions()

# Plot forward integration
fig, ax = plt.subplots(10,2,figsize=(4,10))    

for i in range(par.tlen-1): 
    # Get variables at current time step
    par.Cm[i,:]         = get_Cm()
    par.Cu[i,:]         = get_Cu()
    par.Ce[i,:]         = get_Ce()
    par.Fnet[i,:]       = get_Fnet()
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
    
    if i==(par.tlen-2):
        # Plot at end of the simulation
        ax[0,0].plot(par.x/par.rp,par.Ts[i,:]) 
        ax[1,0].plot(par.x/par.rp,par.ps[i,:]) 
        ax[2,0].plot(par.x/par.rp,par.rhodel[i,:]) 
        ax[3,0].plot(par.x/par.rp,par.u[i,:]) 
        ax[4,0].plot(par.x/par.rp,par.e[i,:]) 
        ax[5,0].plot(par.x/par.rp,par.T[i,:]) 
        ax[6,0].plot(par.x/par.rp,par.p[i,:]) 
        ax[7,0].plot(par.x/par.rp,par.rho[i,:]) 
        ax[8,0].plot(par.x/par.rp,par.delta[i,:])
        ax[9,0].plot(par.x/par.rp,par.xidot[i,:])

        ax[0,1].plot(par.x/par.rp,par.Ts[i+1,:]) 
        ax[1,1].plot(par.x/par.rp,par.ps[i+1,:]) 
        ax[2,1].plot(par.x/par.rp,par.rhodel[i+1,:]) 
        ax[3,1].plot(par.x/par.rp,par.u[i+1,:]) 
        ax[4,1].plot(par.x/par.rp,par.e[i+1,:]) 
        ax[5,1].plot(par.x/par.rp,par.T[i+1,:]) 
        ax[6,1].plot(par.x/par.rp,par.p[i+1,:]) 
        ax[7,1].plot(par.x/par.rp,par.rho[i+1,:]) 
        ax[8,1].plot(par.x/par.rp,par.delta[i+1,:])
        ax[9,1].plot(par.x/par.rp,par.xidot[i+1,:])

for c in range(0,2):    
    ax[0,c].set_xlabel('Ts (K)')
    ax[1,c].set_xlabel('ps (Pa)')
    ax[2,c].set_xlabel('rhodel')
    ax[3,c].set_xlabel('u (m/s)')
    ax[4,c].set_xlabel('e')
    ax[5,c].set_xlabel('T (K)')
    ax[6,c].set_xlabel('p (Pa)')
    ax[7,c].set_xlabel('rho')
    ax[8,c].set_xlabel('delta (m)')
    ax[9,c].set_xlabel('xidot')
plt.tight_layout()
plt.show()
    
    
