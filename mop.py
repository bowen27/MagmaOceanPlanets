# Goal: simulate a 1D boundary layer of a magma ocean planet
# ** denotes a problem or code requiring verification

# To do:
# get_e , get_T, get_p, get_rho, get_delta, get_xidot
# get_x_bc for each variable (assuming no flux bc)
# get_Cm, get_Cu, get_Ce

# Import python libraries
import numpy as np

# Create parameter class

class parameters:
    # Spatial grid
    rp   = 6314e3 
    ndeg = 180    
    ddeg = 1
    dx   = rp*np.pi/ndeg
    x    = np.arange(0,ndeg+ddeg,ddeg)*dx
    xlen = len(x)
    # Temporal grid
    dt   = 0.01
    tmin = 0
    tmax = 1
    t    = np.arange(tmin,tmax+dt,dt)
    tlen = len(t)
    # CFL criterion must be met
    if not dt/dx < 0.4:
        print('CFL criterion is not met.')
        quit()
    # Thermodynamics **
    cp = 1000     # J/kg/K
    L  = 2260e3   # J/kg
    M  = 29e-3    # kg/mol
    R  = 8.314/M  # J/kg/K
    # Clausius-Clapeyron **
    pcc = 1e6
    Tcc = 1e6
    # Free Atmosphere
    p0 = 1e5
    g  = 9.81
    
    def __init__(self):
        self.Ts , self.ps   , self.rhodel, self.xidotdel,\
        self.u  , self.e    , self.T     , self.p       ,\
        self.rho, self.delta, self.xidot , self.Fnet    ,\
        self.Cm , self.Cu   , self.Ce                    \
        = (np.zeros([parameters.tlen,parameters.xlen],dtype='float') for i in range(15))
                
# Test: Create a class, update Ts, and then access the last element of Ts on the first time step   
par = parameters()
# Ts = np.arange(par.xlen,dtype='float') 
# par.Ts[0,:] = np.arange(par.xlen,dtype='float') 
# print(par.Ts[0,-1])

# Upwind Scheme ** 
def dfdx(f,x):
    # Derivatives are defined at same gridpoints as other data.
    # 
    # ** left boundary condition is needed
    dfdx = np.zeros([par.xlen],dtype='float')
    for j in range(1,par.xlen):
        df   = f[j] - f[j-1]
        dx   = x[j] - x[j-1]
        dfdx = df/dx
    return dfdx

# Saturation Vapor Pressure **
def esat(x):
    # Input x could be T or Ts
    y = np.zeros([par.xlen],dtype='float')
    y = par.pcc*np.exp(-(par.L/par.R)*(1/x - 1/par.Tcc))
    return y 

def get_Ts(): # ** verify
    # Prognostic
    Ts = np.zeros([par.xlen],dtype='float')
    Ts = par.Ts[i,:] + par.dt/par.cp*\
            (\
             par.Fnet[i,:]-par.L*par.xidot[i,:]*par.delta[i,:]\
            )     
    return Ts

def get_ps(): # ** verify
    # Diagnostic
    ps = np.zeros([par.xlen],dtype='float')
    ps = esat(par.Ts[i,:])
    return ps

def get_rhodel(): # ** verify
    # Diagnostic
    rhodel = np.zeros([par.xlen],dtype='float')
    rhodel = (par.ps[i,:]-par.p0)/par.g
    return rhodel

def get_xidotdel(): # ** verify
    # Prognostic 
    y = np.zeros([par.xlen],dtype='float')
    for j in range(par.xlen):
        if j!=0:
            y[j] = (par.rhodel[i+1,:]-par.rhodel[i,:])/par.dt + \
                   (par.u[i,j]*par.rhodel[i,j]-par.u[i,j-1]*par.rhodel[i,j-1])/par.dx
        else: 
            # boundary condition at left wall
    return y

def get_u(): # ** verify
    # Prognostic
    rhodelu = np.zeros([par.xlen],dtype='float') # (t+1)
    u = np.zeros([par.xlen],dtype='float')       # (t+1)
    
    for j in range(par.xlen):
        if j!=0:
            rhodelu[j] = par.rho[i,j]*par.delta[i,j]*par.u[i,j] +                                     \
                         par.dt*(par.Cu[i,j] -                                                        \
                                     (                                                                \
                                     par.delta[i,j]  *(par.rho[i,j]  *par.u[i,j]**2   + par.p[i,j])-  \ 
                                     par.delta[i,j-1]*(par.rho[i,j-1]*par.u[i,j-1]**2 + par.p[i,j-1]) \
                                     )/                                                               \
                                     (par.x[j]-par.x[j-1])                                            \
                                 )
        else: 
            # boundary condition at left wall
    u = rhodelu[:]/par.rhodel[i+1,:]
    return u

def get_e():
    
# Time Integration
# ** initial conditions must be set
for i in range(par.tlen):
    
    # Get variables at the next time step
    par.Ts[i+1,:]       = get_Ts()
    par.ps[i+1,:]       = get_ps()
    par.rhodel[i+1,:]   = get_rhodel()
    par.xidotdel[i+1,:] = get_xidotdel()
    par.u[i+1,:]        = get_u()
    par.e[i+1,:]        = get_e()
    par.T[i+1,:]        = get_T()
    par.p[i+1,:]        = get_p()
    par.rho[i+1,:]      = get_rho()
    par.delta[i+1,:]    = get_delta()
    par.xidot[i+1,:]    = get_xidot()
