# Goal: simulate a 1D boundary layer of a magma ocean planet
# ** denotes a problem or code requiring verification

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

# Upwind Scheme
def dfdx(f,x):
    # Derivatives are defined at same gridpoints as other data.
    # 
    # ** left boundary condition is needed
    dfdx = np.zeros([par.xlen],dtype='float')
    for i in range(1,par.xlen):
        df   = f[i] - f[i-1]
        dx   = x[i] - x[i-1]
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
    Ts = par.Ts[0][:] + par.dt/par.cp*\
            (\
             par.Fnet[0][:]-par.L*par.xidot[0][:]*par.delta[0][:]\
            )     
    return Ts

def get_ps(): # ** verify
    # Diagnostic
    ps = np.zeros([par.xlen],dtype='float')
    ps = esat(par.Ts[0][:])
    return ps

def get_rhodel(): # ** verify
    # Diagnostic
    rhodel = np.zeros([par.xlen],dtype='float')
    rhodel = (par.ps[0][:]-par.p0)/par.g
    return rhodel

def get_xidotdel(): # ** verify
    # Prognostic 
    xidotdel = np.zeros([par.xlen],dtype='float')
    xidotdel = par.rhodel[0][:]
    
# Time Integration
for i in range(parameters.tlen):
    t = parameters.t[i]
    