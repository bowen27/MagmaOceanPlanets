# Goal: simulate a 1D boundary layer of a magma ocean planet

# Import python libraries
import numpy as np

# Create parameter class

class parameters:
    # Spatial grid
    Rp   = 6314e3
    ndeg = 180
    ddeg = 1
    dx   = Rp*np.pi/ndeg
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
    
    def __init__(self):
        self.Ts = []
        self.ps = []
        self.rhodelta = []
        self.xidotdelta = []
        self.u = []
        self.e = []
        self.T = []
        self.p = []
        self.rho = []
        self.delta = []
        self.xidot = []
  
    def set_Ts(self, x): 
        # see python list methods
        self.Ts.clear()
        self.Ts.append(x)
    
# Test: Create a class, setup Ts as any array, and then access the last element of Ts   
par = parameters()
Ts = np.arange(parameters.xlen,dtype='float') 
par.set_Ts(Ts) 
# print(par.Ts[0])      # returns the entire array
# print(par.Ts[0][-1])  # returns the last value of the array 


      
