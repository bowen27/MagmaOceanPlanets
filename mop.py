# Goal: simulate a 1D boundary layer of a magma ocean planet

# Import python libraries
import numpy as np

# Create parameter class

class parameters:
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
        self.Ts.append(x) # see python list methods
    
par = parameters()
x = np.ones([3],dtype='float') # suppose this is the new array for Ts
par.set_Ts(x) # we set Ts

# To access the array, print(par.Ts[0])
# However, how do you access a single value of the np array?
print(par.Ts[0])    # returns [1. 1. 1.]
print(par.Ts[0][0]) # returns 1.0, but is not very clean


      
