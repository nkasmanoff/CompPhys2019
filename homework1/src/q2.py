import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set_style('whitegrid')
sns.set_palette('husl')
from numpy import exp



eminust = lambda t: exp(-t) #defined function for this problem 
def integrate(f,a,b,N,method = 'all'):
    """
    Integrates a function f from a to b, given N estimation bins, with the 
    provided method. 
    
    
    Parameters
    ----------
    
    f : function
        Function to integrate. 
    a : float
        Starting point of integral.
    b : float
        End point of integral. 
    N : int
        Number of bins used to estimate the function's area under curve. 
    method : str
        What computational method to use, be it midpoint, trapezoid, simpson'. 
        
    Returns
    -------
    
    
    I_method : float
        Integral sum with that method's solution
    """
    #trapezoid method here. 
    h = (b - a) / (N)
    
    midpoint_sum = sum([h*f((k+.5)*h) for k in range(1,N)])  #   expected error ~h^1
    I_midpoint = midpoint_sum
            
    trap_sum = sum([f(a + k*h) for k in range(1,N)])
    I_trapezoid = h*(f(a) / 2 + f(b) / 2 + trap_sum)  #   expected error ~h^2
              
        
    s = f(a) + f(b) + 4*f(b-h)
    for k in range(1,N//2):
        s += 4*f(a + (2*k-1)*h) + 2*f(a+2*k*h)
    I_simpson = (h/3)*s  #   expected error ~h^4
        
    if method == 'midpoint':
        return np.float32(I_midpoint)
    if method == 'trapezoid':
        return np.float32(I_trapezoid)
    #still need to work on trapezoid rule .
    if method == 'simpson':
        return I_simpson
    
    if method == 'all':
        return [I_midpoint,I_trapezoid,I_simpson]
    


f = eminust
true_integral = -eminust(1) +eminust(0) #this integral doesn't have an exact solution? 
a = 0
b = 1

Ns = 10**np.arange(0,7)
calculated_integral = []

for N in Ns:
    
    calculated_integral.append(integrate(f=eminust,a=a,b=b,N=N,method='all'))
error = abs(np.array(calculated_integral)  - true_integral)
plt.figure(figsize=(10,10))

#%%
plt.plot(Ns,error[:,0],'bo-',label = 'midpoint')
plt.plot(Ns,error[:,1],'ro-',label = 'trapezoid')
plt.plot(Ns,error[:,2],'go-',label = 'simpson')
#plt.xlim([min(Ns),max(Ns)])
plt.legend()
plt.title("Absolute Error vs Bin Count of e^-t" + " from 0 to 1"  +" (log-log)",fontsize = 25)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('N (total bins used)',fontsize=15)
plt.ylabel('ε (relative error)',fontsize = 15)


plt.savefig('../bin/q2plots.png')