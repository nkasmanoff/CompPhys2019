"""

This is just a bonus script to show how without office hours and I tried to 
apply gradient descent to the Schecter function's parameters in a linear fashion,

it fails terribly. Imagine for instance the inital step size is .01. 

Mstar0 = 1e11, and Mstar1 = 1e11 + .01. That is such a tiny change!

if you execute this function as python linear_grad_descent.py, 
you will not obtain an answer.




"""
import numpy as np
import pandas as pd
cosmos = pd.read_table("../dat/smf_cosmos.dat", sep="\s+",header=None)
cosmos.rename(columns={0:"log M_gal [dex]",1:"n(M_gal) [1/dex/Volume]",2:"error in n(M_gal)"},inplace=True)


def schecter(M,phi,Mstar,alpha):
    """Calculates the schecter function. 
    
    """
    from numpy import log,exp,power
#    M = power(10.0,M)
 #   Mstar = power(10.0,Mstar)
 #   phi = power(10.0,phi)
    return phi *  exp(-M/Mstar) * power(M/Mstar,alpha+1) * log(10) #log 10 smooths binning


def chisquared(phi,Mstar,alpha,f=schecter,cosmos = cosmos):
    
    M = cosmos['log M_gal [dex]']
  #  print("Mstar: ", Mstar)
    observed = f(M,phi,Mstar,alpha)
    
    expected = cosmos['n(M_gal) [1/dex/Volume]']
    
    error = cosmos['error in n(M_gal)']
    
    return sum((observed - expected)**2 / error **2) 

f = chisquared

def grad_desc(f,a0,b0,c0,gamma=.01,talk=True,initial_step = 1e-3):
    
    tol = gamma
    a = a0 + initial_step
    b = b0 + initial_step
    c = c0 + initial_step
    chi_squareds = []
    while True:
        
        chi_squareds.append(f(a,b0,c0)) #obtain the chi squared term
        deriv_a = (f(a,b0,c0) - f(a0,b0,c0)) / (a - a0)
        
        deriv_b = (f(a0,b,c0) - f(a0,b0,c0)) / (b - b0)

        deriv_c = (f(a0,b0,c) - f(a0,b0,c0)) / (c - c0)
        
        a1 = a - gamma * deriv_a
        b1 = b - gamma * deriv_b
        c1 = c - gamma * deriv_c
        if talk:
            print("chi squared =",f(a,b0,c0))
        if abs(a1 - a) + abs(b1-b) + abs(c1 - c) < tol:
            return a1,b1,c1,chi_squareds
        
        if type(f(a,b0,c0)) == type(np.nan):
            return a1,b1,c1,chi_squareds

        a0 = a
        a = a1
        
        b0 = b
        b = b1
        
        c0 = c
        c = c1
        
a1,b1,c1,chi_squareds = grad_desc(chisquared,-.01,10*7,-1.01,gamma=1e-24,initial_step=10000)