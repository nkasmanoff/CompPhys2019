"""
Here is the code used to answer question 2 of homework 2. 

"""
import numpy as np


def binary(f,a,b,tol = 1e-6):
    """

    Find the root of function f give the starting points a and b. 
    
    The goal is for f(a) and f(b) to have opposite signs, and if they don't redo this step.
    
    But if they do have opposite signs, find the midpoint btw a and b, and evaluate. if this f(midpoint) is below tolerance,
    return.
    
    
    if not, determine the sign of this midpoint, and repeat the search. 

    """
    
    from numpy import sign #tells you if a positive or negative
    
    if f(a) * f(b) > 0:
        return "The selected endpoints are the same sign, please try different inputs!"

    
    c = (a + b)/2
    while abs(f(c)) > tol:
  
        if f(a) * f(c) > 0: 
            a = c
          #  print("Update a")
        elif f(a) * f(c) < 0:
            b = c
          #  print("Update b")
        else:
            return c
        
        c = (a + b)/2

    return c#,f(c)


f = lambda x: 5*np.exp(-x) + x - 5#5e−x + x − 5


minima = binary(f=f,a=3,b =6 ,tol = 1e-6)

print("The minima of 5e^(−x) + x − 5 is ", minima)
l = 5.02 *10**-7

kb = 1.380649e-23 #×10−23 J⋅K−1
#T is kelvin. 
h =  6.62607015e-34#−34  J⋅s
c = 2.99e8 #m/s
x = minima

b = (h*c)/(kb*x)

print("Wien's Displacement Constant b = ", b)


T = b/l

print("The surface temperature of the Sun is approximately ", T, ' K.')