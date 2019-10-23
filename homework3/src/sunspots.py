"""

Problem 7.2 load in the sunspots.txt file, visually examine the data to guess the length of a sunspot cycle, 
then obtain the power spectrum of the dft of the data and convert this peak wavenumber to wavelength and compare. 



"""

import numpy as np
import pandas as pd
import os 
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_palette('husl')
sns.set_style('whitegrid')


data = np.loadtxt('../dat/sunspots.txt',float)
y = data[:,1]

plt.figure()
plt.plot(y)
plt.ylabel("Sunspots",fontsize=25)
plt.xlabel('Month',fontsize=25)
plt.savefig('../bin/sunspots.png')
plt.figure()

plt.plot(y[0:200])
plt.ylabel("Sunspots",fontsize=25)
plt.xlabel('Month',fontsize=25)
plt.savefig('../bin/zoominsunspots.png')


def dft(y):
    """
    Code to calculate the discrete fourier transform c_k.
    
    Note that there will both a real and imaginary value for every c, 
    
    so to find the maximum k we need to calculate the power spectrum .
    
    """
    N = len(y) #we create the interval as the length of the dataset.. 
    
    c = np.zeros(N,complex) #list of fourier coef.

    for k in range(N): # Iterate over possible k's
        for n in range(N): #Iterate over all n's
            c[k] += y[n]*np.exp(-2j*np.pi*k*n/N) #create a sum to solve for c[k]
    return c

c = dft(y)

power = np.real(c)**2 + np.imag(c)**2

plt.figure()
plt.plot(power) #k is just the index number. 

plt.ylabel("$|c(k)|^2$",fontsize=20)
plt.xlabel('$k$',fontsize=20)
plt.savefig('../bin/powerspec.png')



plt.figure()
plt.plot(power[5:30]) #k is just the index number. 

plt.ylabel("$|c(k)|^2$",fontsize=20)
plt.xlabel('$k$',fontsize=20)
plt.savefig('../bin/zoomedinpowerspec.png')

k_max = list(power[5:30]).index(max(power[5:30])) + len(power[0:5])

print("The peak wavelength is " , k_max)