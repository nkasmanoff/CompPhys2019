import numpy as np
import matplotlib.pyplot as plt


N = int(1e4) #10000 points

a = 1664525
c = 1013904223
m = 4294967296
x = 1
y = 42
U1 = []
U2 = []
for i in range(N):
    x = (a*x+c)%m
    U1.append(x)
    y = (a*y+c)%m
    U2.append(y)
 
U1 = np.divide(U1,m)
U2 = np.divide(U2,m)   


from numpy import sqrt,log,cos,pi,sin


# transformation function
def box_muller(u1,u2):
    z1 = sqrt(-2*log(u1))*cos(2*pi*u2)
    z2 = sqrt(-2*log(u1))*sin(2*pi*u2)
    return z1,z2


z1, z2 = box_muller(U1,U2)

from numpy.random import randn


standard_norm = []
i = 0
while i < N:

    standard_norm.append(randn())
    i +=1

plt.figure()
plt.hist(z1,bins=100
,label = 'Noah')

plt.hist(standard_norm,bins=100
         ,alpha = .75,label = 'numpy')
plt.yscale('log')
plt.legend()
plt.savefig('../bin/noahvsnumpy.png')



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

c = dft(z1)

plt.figure()
plt.plot(z1,'.--',linewidth = .1
    )

plt.savefig('../bin/whitenoise.png')



power = np.real(c)**2 + np.imag(c)**2
plt.figure()
plt.plot(power,'.') #k is just the index number. 

plt.ylabel("$log(|c(k)|^2)$",fontsize=20)
plt.xlabel('$log(k)$',fontsize=20)
print("Power Spectrum of Fourier Transform of Z1")
plt.xscale('log')
plt.yscale('log')
plt.savefig('../bin/gauss_powspec.png')


plt.figure()
plt.plot(random_walk)
plt.plot(z2.cumsum())
plt.xlabel("Step i", fontsize = 25)
plt.ylabel("$X(i)$", fontsize = 25)
plt.savefig('../bin/gauss_RW.png')


c = dft(random_walk)


power = np.real(c)**2 + np.imag(c)**2
plt.figure()
plt.plot(power) #k is just the index number. 
plt.grid()
plt.ylabel("$log(|c(k)|^2)$",fontsize=20)
plt.xlabel('$log(k)$',fontsize=20)
print("Power Spectrum of FT Gaussian Random Walk")
plt.xscale('log')
plt.savefig('../bin/gauss_RW_powspec.png')
plt.yscale('log')
