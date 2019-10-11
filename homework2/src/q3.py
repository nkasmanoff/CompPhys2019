"""
Here is the code used to answer question 3 of homework 2. 

"""


#first up is part a, applying gradient descent to a simple function.

import numpy as numpy 
import pandas as pd
import matplotlib.pyplot as plt
print("Part I: applying grad descent to minimize f(x,y) = (x-2)^2 + (y-2)^2")



def easy_2D_grad_descent(f,x0,y0,tol=1e-3,h=1e-6,gamma=.1):
    """Gradient descent function for the easy to evaluate function given as the starting point of q3. 
    """
    
    #first, define the initial values given as x and y
    
    x = x0
    y = y0
    fs = []
    while True:
        fs.append(f(x,y))
        fprimex = (f(x=x + h,y=y ) - f(x=x,y=y )) / (h)
        fprimey = (f(x=x,y=y+h ) - f(x=x,y=y ) )/ (h)
   #     fprimez = f(x=x,y=y +h) - f(x=x,y=y ) / (h)

        xnew = x- gamma*fprimex
        ynew = y - gamma*fprimey
    #    znew = z - gamma*fprimez
       # print('x',x)
       # print('xnew',xnew)
        if abs(xnew-x) + abs(ynew - y)<tol: #+ abs(znew - z) < tol:
            break
        x = xnew
        y = ynew
    return x,y, fs



f = lambda x,y: (x-2)**2 + (y-2)**2
x,y,fs = easy_2D_grad_descent(f,20,-30)

print("The minima of f is ", fs[-1], " at the value of x = ",x, ' and y=',y)
print("Sending plot of value over iteration to bin")
plt.figure(figsize=(10,8))
plt.plot(fs,linewidth=5)
plt.ylim([0,1000])
plt.xlabel("Step # (log)",fontsize=20)
plt.xscale('log')
plt.ylabel("$f(x,y)$",fontsize=20)
plt.savefig('../bin/q3plot1.png')


print("Now Schecter Function Time... ")

cosmos = pd.read_table("../dat/smf_cosmos.dat", sep="\s+",header=None)
cosmos.rename(columns={0:"log M_gal [dex]",1:"n(M_gal) [1/dex/Volume]",2:"error in n(M_gal)"},inplace=True)


def schecter(M,phi,Mstar,alpha):
    """Calculates the schecter function. 
    
    """
    from numpy import log,exp,power
    M = power(10.0,M)
    Mstar = power(10.0,Mstar)
    phi = power(10.0,phi)
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
            print(f(a,b0,c0))
        if abs(a1 - a) + abs(b1-b) + abs(c1 - c) < tol:
            return a1,b1,c1,chi_squareds
        
        a0 = a
        a = a1
        
        b0 = b
        b = b1
        
        c0 = c
        c = c1
        
        
a1,b1,c1,chi_squareds = grad_desc(chisquared,-2.5,10.5,-1.01,gamma=1e-5,initial_step=.01)


plt.figure(figsize=(10,8))
#Plot to show how well fit. 
plt.plot(cosmos['log M_gal [dex]'], schecter(cosmos['log M_gal [dex]'],phi=a1,Mstar=b1,alpha=c1),label='Schecter Best Fit')
#plt.plot(cosmos['M_gal'], schecter(cosmos['M_gal'],phi=,Mstar=1e11,alpha=-1.01))
#plt.plot(df['log M_gal [dex]'],df['n(M_gal) [1/dex/Volume]'],'bo')
plt.errorbar(x = cosmos['log M_gal [dex]'],y=(cosmos['n(M_gal) [1/dex/Volume]']),
             yerr=cosmos['error in n(M_gal)'],color='k',marker='o',ls='',label='COSMOS Data')
#plt.title("Fitted Schecter Function ($\chi^2 = 2.9$)",)
plt.xlabel("$log(M_{gal})$ [dex]",fontsize =20)
plt.ylabel('$n(M_{gal})$ [1/dex/Volume]',fontsize =20)
plt.yscale('log')
plt.legend(fontsize=20)
#plt.xscale('log')
plt.savefig('../bin/q3plot2.png')


#given different step sizes:
plt.figure(figsize=(10,8))

initial_steps = [.1,.01,.001]
print("Now sampling for different initial step sizes.")
for initial_step in initial_steps:
    a1,b1,c1,chi_squareds = grad_desc(chisquared,a0=-2.5,b0=10,c0=-1.1,gamma=1e-5,initial_step = initial_step,talk=False)
    print(len(chi_squareds))
    plt.plot(chi_squareds,label = 'Initial Step Size = ' + str(initial_step))
    plt.legend(fontsize = 20)
    plt.xlabel('Step Count (log)',fontsize = 20)
    plt.ylabel('$\chi^2$',fontsize = 20)
    
plt.xscale('log')
plt.savefig('../bin/q3plot3.png')



#given different step sizes:
plt.figure(figsize=(10,8))

initial_phi_mags = [-1,-2,-3]
print("Now sampling for different initial starting parameters.")
for initial_phi_mag in initial_phi_mags:
    a1,b1,c1,chi_squareds = grad_desc(chisquared,a0=initial_phi_mag,b0=10,c0=-1.1,gamma=1e-5,initial_step = .01,talk=False)
    print(len(chi_squareds))
    plt.plot(chi_squareds,label = 'Initial $log10(\phi_{*})$ = ' + str(initial_phi_mag))
    plt.legend(fontsize = 20)
    plt.xlabel('Step Count (log)',fontsize = 20)
    plt.ylabel('$\chi^2$',fontsize = 20)
    plt.xlim([0,100])
plt.xscale('log')
plt.xlim([0,100])

plt.savefig('../bin/q3plot4.png')


print("Done! ")



