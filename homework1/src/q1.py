import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set_style('whitegrid')
sns.set_palette('husl')


def derivs(f,x,h,method = 'all'):
    """
    Using forward differencing, central differencing, or richardson extraplation differencing,
    calulate the derivative of a given function f. 
    
    
    We seek to model the error behavior of these various methods for analytically solvable functions like cosine and
    exponentials, and from there interpret the manifestation of these errors over different regimes. 
    

    Parameters
    ----------
    
    f : function
        function of interest, i.e np.cos, np.exp, etc. 
    x : float
        x value of the function to compute the derivative at. 
    h : float
        Step size for the derivative calulcation, tends towards 0. 
    method : str
        desired method of differentiation. Default setting is all, and will compute all three derivative calculations at one. 
        
    Returns
    -------
    
    f or c + e + deriv : float
        Derivative of function f at position x given step size h, for the provided method. 
        Note if method is set to all, will return a list of all three in ascending order of precision. 
        f = forward, c = central, e = extrapolation
        
    """
    #make sure to convert all values to float 32
    fderiv = (np.float32(f(np.float32(x) + np.float32(h))) - np.float32(f(np.float32(x)))) / np.float32(h)
    cderiv = np.float32((np.float32(f(np.float32(x) + np.float32(h) )) - np.float32(f(np.float32(x) - np.float32(h)))) / (np.float32(2)*np.float32(h)))
    ederiv = ( np.float32(-f(np.float32(x) + np.float32(2)*np.float32(h))) +
              np.float32(8) * np.float32(f(np.float32(x) + np.float32(h))) - np.float32(8) * 
    np.float32(f(np.float32(x) - np.float32(h))) + np.float32(f(np.float32(x) - np.float32(2)*np.float32(h))) )  /(np.float32(12)*np.float32(h))

    if method == 'forward': #   expected error ~h^1
        return fderiv
    if method == 'central': #   expected error ~h^2
        return cderiv
    if method == 'extrapolation': #   expected error ~h^4
        return ederiv
    if method == 'all':
        return [fderiv,cderiv,ederiv]


from numpy import cos,exp
values = [.1,1] #+ 2*np.pi Is the error the same? 
functions = [cos,exp]

hvals = .1**np.arange(0,9,dtype='single')
plt.figure(figsize=(10,12))
i = 1
for _, value in enumerate(values):
    for _, function in enumerate(functions):
        func_str = str(function).split(' ')[1].split("'")[1]
        calculated_derivative = []

        for h in hvals:
            calculated_derivative.append(derivs(function,value,h,method='all'))

        if function == cos:
            true_deriv = -np.sin(value)
        if function == exp:
            true_deriv = exp(value)
        true_deriv = np.float32(true_deriv)
        error = abs(calculated_derivative - true_deriv) #absolute error
       # print("i+1: ", i)
        plt.subplot(2, 2, i)
        
        plt.plot(hvals,error[:,0],'b-o',label = 'forward')
        plt.plot(hvals,error[:,1],'r-o',label = 'central')
        plt.plot(hvals,error[:,2],'g-o',label = 'extrapolation')
        #plt.xlim([max(hvals),min(hvals)])
        plt.legend()
        plt.title("Absolute Error vs Step Size of "+ func_str + " @ x = " + str(value)  +" (log-log)")#,fontsize = 10)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('h (step size)')
        plt.ylabel('ε (absolute error)')

        i += 1


plt.savefig('../bin/q1plots.png')