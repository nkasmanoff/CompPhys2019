

import numpy as np

import matplotlib.pyplot as plt

%matplotlib inline


# Given Constants
a = 1.
b = 3.
x_0 = 0.
y_0 = 0.
delta = 10 ** -10  # target accuracy per unit time
t_0 = 0.
t_f = 20.

def f(r):
    """
    
    Brussellator Equation
    
    
    """
    x = r[0]
    y = r[1]
    
    dx = 1 - (b+1)*x + a*x*x*y
    
    dy = b*x - (a * x*x*y)
    
    return np.array([dx,dy])


from scipy import copy
def solution(t_0,t_f):
    r = np.array([x_0,y_0])
    tpoints = [t_0]
    xpoints = [r[0]]
    ypoints = [r[1]]
    
    def BS_step(r,t,H):
        """
        Computes a BS step  with richardson extrapolation until 
        
        target accuracy is reached. 
        
        """
        
        def modified_midpoint(r,n):
            
            r = copy(r)
            
            h = H / n
            k = r + .5 * h * f(r)
            r += h * f(k)
            
            for i in range(n - 1):
                k += h * f(r)
                r += h*f(k)
                
            return 0.5 * (r + k + 0.5 * h * f(r))
        
        
        def compute_R_n(R1,n):
            """
            Calculate the nth row of richardson extrapolation given R1 is inputed, 
            and is simply the modified midpoint method. 
            
            By emplying recursion, we include in this function 
            the function for calculating R_n_m
            
            """
            

            def R_n_m(m):
                """
                Computes R_n,m
                :param m: integer <= n
                :return: the vector R_n,m
                """
                return R2[m - 2] + (R2[m - 2] - R1[m - 2]) / ((n / (n - 1)) ** (2 * (m - 1)) - 1)

            #stops after n > 8
        
            if n > 8:
                #recursion
                r1 = BS_step(r,t,H / 2)
                return BS_step(r1, t + H / 2, H / 2)
            
            else: #if not, perform the extrapolation up to n = 8. 
                # Compute R_n,1
                R2 = [modified_midpoint(r, n)]
                # Compute the rest of the row
                for m in range(2, n + 1):
                    R2.append(R_n_m(m))
                    
                    
                #
                R2 = np.array(R2)
                
                #calculation of error for this calculation
                error_vec = (R2[n - 2] - R1[n - 2]) / ((n / (n - 1)) ** (2 * (n - 1)) - 1)
                
                error  = np.sqrt(error_vec[0] ** 2 + error_vec[1] ** 2)
                
                
                # if error is less than goal accuracy, stop. Otherwise, repeat with more steps!
                
                target = H * delta
                
                if error < target:
                    tpoints.append(t + H)
                    xpoints.append(R2[n - 1][0])
                    ypoints.append(R2[n - 1][1])
                    return R2[n - 1]
                else:
                    return compute_R_n(R2, n + 1)


        return compute_R_n(np.array([modified_midpoint(r, 1)], float), 2)

    BS_step(r, t_0, t_f - t_0)
    return tpoints, xpoints, ypoints

t, x, y = solution(t_0, t_f)
plt.plot(t, x, 'b.-',label = 'x')
plt.plot(t, y, 'r.-',label = 'y')
plt.xlabel('$t$',fontsize=20)
plt.legend()
plt.ylabel('Concentrations',fontsize=20)
plt.savefig('../bin/adaptivebs.png')