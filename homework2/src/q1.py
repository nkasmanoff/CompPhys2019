"""
Here is the code used to answer question 1 of homework 2. 

"""


#part 1, question 6.10 from Newman, defining the ordinary relaxation method. 



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def relaxation_method(x0,f,etol):
    """Starting at some initial value x0, 
    apply the relaxation method until convergence (below etol) to the function f. 
    
    """
    x = x0
    steps = 0
    while abs(x - f(x)) > etol:

        x = f(x)
      #  print(x,f(x))

        steps +=1
    return x,steps


c = 2
f = lambda x: 1-np.exp(-c*x)

x, N_steps = relaxation_method(10,f,1e-6)

print("Using simple relaxation, we recover x = ",x, " in ", N_steps, ' for f = ', f )


#part b of 6.10

cs = np.arange(0,3,.01)
xs = []
for c in cs:
    f = lambda x: 1-np.exp(-c*x)
    x = relaxation_method(10,f,1e-6)[0]
    xs.append(x)


plt.figure(figsize=(6,4))
plt.plot(cs,xs,linewidth = 3)
plt.xlabel('$c$',fontsize = 20)
plt.ylabel('$x$',fontsize = 20)

plt.savefig('../bin/q1plot1.png')




#now the actual code for question 6.11, the main part of this question!
def overrelaxation_method(x0,f,w,etol):
    """Starting at some initial value x0, 
    apply the relaxation method until convergence (below etol) to the function f. 
    
    """
    steps = 0
    x = x0
    while abs(x - f(x)) > etol:
      #  delta_x = x - f(x) #this is the current x
      #  x = x + (1 + w)*delta_x
        x = (1+w)*f(x) - w*x
   #     print(x,f(x))
        steps +=1
        if type(x/np.inf) == type(np.nan):  #for when this function blows up 
            return 
  #  print("Steps to converge = ", steps)
    return x,steps





x, N_steps = relaxation_method(10,f,1e-6)

print("Using simple relaxation, we recover x = ",x, " in ", N_steps, ' for f = ', f )

c = 2
f = lambda x: 1-np.exp(-c*x)

print("The overrelaxation method takes ", overrelaxation_method(1,f,.5,1e-6)[1], ' steps to converge.')
print("By comparison, the relaxation method took ", relaxation_method(1,f,1e-6)[1], " steps.")


#part c


ws = np.arange(0,1,.1)
convergence_steps = []
for w in ws:
    convergence_steps.append(overrelaxation_method(1,f,w,1e-6)[1])


plt.figure(figsize=(6,4))
plt.plot(ws,convergence_steps,linewidth = 3)
#plt.title("Steps to Convergence vs. w")
plt.xlabel('$\omega$',fontsize = 20)
plt.ylabel('Steps to Converge',fontsize=20)
plt.savefig('../bin/q1plot2.png')



