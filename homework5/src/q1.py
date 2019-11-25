"""
Solution to question 1, 9.7 in Newman, using relaxation method to solve kinematic equation for a
ball thrown in the air and landing in the same position ten seconds later. 

Bonus output at the end which calculated the initial condition (velocity at which ball is thrown up)
"""


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('white')

#constants 
g = 9.8 #gravity
M = 100 # number of time slots
tf = 10 #total amount of time

h = tf/ M # time step to satisfy grid dimensions 

target  = 1e-6 #target acc

#initalize array to save x position
x = np.zeros([M+1],float)
xprime = np.zeros([M+1],float)


#main loop

delta = 1.0

while delta > target:
    for i in range(1,M):
        if i == 0 or i == M:
            xprime[i] = x[i]  #set by boundary conditions
        else: #going to have to apply jacobi / relaxation method
            xprime[i] = (g*h*h + x[i+1] + x[i-1]) / 2
        #calculate max difference        
        delta = max(abs(x - xprime))
        #swap the two arrays
        
        x[i],xprime[i] = xprime[i],x[i]


plt.plot(x,'ro--')
plt.ylabel('$x(t) $ (m)',fontsize=25)
plt.xlabel('$t$ (deci-s) ',fontsize=25)
plt.savefig('../bin/balltrajectory.png')

print("Now solving for initial conditions...")
v0 = (x[1] - x[0]) / h
vf = (x[-1] - x[-2]) /h




print("The initial velocity of this problem is ", v0)
print("The final velocity of this system is ", vf)