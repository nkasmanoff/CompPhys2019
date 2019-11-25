"""

-Create a grid to house the particles.dat file 

- Use relaxation method to solve for potential 
	- wait 15 minutes ...

- Try over-relaxed gauss seidel, and use golden ratio for finding optimal omega. 



"""


### A ###

import numpy as np
import matplotlib.pyplot as plt

particles = np.loadtxt('../dat/particles.dat')
plt.figure()
plt.scatter(particles[:,0],particles[:,1])
plt.xlabel("$x$",fontsize=25)
plt.ylabel("$y$",fontsize=25)
plt.savefig('../bin/pointcloud.png')


class Grid():
    def __init__(self,length,width):
        self.length = length
        self.width = width
        self.box = np.zeros([length,width])
    
    def CIC_interp(self,particles):
        """
        
        """
        for particle in particles:
            #closest values in X and y
            NNx = int(np.floor(particle[0]) ), int(np.ceil(particle[0]) )
            NNy = int(np.floor(particle[1])), int(np.ceil(particle[1]))
            
            contrib_lower_X = abs(1 - particle[0] + NNx[0])
            contrib_upper_X = 1 - contrib_lower_X

            contrib_lower_Y = abs(1 - particle[1] + NNy[0])
            contrib_upper_Y = 1 - contrib_lower_Y
            
            self.box[NNx[0]-1,NNy[0]-1] += contrib_lower_X*contrib_lower_Y
            self.box[NNx[0]-1,NNy[1]-1] += contrib_lower_X*contrib_upper_Y
            self.box[NNx[1]-1,NNy[0]-1] += contrib_upper_X*contrib_lower_Y
            self.box[NNx[1]-1,NNy[1]-1] += contrib_upper_X*contrib_upper_Y
            
        return self.box


grid = Grid(100,100)
plt.figure()
plt.imshow(grid.CIC_interp(particles))
plt.xlabel("$x$",fontsize=25)
plt.ylabel("$y$",fontsize=25)
plt.savefig('../bin/cic_grid.png')



###### B

cell = grid.CIC_interp(particles)
# Constants
M = 100
target = 1e-10   # Target accuracy
epsilon0 = 1 #8.85e-12

# Create arrays to hold potential values
phi = np.zeros([M,M],float)
phiprime = np.empty([M,M],float)

a = 1
# Main loop
delta = 1.0

rho = grid.CIC_interp(particles)
steps = 0
print("Now commencing ordinary relaxation method . ")
while delta > target:
    for i in range(0,M-1):
        for j in range(0,M-1):
        #    print(delta)
            phiprime[i,j] = (phi[i+1,j] + phi[i-1,j] + phi[i,j+1] + phi[i,j-1])/4 + (1/4 / epsilon0) * rho[i,j] 
            #calculate max difference        
            #swap the two arrays
            phi[i,j],phiprime[i,j] = phiprime[i,j],phi[i,j]
      #  print(np.mean(phi[i,:]))
    delta = np.max(abs(phi - phiprime))
    steps +=1
    
    if steps % 100 == 0:
        print("Steps = ", steps )
        print("Delta = ", delta)

print("Done! Solving for potential now ")
plt.figure()
plt.imshow(phi)
plt.xlabel("$x$",fontsize=25)
plt.ylabel("$y$",fontsize=25)
plt.savefig('../bin/phi_relaxed.png')


#### C 

print("Now trying gauss seidel")


def gauss_seidel_overrelaxation(omega,returnphi=False):
    #Apply gauss seidel over relaxation method for different values of omega 
    
    delta =  1.0 #make delta larger than target, allows loop to begin.
    M = 100
    a = 1
    
    target = 1e-10  # Target accuracy
    epsilon0 = 1 #8.85e-12
    # Create array[ no s since gauss seidel!!!] to hold potential values
    phi = np.zeros([M,M],float)
    rho = grid.CIC_interp(particles)
    # rho = -1.6e-19 * rho #unit charge, if useful! 
    steps = 0
    while delta > target:
        # Main loop
        delta = 0.0
        for i in range(0,M-1):
            for j in range(0,M-1):
                
                phi_old = phi[i,j]
                
                phi_new = (1 + omega) * .25 * (phi[i-1,j] + phi[i+1,j] + phi[i,j-1] + phi[i,j+1] + rho[i,j]  * a*a / epsilon0)  - omega*phi_old
                
                phi[i,j] = phi_new 
                
                difference = abs(phi_new - phi_old)
                
                if difference > delta: 
                    delta = difference
        steps+=1
        
    if returnphi:
        return steps,phi
    
    return steps

z = (1+np.sqrt(5))/2       # Golden ratio
accuracy= .001

# Function to calculate
omegamins = []
stepcounts = []
f = gauss_seidel_overrelaxation
# Initial positions of the four points
w1 = .5
w4 = .999
w2 = w4 - (w4-w1)/z
w3 = w1 + (w4-w1)/z

# Initial values of the function at the four points
f1 = f(w1)
f2 = f(w2)
f3 = f(w3)
f4 = f(w4)

# Main loop of the search process
while w4-w1>accuracy:
    print(w1,w2,w3,w4)
    if f2<f3:
        w4,f4 = w3,f3
        w3,f3 = w2,f2
        w2 = w4 - (w4-w1)/z
        f2 = f(w2)
    else:
        w1,f1 = w2,f2
        w2,f2 = w3,f3
        w3 = w1 + (w4-w1)/z
        f3 = f(w3)
        
    stepounts.append(f(0.5*(w1+w4)))
    omegamins.append(0.5*(w1+w4))
# Print the result
print("The minimum falls at",0.5*(w1+w4))


plt.figure()

plt.plot(omegamins,stepcounts)
plt.xlabel("omega",fontsize = 25)
plt.ylabel("Steps",fontsize = 25)
plt.savefig('../bin/stepsomega.png')



