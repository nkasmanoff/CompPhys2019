import numpy as np
import matplotlib.pyplot as plt
from numpy.random import randint


#first, initialize the array
N = 20
lattice = np.empty((20,20),int)

for i in range(N):
    for j in range(N):
        #initialize randomly as 1 or -1
        if np.random.random() <= 0.5:
            lattice[i,j] = 1
        else: 
            lattice[i,j] = -1


plt.figure()
plt.imshow(lattice,cmap = 'jet')
plt.savefig('../bin/initial_lattice.png')


def energy(s):
    J = 1
    s1 = s[:-1,:]*s[1:,:]
    s2 = s[:,:-1]*s[:,1:]

    E = -J*(s1.sum() + s2.sum())

    return E



print("Initial energy of lattice", energy(lattice))


def metropolize(s,steps,T = 1):
    """
    Returns a time evolved ising grid, after the designated number of steps. 
    """
    J = 1
    kb = 1
    beta = 1

    # now do me
    eplot = []
    Mplot = []
    E1 = energy(s)
    M = s.sum()
    for k in range(steps): # over total steps
        i = randint(N) #select a random point
        j = randint(N)
        s[i,j] *=-1 #propose an energy change,   by -1 will immediately flip it!

        E2 = energy(s) #compute new energy

        dE = E2 - E1 #see what the change in energy is 

        if dE>0:  #if the system's total energy was not greater, 

            if np.random.random()<np.exp(-beta*dE):  #accept the flip with probability according to boltzlmann dist. 
                E1 = E2 #flip is acepted
                M = s.sum()   #calculate new magnetization 
            else:  
                s[i,j]*=-1  #if it wasn't accepted, revert back like this. 

        else:
            E1 = E2 #flip is accepted because energy falls
            M = s.sum()
        eplot.append(E1)
        Mplot.append(M)
    print("Done!")
    
    
    return s, eplot,Mplot


i = 0
Mplots = []
while i < 5: 
    s = make_lattice(20)
    snew,eplot,Mplot = metropolize(s,int(1e6))

                                       

    
    Mplots.append(Mplot)
    
    
    i +=1

plt.figure()
[plt.plot(Mplot) for Mplot in Mplots]
plt.xlabel('Steps',fontsize=25)
plt.ylabel('M',fontsize=25)
plt.savefig('../bin/magnetizations.png')



s = make_lattice(20)
#start with original
plt.figure()

plt.subplot(151)
plt.title('0 Steps')
plt.imshow(s)
plt.ylabel('$y$',fontsize = 25)



snew,eplot,Mplot = metropolize(s,int(1e3))
plt.subplot(152)
plt.title('1e3 Steps')
plt.imshow(snew)


plt.subplot(153)
snew,eplot,Mplot = metropolize(snew,int(1e4))
plt.title('1e4 Steps')
plt.imshow(snew)
plt.xlabel('$x$',fontsize=25)



plt.subplot(154)
snew,eplot,Mplot = metropolize(snew,int(1e5))
plt.title('1e5 Steps')
plt.imshow(snew)



plt.subplot(155)
snew,eplot,Mplot = metropolize(snew,int(1e6))
plt.title('1e6 Steps')
plt.imshow(snew)


plt.figure('../bin/ising_movie.png')
