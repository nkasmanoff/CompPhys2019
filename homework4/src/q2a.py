"""

Solutions and requested graphs for question 2a, studying the orbital decay of a 
binary black hole system. 


First, find the correct delta by integrating over ~10 orbits with very little energy loss



"""




import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')


def f(r,t):
    """
    For some given function with an input r and timestep t, we output the derivative of every value,
    which is in turn used in the rk4 integration 
    """
    
    x = r[0]
    y = r[1]
    vx = r[2]
    vy = r[3]
    
    Dx = vx
    Dy = vy
    
    R = np.sqrt(x**2 + y**2)
    
    Dvx = -G*M*x/(4*R**3)
    
    Dvy = -G*M*y/(4*R**3)
    
    
    return np.array([Dx,Dy,Dvx,Dvy])


def rk4(h,r,f,t):
    """
    RK4 Integration for some step size of h
    
    """
    
    k1 = h*f(r,t)
    k2 = h*f(r + 0.5*k1,t+0.5*h)
    k3 = h*f(r + 0.5*k2,t+0.5*h)
    k4 = h*f(r + k3,t+h)

    rout = r + (k1 + 2*k2 + 2*k3 + k4)/6 #update occurs here
    return rout


def find_delta(delta):
    
    h = 1e4 # This is the initial stepping size. 

    r0 = np.array([1.0,0.0,0.0,vy],float)

    r_sols = []
    ts = []
    t = 0 #initialize as 0 time 
    int_time = 45 #ends after about 10 orbits. 

    while t < int_time: 

        rtemp = rk4(h, r0, f,t)

        # Error at h and 2*h
        rError1 = rk4(h, rtemp, f,t)
        rError2 = rk4(2*h, r0, f,t)
        xerror = (rError1[0] - rError2[0])/30.  # 0th element -> x component
        yerror = (rError1[1] - rError2[1])/30.
        rho = h*delta/np.sqrt(xerror**2 + yerror**2)

        # If rho > 1, actual accuracy is better than the target accuracy. Keep it.
        if rho > 1:
          #  print("worked")
            t += h
            r0 = rtemp
            h = h*rho**(1/4) # Make it bigger since rho^1/4 > 1

            #save values and time step below
            r_sols.append(rtemp)
            ts.append(t)

        elif rho < 1:
         #   print('adapted')
            h = h * rho**(1/4)



    xs = [x[0] for x in r_sols]
    vxs = [x[2] for x in r_sols]


    ys = [y[1] for y in r_sols]
    vys = [y[3] for y in r_sols]

    rs = np.sqrt(np.power(xs,2) + np.power(ys,2))
    return xs,ys,ts,rs

G = M  = 1
rapo = 1
a = (1 + 1e-7)/2


vy = np.sqrt((G*M/4) *(2/rapo - 1/a))

a = (1 + 1e-7)/2


T = np.sqrt(16*np.pi**2 * a**3)

print("The time to take 10 orbits is approximately ", np.ceil(10*T))

print("The starting velocity at aphelion is ", vy)




deltas = [1e-4,1e-5,1e-6]
plt.figure(figsize=(8,8))
for d in deltas:
    xs, ys, ts, rs = find_delta(delta=d)
    
    plt.plot(ts,rs,linewidth=5,label = str(d))
    
plt.xlabel('$t$',fontsize=25)
plt.ylabel('$r$',fontsize=25)
plt.legend()
plt.ylim([0.9,1.1])
plt.savefig('../bin/deltas.png')