
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import interpolate
sns.set_style('whitegrid')
sns.set_palette('husl')


Pk = pd.read_csv('../dat/lcdm_z0.matter_pk',delimiter=' ',names=['k','P(k)','idk1',
                                                                'idk2'],header=None)

Pk = Pk[Pk.columns[0:2]]

def integrate(f,a,b,N):

    h = (b - a) / (N)
        
    s = f(a) + f(b) + 4*f(b-h)
    for k in range(1,N//2):
        s += 4*f(a + (2*k-1)*h) + 2*f(a+2*k*h)
    I_simpson = (h/3)*s  #   expected error ~h^4
        

    return I_simpson


#cubic spline pk, then multiply 

corrs = []
rs = np.linspace(2,120,60)
for r in rs:
    cs = interpolate.CubicSpline(x=Pk['k'],y=Pk['P(k)'])
    f3 = lambda k: ((k*k*cs(k)*np.sin(k*r)) / (k*r))
    integral = integrate(f = f3,a = min(Pk['k'])+.00001, b = max(Pk['k']) - .00001 , N=int(1e5)) / (2*np.pi**2)#THIS WAS ITTT!!
    corrs.append(integral)
    
  #  print(r,integral)


plt.figure(figsize=(10,10))

corrs = np.array(corrs)
corrs80 = corrs[rs>80]
rs80 = rs[rs>80]
plt.plot(rs,corrs,'b--')
plt.plot(rs80[corrs80.argmax()],corrs80[corrs80.argmax()],'y*',markersize = 20,label = 'BAO Bump')
plt.xlabel("Separation r [Mpc/h]",fontsize = 15)
plt.ylabel('Correlation Function ξ (r)',fontsize = 15)
plt.title("Correlation Function vs. Separation ",fontsize = 20)
plt.legend()

plt.savefig('../bin/q3plot1.png')



plt.figure(figsize=(10,10))
r2corrs = np.multiply(np.power(rs,2),corrs)
r2corrs80 = r2corrs[rs > 80]
plt.plot(rs,r2corrs,'b--')
plt.plot(rs80[r2corrs80.argmax()],r2corrs80[r2corrs80.argmax()],'y*',markersize=20,label = 'BAO Bump')
plt.xlabel("Separation r [Mpc/h]",fontsize = 15)
plt.ylabel('r^2 * ξ (r)',fontsize = 15)
plt.title("r^2 * ξ (r) vs. r  ",fontsize = 20)
plt.legend()

plt.savefig('../bin/q3plot2.png')
