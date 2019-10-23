import matplotlib.pyplot as plt
import numpy as np

blur = np.loadtxt('../dat/blur.txt')

plt.figure()
plt.imshow(blur,cmap='gray')
plt.savefig('../bin/blur.png')


def gaussian_2D(x,y,sigma):
    """
    Returns a 2d gaussian of the form f(x,y) = e^-(x2+y2 / 2s2
    
    Center this at x0 y0 by shifting 

    """
    
    from numpy import exp,power
    
    return exp(-(power(x,2)+power(y,2))/(2*power(sigma,2)))
    

xs = np.arange(0,1024,1)#np.concatenate([np.arange(0,512,1),np.arange(0,512,1)[::-1]]) #reversed ? 
ys = np.arange(0,1024,1)#np.concatenate([np.arange(0,512,1),np.arange(0,512,1)[::-1]]) #1 to 50 for each

psf_val = []
n = 1024
for y in ys: 
    for x in xs:
        if x > n/2:
            x = n - x
        if y > n/2:
            y = n - y
    
        psf_val.append(gaussian_2D(x,y,sigma=25))
    
psf_val = np.array(psf_val)


psf = psf_val.reshape(1024,1024)
plt.figure()
plt.imshow(psf,cmap='gray',vmin = 0.01)
plt.savefig('../bin/psf.png')

# get ft with 

from numpy.fft import rfft2 , irfft2


blur_fft = rfft2(blur)
psf_fft = rfft2(psf)


#remove all small points. 

iz = np.where(psf_fft < 1e-3)
psf_fft[iz] = 1

plt.figure()
plt.imshow(irfft2(np.divide(blur_fft,psf_fft)/1024**2),'gray')
plt.savefig('../bin/deblur.png')


