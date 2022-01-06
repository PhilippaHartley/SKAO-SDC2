from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np


import numpy as np


import time
import numpy as np
import os, glob, sys
from astropy.io import fits
import galsim
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm


cube = fits.open('out/sky_HI_1050_v6.fits')

data = cube[0].data[:,3000:4000,: ]


mu = 0
sigma = 5e-6

'''
noisedata = data#[0:1000, 350:350+512,350:350+512]

for i in range(len(noisedata)):
    print (noisedata[i].shape)

    noise = np.random.normal(mu, sigma, noisedata[i].shape)
    noisedata[i]+=noise

#summed_line_flux = np.sum(noisedata, axis = 0) 
#plt.imshow(summed_line_flux)
#plt.show()


#cube[0].data = noisedata



#cube.writeto('sky_HI_1150tests_noise_small.fits', overwrite = True)
data = noisedata
'''


front =np.sum(data, axis = 0)
top = np.sum(data, axis = 1)
side = np.sum(data, axis = 2)



hdu = fits.PrimaryHDU(front)
    
hdulist = fits.HDUList([hdu])






hdulist.writeto('front.fits' ,overwrite=True)


hdu = fits.PrimaryHDU(top)
    
hdulist = fits.HDUList([hdu])






hdulist.writeto('top.fits' ,overwrite=True)

hdu = fits.PrimaryHDU(side)
    
hdulist = fits.HDUList([hdu])






hdulist.writeto('side.fits' ,overwrite=True)






exit()

np.putmask(front, front==0, front[front>0].min())

np.putmask(top, top==0, top[top>0].min())

np.putmask(side, side==0, side[side>0].min())


plt.close('all')
fig = plt.figure()
ax = fig.gca(projection='3d')

X = np.linspace(-5, 5, 43)
Y = np.linspace(-5, 5, 28)
X, Y = np.meshgrid(X, Y)

varone=np.random.rand(75,28,43) * 5.0 - 10.0
Z=varone[0,:,:]

cset = [[],[],[]]

# this is the example that worked for you:
cset[0] = ax.contourf(front, top, size, zdir='z', offset=5,
                      levels=np.linspace(np.min(Z),np.max(Z),30),cmap='jet')

# now, for the x-constant face, assign the contour to the x-plot-variable:
cset[1] = ax.contourf(Z, Y, X, zdir='x', offset=5,
                      levels=np.linspace(np.min(Z),np.max(Z),30),cmap='jet')

# likewise, for the y-constant face, assign the contour to the y-plot-variable:
cset[2] = ax.contourf(X, Z, Y, zdir='y', offset=-5,
                      levels=np.linspace(np.min(Z),np.max(Z),30),cmap='jet')

# setting 3D-axis-limits:    
ax.set_xlim3d(-5,5)
ax.set_ylim3d(-5,5)
ax.set_zlim3d(-5,5)

plt.show()
