from astropy.cosmology import LambdaCDM
import numpy as np
import configparser
import sys

config = configparser.ConfigParser()
config.read(sys.argv[1])

frest = 1.42e9
f1 = 1.15e9
f2 = 0.95e9
df = 30e3
z1a = (frest/f1)-1
z1b = (frest/(f1-df))-1
z2a = (frest/f2)-1
z2b = (frest/(f2-df))-1
print (z1a, z2a)
pix_size = 2.5/206265



H = config.getfloat('cosmology', 'H')
M = config.getfloat('cosmology', 'M')
L = config.getfloat('cosmology', 'L')
c = config.getfloat('cosmology', 'c')
G = config.getfloat('cosmology', 'G')

cosmo = LambdaCDM(H0 = H, Om0 = M, Ode0 = L)

dDA1 = cosmo.angular_diameter_distance(z1a).value-cosmo.angular_diameter_distance(z1b).value # Mpc
dDA2 = cosmo.angular_diameter_distance(z2a).value-cosmo.angular_diameter_distance(z2b).value # Mpc



dRADec1 = cosmo.angular_diameter_distance(z1a).value*pix_size # Mpc

dRADec2 = cosmo.angular_diameter_distance(z2a).value*pix_size # Mpc
print (dDA1, dDA2)
print (dRADec1, dRADec2)

