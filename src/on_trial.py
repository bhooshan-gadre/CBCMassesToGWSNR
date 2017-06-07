import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

a = np.arange(1., 5., 1.)
b = np.arange(10., 50., 10.)

x, y = np.meshgrid(a, b)
x = x.reshape(len(x)**2)
y = y.reshape(len(y)**2)
c = x ** 2 + y ** 3

inter = interpolate.interp2d(x, y, c)


sky_lookup = np.genfromtxt('/home/shreejit/H1_allSky_antenna_patterns.dat')
#find_F_plus = interpolate.interp2d(sky_lookup[:,0], sky_lookup[:,1], sky_lookup[:,2])	# Interpolating the right ascension and 
#find_F_cross = interpolate.interp2d(sky_lookup[:,0], sky_lookup[:,1], sky_lookup[:,3])	# declinations with F_plus and F_cross
