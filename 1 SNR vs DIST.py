import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

data = np.genfromtxt('/home/shreejit/asd.txt')		# Reading the amplitude spectral density
asd = interpolate.interp1d(data[:, 0], data[:, 1])	# Interpolating to obtain a continuous distribution

# Defining a function to integrate frequency dependent part in SNR calculation, from 9 Hz to ISCO freq
def sm(MTOT):
	f_isco = c ** 3. / (6. ** 1.5 * np.pi * G * MTOT)
	freq = np.arange(9., f_isco, df)
	f1 = freq ** (-7. / 3.)
	y = asd(freq)
	sm1 = np.sum(f1 / y ** 2.) * df
	return sm1
	
# Defining a function which takes masses of BHs and gives SNR
def find_SNR(M1, M2):
	M = M2 + M1
	mu = (M1 * M2) / M
	n = mu / M
	chirpM = M * np.power(n, (3. / 5.))
	Const = G * chirpM / np.power(c, 3.)

	#Calculation of SNR
	h1 = np.sqrt(5. / np.pi / 24.) * (Const * c)
	h1 /= r
	h1 *= np.power((Const * np.pi), -1. / 6.)
	h1 **= 2.

	SNR = np.sqrt(4. * h1 * sm(M))	# Integrated term 'sm(M)' multiplied
	return SNR
	
# This function finds out the value of distance at which we get a specific value of SNR (5 or 60 in our case)
# Gets input in the form (tolerance in SNR, value of SNR at which r is to be found out, array of r, array of SNR)
def r_for_SNR(tol, val, x, y):
	yy = y[np.where(y>(val-tol))]
	# print np.where(yy<(val+tol))
	return x[np.where(yy<(val+tol))]
	
# Main function which finds SNR, cutoff on distance and plots the graphs
# Takes input as a matrix containing info: [SNRcutoff, M_BH, low_lim_x_axis, uppr_lim_x_axis, uppr_lim_y_axis, tolerance, offset for annotation in x, in y, annotation]
def calc_plot(Mat):
	cut = Mat[0]
	ML = Mat[1]
	x_lL = Mat[2]
	x_L = Mat[3]
	y_L = Mat[4]
	tolrnc = Mat[5]

	M1 = ML * Msun
	M2 = ML * Msun

	SNR = find_SNR(M1, M2)

	rL = r_for_SNR(tolrnc, cut, r, SNR)/Mpc
	print '(', rL, ', ', cut, ')'

	# Plotting
	fig = plt.figure()
	sb = fig.add_subplot(111)
	sb2 = fig.add_subplot(111)

	sb.plot(r/Mpc, SNR, '-')
	sb2.plot([rL], [cut], 'o')
	sb2.annotate('  (%1.1f Mpc,  %1.1f)' %(rL, cut), (rL + Mat[6], cut + Mat[7]))
	sb.grid(True)
	plt.xlabel('Distance (Mpc)')
	plt.ylabel('SNR')
	plt.title(Mat[8])
	plt.axis([x_lL, x_L, 0., y_L])


# Constants
Msun = 1.989e30
G = 6.67259e-11 
c = 2.99792458e8 
Mpc = 3.08568025e22
df = 0.01
i = 0

# Defining numpy array of r
r = np.linspace(100, 5000, 1001)
r *= Mpc

Mats = [[60., 50., 0., 3000., 300., 0.2, 50., 4., 'Lower Cutoff'], [5., 5., 0., 5000., 100., 0.005, -700., 4., 'Higher Cutoff']]

calc_plot(Mats[0])
calc_plot(Mats[1])

plt.show()

