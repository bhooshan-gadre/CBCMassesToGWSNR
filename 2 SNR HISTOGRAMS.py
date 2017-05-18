import numpy as np
import random
import matplotlib.pyplot as plt
from scipy import interpolate

# Constants
Msun = 1.989e30
G = 6.67259e-11 
c = 2.99792458e8 
Mpc = 3.08568025e22
df = 0.01

data = np.genfromtxt('/home/shreejit/asd.txt')		# Reading the amplitude spectral density
asd = interpolate.interp1d(data[:, 0], data[:, 1])	# Interpolating to obtain a continuous distribution

# Defining a function to integrate frequency dependent part in SNR calculation, from 9 Hz to ISCO freq
# input MTOT and output are numpy arrays
def sm(MTOT):
	fr_sum = np.zeros(len(MTOT))		# making an array of zeroes of dim of MTOT
	for ii in range(len(MTOT)):		# This loop replaces each zero with respective integration value (integrated from 9Hz to ISCO freq)
		f_isco = c ** 3. / (6. ** 1.5 * np.pi * G * MTOT[ii])
		freq = np.arange(9., f_isco, df)
		f1 = freq ** (-7. / 3.)
		y = asd(freq)
		sm1 = np.sum(f1 / y ** 2.)
		fr_sum[ii] = sm1 * df
	return fr_sum				# Array of interated values is returned

# Defining a function which takes (numpy arrays of) masses of BHs and gives (numpy array of) SNR
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

	SNR = np.sqrt(4. * h1 * sm(M))		# Integrated term 'sm(M)' multiplied
	return SNR


# Input
# Lower and upper limits on space
rl = 648.8 * Mpc
ru = 2030.6 * Mpc
Md = 1000
M1 = np.random.uniform(5. * Msun, 50. * Msun, Md)	# 'Md' no of injections (of random masses in 5 - 50 Msun) made in the vol bound by 'rl' and 'ru'
M2 = np.random.uniform(5. * Msun, 50. * Msun, Md)

data = np.genfromtxt('/home/shreejit/asd.txt')		# Reading the amplitude spectral density
asd = interpolate.interp1d(data[:, 0], data[:, 1])	# Interpolating to obtain a continuous distribution

vol = 4. / 3. * np.pi * (ru ** 3. - rl ** 3.)		# Vol of visible space

# Descetizing 'vol' in equal parts
v = np.linspace(0., vol, Md)
# Bringing the 'vol' descretization into 'r'
r = v
r *= 3. / (4. * np.pi)
r += rl ** 3.
r **= (1. / 3.)

SNR = find_SNR(M1, M2)

# Inj = len(SNR)

# Plotting
plt.figure()		# Histogram of SNR
n, bins, patches = plt.hist(SNR, 50, color='g', normed=1, alpha=0.80)
plt.xlabel('SNR')
plt.ylabel('Likelihood of BBH Signals')
plt.title('Likelihood of BBH Signals Vs. SNR')
plt.grid(True)

plt.figure()		# Cumulative Histogram of SNR
plt.hist(SNR, bins=50, color='b', normed=1, cumulative=-1, alpha=0.80)
plt.xlabel('SNR')
plt.ylabel('Cumulative no of BBH Signals')
plt.title('Cumulative Number of BBH Signals Vs. SNR')
plt.grid(True)
plt.axis([0., 65., 0., 1.0])
plt.show()

