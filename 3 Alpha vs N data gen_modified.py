import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy import interpolate
import time

start_time = time.time()

print("Start Time: %s seconds" % (time.time() - start_time))

data = np.genfromtxt('/home/shreejit/asd.txt')		# Reading the amplitude spectral density
asd = interpolate.interp1d(data[:, 0], data[:, 1])	# Interpolating to obtain a continuous distribution

# Defining a function which takes (numpy arrays of) masses of BHs and gives (numpy array of) SNR
def find_SNR(M1, M2):
	M = M2 + M1
	mu = (M1 * M2) / M
	n = mu / M
	chirpM = M * np.power(n, (3. / 5.))
	Const = G * chirpM / np.power(c, 3.)
	wt_M = chirpM ** alpha

	f_isco = c ** 3. / (6. ** 1.5 * np.pi * G * M)
	sm1 = np.ones(Md)

	timer = time.time()
	for ii in range(Md):
		sm1[ii] = np.sum(integrand[np.where(freq < f_isco[ii])])
	print time.time() - timer

	#Calculation of SNR
	h1 = np.sqrt(5. / np.pi / 24.) * (Const * c)
	h1 /= r
	h1 *= np.power((Const * np.pi), -1. / 6.)
	h1 **= 2.

	SNR = np.sqrt(4. * h1 * sm1)		# Integrated term 'sm(M)' multiplied
	return SNR, wt_M

# Constants
Msun = 1.989e30
G = 6.67259e-11 
c = 2.99792458e8 
Mpc = 3.08568025e22
df = 0.01

# Input
rl = 648.8 * Mpc
ru = 2030.6 * Mpc
Md = 1000
vol = 4. / 3. * np.pi * (ru ** 3. - rl ** 3.)

# Ns = np.random.poisson(lam = 17.6, size = 2000)
Ns = np.array([10.])    #np.arange(5., 55., 1.)
alphas = np.array([-1.])    #np.arange(-3., 3., 0.1)
Iters = 20

# Arrays
v = np.linspace(0., vol, Md)
r = v
r *= 3. / (4. * np.pi)
r += rl ** 3.
r **= (1. / 3.)

title = np.array(['N','alpha', 'frac det>SNR8', 'frac det>SNR9', 'frac det>SNR10', 'frac det>SNR11', 'frac det>SNR12', 'frac det>SNR13', 'det>SNR13'])
rf = open('Data_N_alpha_SNR8_9_10_11_12_13_det.csv', 'a')
# np.savetxt(rf, np.transpose(title), fmt='%s', delimiter=',') 	# Last column will be detections above SNR 13. We have to compare it to the expected value 6 per year

freq = np.arange(9., 8000., df)
y = asd(freq)
integrand = freq ** (-7. / 3.) * df / y ** 2.

counter = 1

for N in Ns:

    for alpha in alphas:

        a = np.ones([Iters, 9])

        for j in range(Iters):
            # timer2 = time.time()
            M1 = np.random.uniform(5. * Msun, 50. * Msun, Md)
            M2 = np.random.uniform(5. * Msun, 50. * Msun, Md)
            SNR, wt_M = find_SNR(M1, M2)
            
            n, bins = np.histogram(SNR, bins=np.linspace(0., 60., 61), weights=wt_M, density=True)

            # Calculating the fraction of area above a specific SNR value
            a[j][0] = N
            a[j][1] = alpha
            a[j][2] = np.sum(n[8:1000])
            a[j][3] = np.sum(n[9:1000])
            a[j][4] = np.sum(n[10:1000])
            a[j][5] = np.sum(n[11:1000])
            a[j][6] = np.sum(n[12:1000])
            a[j][7] = np.sum(n[13:1000])
            a[j][8] = a[j][7] * N
            
            print 'row: ', counter
            counter += 1
            # print 'Time: ', time.time() - timer2
		
        
        np.savetxt(rf, a, delimiter=',')
        
rf.close()
print("End Time: %s seconds" % (time.time() - start_time))
