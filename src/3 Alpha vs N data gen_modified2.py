import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy import interpolate
import time

start_time = time.time()

print("Start Time: %s seconds" % (time.time() - start_time))

data = np.genfromtxt('/home/shreejit/asd.txt')      # Reading the amplitude spectral density
asd = interpolate.interp1d(data[:, 0], data[:, 1])  # Interpolating to obtain a continuous distribution

# Defining a function to integrate frequency dependent part in SNR calculation, from 9 Hz to ISCO freq
# input MTOT and output are numpy arrays
def sm(MTOT):
    f_isco = c ** 3. / (6. ** 1.5 * np.pi * G * MTOT)
    freq = np.zeros(Md * 1000).reshape(Md, 1000)
    
    for ii in range(Md):
        freq[ii] = np.linspace(9., f_isco[ii], 1000)

    df1 = np.transpose(freq)[1] - np.transpose(freq)[0]
    
    f1 = freq ** (-7. / 3.)
    y = asd(freq)
    sm1 = np.sum((f1 / y ** 2.), axis=1) * df1
    return sm1               # Array of interated values is returned

# Defining a function which takes (numpy arrays of) masses of BHs and gives (numpy array of) SNR
def find_SNR(M1, M2):
    M = M2 + M1
    mu = (M1 * M2) / M
    n = mu / M
    chirpM = M * np.power(n, (3. / 5.))
    Const = G * chirpM / np.power(c, 3.)
    wt_M = chirpM ** alpha

    #Calculation of SNR
    h1 = np.sqrt(5. / np.pi / 24.) * (Const * c)
    h1 /= r
    h1 *= np.power((Const * np.pi), -1. / 6.)
    h1 **= 2.

    SNR = np.sqrt(4. * h1 * sm(M))      # Integrated term 'sm(M)' multiplied
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
Ns = np.array([10., 11.])    #np.arange(5., 55., 1.)
alphas = np.array([-1., 0.])    #np.arange(-3., 3., 0.1)
Iters = 100

# Arrays
v = np.linspace(0., vol, Md)
r = v
r *= 3. / (4. * np.pi)
r += rl ** 3.
r **= (1. / 3.)

rf = open('Data_N_alpha_SNR8_9_10_11_12_13_det.csv', 'w')
wrt = csv.writer(rf)
wrt.writerow(['N','alpha', 'frac det>SNR8', 'frac det>SNR9', 'frac det>SNR10', 'frac det>SNR11', 'frac det>SNR12', 'frac det>SNR13', 'det>SNR13'])  # Last column will be detections above SNR 13. We have to compare it to the expected value 6 per year

counter = 1

for N in Ns:

    for alpha in alphas:

        for j in range(Iters):
            # timer = time.time()
            a = []
            M1 = np.random.uniform(5. * Msun, 50. * Msun, Md)
            M2 = np.random.uniform(5. * Msun, 50. * Msun, Md)

            SNR, wt_M = find_SNR(M1, M2)
                        
            n, bins = np.histogram(SNR, bins=np.linspace(0., 60., 61), weights=wt_M, density=True)

            # Calculating the fraction of area above a specific SNR value
            frac_8 = np.sum(n[8:1000])
            frac_9 = np.sum(n[9:1000])
            frac_10 = np.sum(n[10:1000])
            frac_11 = np.sum(n[11:1000])
            frac_12 = np.sum(n[12:1000])
            frac_13 = np.sum(n[13:1000])

            a = [N, alpha, frac_8, frac_9, frac_10, frac_11, frac_12, frac_13, frac_13*N]   # Array to be written as row
            wrt.writerow(a)
            print 'row: ', counter
            counter += 1
            # print time.time() - timer
        
rf.close()
print("End Time: %s seconds" % (time.time() - start_time))
