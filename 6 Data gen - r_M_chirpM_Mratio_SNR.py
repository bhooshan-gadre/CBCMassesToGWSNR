import numpy as np
import csv
from scipy import interpolate
import time

start_time = time.time()

print("Start Time: %s seconds" % (time.time() - start_time))

data = np.genfromtxt('/home/shreejit/asd.txt')
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
    return sm1               # Array of integrated values is returned

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

    SNR = np.sqrt(4. * h1 * sm(M))      # Integrated term 'sm(M)' multiplied
    return SNR, chirpM, M

#######################################################################


# Constants
Msun = 1.989e30
G = 6.67259e-11 
c = 2.99792458e8 
Mpc = 3.08568025e22

# Input
rl = 648.8 * Mpc
ru = 2030.6 * Mpc
Md = 1000
df = 0.01
count = 1

vol = 4. / 3. * np.pi * (ru ** 3. - rl ** 3.)

# Arrays
v = np.linspace(0., vol, Md)
r = v
r *= 3. / (4. * np.pi)
r += rl ** 3.
r **= (1. / 3.)

rf = open('Data_r_M_chirpM_Mratio_SNR.csv', 'w')
wrt = csv.writer(rf)
wrt.writerow(['Distance (Mpc)', 'Total Mass (Msun)', 'Chirp Mass (Msun)', 'Mass Ratio', 'SNR'])

for count in range(100):
    M1 = np.random.uniform(5. * Msun, 50. * Msun, Md)
    M2 = np.random.uniform(5. * Msun, 50. * Msun, Md)

    Mratio = M1 / M2
    SNR, chirpM, M = find_SNR(M1, M2)

    for count2 in range(Md):
    	# Getting Mass Ratio in range 0 - 1.
        if Mratio[count2] > 1.: Mratio2 = pow(Mratio[count2], -1.)
        else: Mratio2 = Mratio[count2]

        wrt.writerow([r[count2]/Mpc, M[count2]/Msun, chirpM[count2]/Msun, Mratio2, SNR[count2]])
        # print count
        # count += 1

rf.close()
print("End Time: %s seconds" % (time.time() - start_time))
