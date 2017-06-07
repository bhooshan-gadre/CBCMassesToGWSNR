from function_SNR import *
import time

start_time = time.time()
print("Start Time: %s seconds" % (time.time() - start_time))

# Input
rl = 648.8 * Mpc
ru = 2030.6 * Mpc
Md = 1000
Iter = 1000

vol = 4. / 3. * np.pi * (ru ** 3. - rl ** 3.)

# Arrays
r = np.linspace(0., vol, Md)
r *= 3. / (4. * np.pi)
r += rl ** 3.
r **= (1. / 3.)

rf = open('Data_M1_M2_M_r_alpha_delta_iota_SNR.csv', 'w')
title = 'BH1 Mass (Msun), BH2 Mass (Msun), Total Mass (Msun), Distance (Mpc), alpha, delta, iota, SNR'

data = np.zeros(Iter*Md*8).reshape(Iter*Md, 8)
index = 0
for count in range(Iter):
    M1 = np.random.uniform(5. * Msun, 50. * Msun, Md)
    M2 = np.random.uniform(5. * Msun, 50. * Msun, Md)
    alpha = np.random.uniform(0., 2.*np.pi, Md)
    delta = np.random.uniform(-np.pi/2., np.pi/2., Md)
    iota = np.random.uniform(0., np.pi, Md)
    M = M1 + M2
    Mratio = M1 / M2

    SNR = find_SNR(M1, M2, r, alpha, delta, iota)

    for count2 in range(Md):
        data[index] = np.array([M1[count2], M2[count2], M[count2], r[count2], alpha[count2], delta[count2], iota[count2], SNR[count2]])
        print index
        index += 1

# We sort the data such that the rows are sorted with ascending value of total masses of binaries
data = data[data[:,2].argsort()]

np.savetxt(rf, data, delimiter=',', header=title, newline='\n')
rf.close()
print("End Time: %s seconds" % (time.time() - start_time))
