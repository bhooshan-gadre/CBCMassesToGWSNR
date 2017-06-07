from function_SNR import *
import csv
import matplotlib.pyplot as plt
import time

start_time = time.time()

print("Start Time: %s seconds" % (time.time() - start_time))

# Input
rl =  * Mpc
ru =  * Mpc
Md = 1000
vol = 4. / 3. * np.pi * (ru ** 3. - rl ** 3.)

# Arrays
r = np.linspace(0., vol, Md)
r *= 3. / (4. * np.pi)
r += rl ** 3.
r **= (1. / 3.)

# Last column will be detections above SNR 13. We have to compare it to the expected value 6 per year
title = 'N, gamma, frac det>SNR8, frac det>SNR9, frac det>SNR10, frac det>SNR11, frac det>SNR12, frac det>SNR13, det>SNR13'   
rf = open('Data_N_gamma_SNR8_9_10_11_12_13_det.csv', 'w')   

Ns = np.array([10., 11.])    #np.arange(5., 55., 1.)
gammas = np.array([-1., 0.])    #np.arange(-3., 3., 0.1)
Iters = 100

a = np.zeros([Iters*len(Ns)*len(gammas), 9])
counter = 0

for N in Ns:

      for gamma in gammas:

            for j in range(Iters):
                  # timer2 = time.time()
                  M1 = np.random.uniform(5. * Msun, 50. * Msun, Md)
                  M2 = np.random.uniform(5. * Msun, 50. * Msun, Md)
                  alpha = np.random.uniform(0., 2.*np.pi, Md)
                  delta = np.random.uniform(-np.pi/2., np.pi/2., Md)
                  iota = np.random.uniform(0., np.pi, Md)
                  M = M1 + M2

                  wt_M = M ** gamma # This is the array defining weightages

                  SNR = find_SNR(M1, M2, r, alpha, delta, iota)

                  n, bins = np.histogram(SNR, bins=np.linspace(0., 60., 61), weights=wt_M, density=True)

                  # Calculating the fraction of area above a specific SNR value
                  a[counter][0] = N
                  a[counter][1] = gamma
                  a[counter][2] = np.sum(n[8:1000])
                  a[counter][3] = np.sum(n[9:1000])
                  a[counter][4] = np.sum(n[10:1000])
                  a[counter][5] = np.sum(n[11:1000])
                  a[counter][6] = np.sum(n[12:1000])
                  a[counter][7] = np.sum(n[13:1000])
                  a[counter][8] = a[j][7] * N

                  print 'row: ', counter + 1
                  counter += 1
                  # print 'Time: ', time.time() - timer2

np.savetxt(rf, a, delimiter=',', header=title, newline='\n')
rf.close()
print("End Time: %s seconds" % (time.time() - start_time))
