from function_SNR import *
import csv
import matplotlib.pyplot as plt
import time

start_time1 = time.time()

# Input
rl = 648.8 * Mpc
ru = 2030.6 * Mpc
Md = 1000
vol = 4. / 3. * np.pi * (ru ** 3. - rl ** 3.)

# Arrays
r = np.linspace(0., vol, Md)
r *= 3. / (4. * np.pi)
r += rl ** 3.
r **= (1. / 3.)

# Last column will be detections above SNR 13. We have to compare it to the expected value 6 per year
title = 'N, ALPHA, frac det>SNR8, frac det>SNR9, frac det>SNR10, frac det>SNR11, frac det>SNR12, frac det>SNR13, det>SNR13'   
file1 = open('Data_N_ALPHA_SNR8_9_10_11_12_13_det.csv', 'w')

Ns = np.arange(5., 50., 1.)
ALPHAs = np.arange(-3., 1., 0.1)
Iters = 100

data = np.zeros([Iters*len(Ns)*len(ALPHAs), 9])
counter = 0
Tot_its = len(Ns) * len(ALPHAs) * Iters

for N in Ns:

      for ALPHA in ALPHAs:

            for j in range(Iters):

                  if counter < 1: timer = time.time()

                  timer2 = time.time()
                  M1 = np.random.uniform(5. * Msun, 50. * Msun, Md)
                  M2 = np.random.uniform(5. * Msun, 50. * Msun, Md)
                  alpha = np.random.uniform(0., 2.*np.pi, Md)
                  delta = np.random.uniform(-np.pi/2., np.pi/2., Md)
                  iota = np.random.uniform(0., np.pi, Md)
                  M = M1 + M2

                  wt_M = M ** ALPHA # This is the array defining weightages

                  SNR = find_SNR(M1, M2, r, alpha, delta, iota)

                  n, bins = np.histogram(SNR, bins=np.linspace(0., 60., 61), weights=wt_M, density=True)

                  # Calculating the fraction of area above a specific SNR value
                  data[counter][0] = N
                  data[counter][1] = ALPHA
                  data[counter][2] = np.sum(n[8:1000])
                  data[counter][3] = np.sum(n[9:1000])
                  data[counter][4] = np.sum(n[10:1000])
                  data[counter][5] = np.sum(n[11:1000])
                  data[counter][6] = np.sum(n[12:1000])
                  data[counter][7] = np.sum(n[13:1000])
                  data[counter][8] = data[j][7] * N

                  # print 'row: ', counter + 1
                  counter += 1
                  if counter < 2: it_time = time.time() - timer

                  if counter % 100 == 0:
                        print 'row: ', counter
                        print 'Time remaining: ', (Tot_its - counter) * it_time / 3600, ' Hrs'

np.savetxt(file1, data, delimiter=',', header=title, newline='\n')
file1.close()
print("Time taken for iterations and file saving: %s seconds" % (time.time() - start_time1))
start_time2 = time.time()

title2 = 'N, ALPHA, frac det>SNR8, frac det>SNR9, frac det>SNR10, frac det>SNR11, frac det>SNR12, frac det>SNR13'
file2 = open('dat_6det_yr_N_ALPHA_SNR8_9_10_11_12_13.csv', 'w')

data_6_det = data[np.where(data[:,8]>5.9)]
data_6_det = data_6_det[np.where(data_6_det[:,8]<6.1)]

np.savetxt(file2, data_6_det, delimiter=',', header=title2, newline='\n')
file2.close()
print("Time taken for next file saving: %s seconds" % (time.time() - start_time2))
