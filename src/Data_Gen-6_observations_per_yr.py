import numpy as np
import matplotlib.pyplot as plt
import csv


i = 0

dat = np.genfromtxt('Data_N_alpha_SNR8_9_10_11_12_13_det.csv', delimiter=',')

file1 = open('Data_N_alpha_SNR8_9_10_11_12_13_det.csv', 'r')
file2 = open('dat_6det_yr_N_alpha_SNR8_9_10_11_12_13.csv', 'w')
rfil = csv.writer(file2)

rfil.writerow(['N','alpha', 'frac det>SNR8', 'frac det>SNR9', 'frac det>SNR10', 'frac det>SNR11', 'frac det>SNR12', 'frac det>SNR13'])

AvN = []
lin = file1.readline()

while True:
	lin = file1.readline()
	if not lin:
		break
	
	if dat[i, 8] > 5.9 and dat[i, 8] < 6.1:
		rfil.writerow([dat[i, 0], dat[i, 1], dat[i, 2], dat[i, 3], dat[i, 4], dat[i, 5], dat[i, 6], dat[i, 7]])
	
	i += 1
