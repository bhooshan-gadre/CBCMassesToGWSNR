import numpy as np
import matplotlib.pyplot as plt


i = 0

dat = np.genfromtxt('dat_6det_yr_N_ALPHA_SNR8_9_10_11_12_13.csv', delimiter=',')

file1 = open('dat_6det_yr_N_ALPHA_SNR8_9_10_11_12_13.csv', 'r')

AvN = []
lin = file1.readline()
while True:
	
	lin = file1.readline()
	if not lin:
		break
	
	AvN.append([dat[i+1, 0], dat[i+1, 1]])
	
	i += 1

AvN = np.array(AvN)

X = np.linspace(3, 22, 19)
Y = np.arange(-3.2, 3.4, 0.1)

H, X, Y = np.histogram2d(AvN[:,0], AvN[:,1], bins=(X, Y))
H = H.T

fig = plt.figure(figsize=(40,60))
# plt.plot(AvN[:,0],AvN[:,1], 'bo')

ax = fig.add_subplot(111, title='pcolormesh: actual edges') #, aspect='equal')

Xm, Ym = np.meshgrid(X, Y)
ax.pcolormesh(Xm, Ym, H)

cax = ax.imshow(H, interpolation='none')
fig.colorbar(cax, ax=ax)

plt.xticks(np.linspace(4, 22, 10))
plt.yticks(np.arange(-3.2, 3.4, 0.4))
plt.axis([4, 22, -3.2, 3.2])
plt.xlabel('No of injections (N)')
plt.ylabel('ALPHA')
plt.title('ALPHA vs N (for constraint of 6 detections/yr having SNR>13)')

plt.show()
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
