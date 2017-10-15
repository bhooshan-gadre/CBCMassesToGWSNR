import numpy as np
import matplotlib.pyplot as plt

Msun = 1.989e30

# dat = np.genfromtxt('Data-for-flat-distri_M1_M2_M_chirpM_r_alpha_delta_iota_SNR.csv', delimiter=',')    # Reading data into an array
# dat = np.genfromtxt('Data-for-log-flat-distri_M1_M2_M_chirpM_r_alpha_delta_iota_SNR.csv', delimiter=',')    # Reading data into an array
dat = np.genfromtxt('Data-for-power-law-distri_M1_M2_M_chirpM_r_alpha_delta_iota_SNR.csv', delimiter=',')    # Reading data into an array


def Calc_for(M_l, M_up):
    i = 0
    
    # file1 = open('Data-for-flat-distri_M1_M2_M_chirpM_r_alpha_delta_iota_SNR.csv', 'r')        # Only for checking end of file
    # file1 = open('Data-for-log-flat-distri_M1_M2_M_chirpM_r_alpha_delta_iota_SNR.csv', 'r')        # Only for checking end of file
    file1 = open('Data-for-power-law-distri_M1_M2_M_chirpM_r_alpha_delta_iota_SNR.csv', 'r')        # Only for checking end of file
    
    SNR1 = []
    M1 = []
    lin = file1.readline()
    while True:

        lin = file1.readline()
        if not lin:
            break

        # if (dat[i, 2] > M_l*Msun) and (dat[i, 2] <= M_up*Msun):
        #     M1.append(dat[i, 2])
        #     SNR1.append(dat[i, 8])

        if (dat[i, 3] > M_l*Msun) and (dat[i, 3] <= M_up*Msun):
            M1.append(dat[i, 3])
            SNR1.append(dat[i, 8])        

        i += 1
    
    M1 = np.array(M1)
    wt_M = M1 ** ALPHA
    SNR1 = np.array(SNR1)

    mean1 = np.average(SNR1)
    med1 = np.median(SNR1)
    var1 = np.var(SNR1)

    n1, bin1 = np.histogram(SNR1, bins=np.linspace(5., 60., 56), weights=wt_M)

    return mean1, med1, var1, n1, bin1

def Plot(n1, bin1, M_l, M_up, style):
    # plt.plot(bin1[:-1], n1, style, label='%1.1f Msun < Total Mass <= %1.1f Msun\nMean: %1.1f, Median: %1.1f, Variance: %1.1f' %(M_l, M_up, mean1, med1, var1))
    plt.plot(bin1[:-1], n1, style, label='%1.1f Msun < Chirp Mass <= %1.1f Msun\nMean: %1.1f, Median: %1.1f, Variance: %1.1f' %(M_l, M_up, mean1, med1, var1))

def Normalise(Mat):
    total = 0.
    for ii in range(5):
        total += np.sum(Mat[ii])
    for ii in range(5):
        Mat[ii] = Mat[ii]/total
    return Mat[0], Mat[1], Mat[2], Mat[3], Mat[4]


# ALPHA = -2.3
ALPHA = 0.

# # for total M
# mean1, med1, var1, n1, bin1 = Calc_for(10., 28.)
# mean2, med2, var2, n2, bin2 = Calc_for(28., 46.)
# mean3, med3, var3, n3, bin3 = Calc_for(46., 64.)
# mean4, med4, var4, n4, bin4 = Calc_for(64., 82.)
# mean5, med5, var5, n5, bin5 = Calc_for(82., 100.)

# for chirpM
mean1, med1, var1, n1, bin1 = Calc_for(4., 12.2)
mean2, med2, var2, n2, bin2 = Calc_for(12.2, 20.4)
mean3, med3, var3, n3, bin3 = Calc_for(20.4, 28.6)
mean4, med4, var4, n4, bin4 = Calc_for(28.6, 36.8)
mean5, med5, var5, n5, bin5 = Calc_for(36.8, 45.)

n1, n2, n3, n4, n5 = Normalise([n1, n2, n3, n4, n5])

# frac8 = np.array([np.sum(n1[8:1000]), np.sum(n2[8:1000]), np.sum(n3[8:1000]), np.sum(n4[8:1000]), np.sum(n5[8:1000])])
frac9 = np.array([np.sum(n1[9:1000]), np.sum(n2[9:1000]), np.sum(n3[9:1000]), np.sum(n4[9:1000]), np.sum(n5[9:1000])])
frac13 = np.array([np.sum(n1[13:1000]), np.sum(n2[13:1000]), np.sum(n3[13:1000]), np.sum(n4[13:1000]), np.sum(n5[13:1000])])
frac23 = np.array([np.sum(n1[23:1000]), np.sum(n2[23:1000]), np.sum(n3[23:1000]), np.sum(n4[23:1000]), np.sum(n5[23:1000])])


# frac8 /= np.sum(frac8)
frac9 /= np.sum(frac9)
frac13 /= np.sum(frac13)
frac23 /= np.sum(frac23)


plt.figure()
# plt.plot(np.arange(10., 100., 18.), frac8, 'bo-',  label='SNR>8')

# plt.plot(np.arange(10., 100., 18.), frac9, 'co-',  label='SNR>9')
# plt.plot(np.arange(10., 100., 18.), frac13, 'go-',  label='SNR>13')
# plt.plot(np.arange(10., 100., 18.), frac23, 'mo-', label='SNR>23')

plt.plot(np.arange(4., 45., 8.2), frac9, 'co-',  label='SNR>9')
plt.plot(np.arange(4., 45., 8.2), frac13, 'go-',  label='SNR>13')
plt.plot(np.arange(4., 45., 8.2), frac23, 'mo-', label='SNR>23')


# # By Total Mass of Binary
# # Observed data
# plt.axvline(x=65., ymin=0., ymax = 1., linewidth=2, color='m', alpha=0.5)  #GW150914
# plt.axvline(x=36., ymin=0., ymax = 1., linewidth=2, color='c', alpha=0.5)  #LVT151012
# plt.axvline(x=21.7, ymin=0., ymax = 1., linewidth=2, color='g', alpha=0.5)  #GW151226
# plt.axvline(x=50.6, ymin=0., ymax = 1., linewidth=2, color='g', alpha=0.5)  #GW170104

# # Annotations
# plt.annotate(' GW150914\n (SNR 24) ', (65.+0.4, 0.035), fontsize=12, fontweight='bold', color='m')
# plt.annotate(' LVT151012\n (SNR 9.7) ', (36.+0.4, 0.035), fontsize=12, fontweight='bold', color='c')
# plt.annotate(' GW151226\n (SNR 13) ', (21.7+0.4, 0.035), fontsize=12, fontweight='bold', color='g')
# plt.annotate(' GW170104\n (SNR 13) ', (50.6+0.4, 0.035), fontsize=12, fontweight='bold', color='g')


# By Chirp Mass of Binary
# Observed data
plt.axvline(x=30., ymin=0., ymax = 1., linewidth=2, color='m', alpha=0.5)  #GW150914
plt.axvline(x=15.1, ymin=0., ymax = 1., linewidth=2, color='c', alpha=0.5)  #LVT151012
plt.axvline(x=8.9, ymin=0., ymax = 1., linewidth=2, color='g', alpha=0.5)  #GW151226
plt.axvline(x=21.1, ymin=0., ymax = 1., linewidth=2, color='g', alpha=0.5)  #GW170104

# Annotations
plt.annotate(' GW150914\n (SNR 24) ', (30.+0.4, 0.035), fontsize=12, fontweight='bold', color='m')
plt.annotate(' LVT151012\n (SNR 9.7) ', (15.1+0.4, 0.035), fontsize=12, fontweight='bold', color='c')
plt.annotate(' GW151226\n (SNR 13) ', (8.9+0.4, 0.035), fontsize=12, fontweight='bold', color='g')
plt.annotate(' GW170104\n (SNR 13) ', (21.1+0.4, 0.035), fontsize=12, fontweight='bold', color='g')


plt.ylabel('Probability of Detection (SNR>9, SNR>13 amd SNR>23)')

# plt.xlabel('Total Mass')
# plt.xticks(np.linspace(0., 100., 11))
# plt.title('Probability of Detection Vs. Total Mass (Distribution: weighted)')
# plt.title('Probability of Detection Vs. Total Mass (Distribution: log flat)')
# plt.title('Probability of Detection Vs. Total Mass (Distribution: Power Law)')

plt.xlabel('Chirp Mass')
plt.xticks(np.linspace(0., 45., 10))
# plt.title('Probability of Detection Vs. Chirp Mass (Distribution: weighted)')
# plt.title('Probability of Detection Vs. Chirp Mass (Distribution: log flat)')
plt.title('Probability of Detection Vs. Chirp Mass (Distribution: Power Law)')

plt.grid(True)
#plt.yticks(np.linspace(0., 0.5, 11))
#plt.axis([0., 100., 0., 0.7])
plt.legend(loc='upper right')

plt.show()

plt.figure()
# Plot(n1, bin1, 10., 28., 'ko-')
# Plot(n2, bin2, 28., 46., 'g^-')
# Plot(n3, bin3, 46., 64., 'rx-')
# Plot(n4, bin4, 64., 82., 'y*-')
# Plot(n5, bin5, 82., 100., 'b+-')

Plot(n1, bin1, 4., 12.2, 'ko-')
Plot(n2, bin2, 12.2, 20.4, 'g^-')
Plot(n3, bin3, 20.4, 28.6, 'rx-')
Plot(n4, bin4, 28.6, 36.8, 'y*-')
Plot(n5, bin5, 36.8, 45., 'b+-')

plt.xlabel('SNR')
plt.ylabel('Probability')

# plt.title('Probability Distribution of SNRs for different mass bins (Distribution: weighted)')
# plt.title('Probability Distribution of SNRs for different mass bins (Distribution: log flat)')
# plt.title('Probability Distribution of SNRs for different mass bins (Distribution: Power Law)')

# plt.title('Probability Distribution of SNRs for different chirp-mass bins (Distribution: weighted)')
# plt.title('Probability Distribution of SNRs for different chirp-mass bins (Distribution: log flat)')
plt.title('Probability Distribution of SNRs for different chirp-mass bins (Distribution: Power Law)')

plt.grid(True)
plt.xticks(np.linspace(0., 60., 31))
# plt.yticks(np.linspace(0, 0.3, 7))
# plt.axis([0, 60, 0, 0.3])
plt.legend(loc='upper right')

plt.show()
