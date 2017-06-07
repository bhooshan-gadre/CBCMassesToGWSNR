import numpy as np
import matplotlib.pyplot as plt

Msun = 1.989e30

dat = np.genfromtxt('Data_r_M_chirpM_Mratio_SNR.csv', delimiter=',')    # Reading data into an array

def Calc_for(M_l, M_up):
    i = 1
    file1 = open('Data_r_M_chirpM_Mratio_SNR.csv', 'r')        # Only for checking end of file
    SNR1 = []
    M1 = []
    lin = file1.readline()
    while True:

        lin = file1.readline()
        if not lin:
            break

        if (dat[i, 1] > M_l) and (dat[i, 1] <= M_up):
            M1.append(dat[i, 1])
            SNR1.append(dat[i, 4])

        i += 1
    
    M1 = np.array(M1)
    wt_M = M1 ** gamma
    SNR1 = np.array(SNR1)

    mean1 = np.average(SNR1)
    med1 = np.median(SNR1)
    var1 = np.var(SNR1)

    n1, bin1 = np.histogram(SNR1, bins=np.linspace(5., 60., 56), weights=wt_M)

    return mean1, med1, var1, n1, bin1

def Plot(n1, bin1, M_l, M_up, style):
    plt.plot(bin1[:-1], n1, style, label='%1.1f Msun < Total Mass <= %1.1f Msun\nMean: %1.1f, Median: %1.1f, Variance: %1.1f' %(M_l, M_up, mean1, med1, var1))

def Normalise(Mat):
    total = 0.
    for ii in range(5):
        total += np.sum(Mat[ii])
    for ii in range(5):
        Mat[ii] = Mat[ii]/total
    return Mat[0], Mat[1], Mat[2], Mat[3], Mat[4]


gamma = -5./6.

mean1, med1, var1, n1, bin1 = Calc_for(10., 28.)
mean2, med2, var2, n2, bin2 = Calc_for(28., 46.)
mean3, med3, var3, n3, bin3 = Calc_for(46., 64.)
mean4, med4, var4, n4, bin4 = Calc_for(64., 82.)
mean5, med5, var5, n5, bin5 = Calc_for(82., 100.)

n1, n2, n3, n4, n5 = Normalise([n1, n2, n3, n4, n5])


plt.figure()
Plot(n1, bin1, 10., 28., 'ko-')
Plot(n2, bin2, 28., 46., 'g^-')
Plot(n3, bin3, 46., 64., 'rx-')
Plot(n4, bin4, 64., 82., 'y*-')
Plot(n5, bin5, 82., 100., 'b+-')

plt.xlabel('SNR')
plt.ylabel('Probability')
plt.title('Probability Distribution of SNRs for different mass bins')
plt.grid(True)
plt.xticks(np.linspace(0., 60., 31))
# plt.yticks(np.linspace(0, 0.3, 7))
# plt.axis([0, 60, 0, 0.3])
plt.legend(loc='upper right')

plt.show()
