import numpy as np
import matplotlib.pyplot as plt

# # Collect SNR values in bins of total/chirp mass
# def Calc_for(M_l, M_up, mtyp):
#     # Empty arrays of SNR and M_t (total/chirp mass). Will collect relevant values for which conditions are satisfied
#     SNR1 = []
#     M_t = []
#     for count in range(len(dat[:,0])):
#         if (dat[count, Dict[mtyp]+2] > M_l*Msun) and (dat[count, Dict[mtyp]+2] <= M_up*Msun):
#             M_t.append(dat[count, Dict[mtyp]+2])
#             SNR1.append(dat[count, 8])
#     M_t = np.array(M_t)
#     SNR1 = np.array(SNR1)
#     # Weighing array for uniform distribution as M_t^-ALPHA
#     # ALPHA = 0 for other distributions, ie no Weighing
#     wt_M = M_t ** -ALPHA
#     # Calculating maen, median and variance of SNR array
#     mean1 = np.average(SNR1)
#     med1 = np.median(SNR1)
#     var1 = np.var(SNR1)
#     # Getting n1 as an array of histogram weightages
#     n1, bin1 = np.histogram(SNR1, bins=np.linspace(5., 60., 56), weights=wt_M)
#     return mean1, med1, var1, n1, bin1
# # Plot histograms for different mass bins
# def Plot(n1, bin1, M_l, M_up, style):
#     plt.plot(bin1[:-1], n1, style, label='%1.1f Msun < %s <= %1.1f Msun\nMean: %1.1f, Median: %1.1f, Variance: %1.1f' %(M_l, mtyp, M_up, mean1, med1, var1))
# # Normalise all the histograms collectively
# def Normalise(Mat):
#     total = 0.
#     for ii in range(5):
#         total += np.sum(Mat[ii])
#     for ii in range(5):
#         Mat[ii] = Mat[ii]/total
#     return Mat[0], Mat[1], Mat[2], Mat[3], Mat[4]

# # Mass of Sun
# Msun = 1.989e30
# # Get the type of distribution as an input (Exact name to be input in str format)
# inp_distri = input('Which distribution to generate plots for? (Uniform, log_flat, pow_law)\n: ')
# # Plot against Total Mass / Chirp Mass (Exact name to be input in str format)
# mtyp = input('Which mass to plot against? (Total Mass / Chirp Mass)\n: ')
# # Dictionary to identify the input mtyp
# Dict = {'Total Mass': 0, 'Chirp Mass': 1}
# # Reading resp data file into an array
# dat = np.genfromtxt('./../Data/Data-for-%s-distri_M1_M2_M_chirpM_r_alpha_delta_iota_SNR.csv' %inp_distri, delimiter=',')
# # Weighing factor ALPHA is 2.3 for Uniform distribution and zero for others
# ALPHA = 0.
# if inp_distri == 'Uniform': ALPHA = 2.3
# # Defining bin edges
# m_lims = np.zeros(12).reshape(2,6)
# # In Total Mass
# m_lims[0,:] = np.linspace(10., 100., 6)
# # In Chirp Mass
# m_lims[1,:] = np.linspace(4., 45., 6)
# # Calculating various quantities for all the bins
# mean1, med1, var1, n1, bin1 = Calc_for(m_lims[Dict[mtyp], 0], m_lims[Dict[mtyp], 1], mtyp)
# mean2, med2, var2, n2, bin2 = Calc_for(m_lims[Dict[mtyp], 1], m_lims[Dict[mtyp], 2], mtyp)
# mean3, med3, var3, n3, bin3 = Calc_for(m_lims[Dict[mtyp], 2], m_lims[Dict[mtyp], 3], mtyp)
# mean4, med4, var4, n4, bin4 = Calc_for(m_lims[Dict[mtyp], 3], m_lims[Dict[mtyp], 4], mtyp)
# mean5, med5, var5, n5, bin5 = Calc_for(m_lims[Dict[mtyp], 4], m_lims[Dict[mtyp], 5], mtyp)
# # Normalising the histograms
# n1, n2, n3, n4, n5 = Normalise([n1, n2, n3, n4, n5])

# # Normalised probability of obtaining above SNR values 
# frac8 = np.array([np.sum(n1[8:1000]), np.sum(n2[8:1000]), np.sum(n3[8:1000]), np.sum(n4[8:1000]), np.sum(n5[8:1000])])
# frac9 = np.array([np.sum(n1[9:1000]), np.sum(n2[9:1000]), np.sum(n3[9:1000]), np.sum(n4[9:1000]), np.sum(n5[9:1000])])
# frac13 = np.array([np.sum(n1[13:1000]), np.sum(n2[13:1000]), np.sum(n3[13:1000]), np.sum(n4[13:1000]), np.sum(n5[13:1000])])
# frac23 = np.array([np.sum(n1[23:1000]), np.sum(n2[23:1000]), np.sum(n3[23:1000]), np.sum(n4[23:1000]), np.sum(n5[23:1000])])
# frac8 /= np.sum(frac8)
# frac9 /= np.sum(frac9)
# frac13 /= np.sum(frac13)
# frac23 /= np.sum(frac23)

# # Plot of PDF of mass bins
# plt.figure(figsize=(10,8))
# # lower limit of x-axis
# ml = m_lims[Dict[mtyp], 0]
# # upper limit of x-axis
# mu = m_lims[Dict[mtyp], -1]
# # PDFs for cases: SNR > 8, 9, 13, 23
# plt.plot(np.arange(ml, mu, (mu - ml)/5.), frac8, 'bo-',  label='SNR>8')
# plt.plot(np.arange(ml, mu, (mu - ml)/5.), frac9, 'co-',  label='SNR>9')
# plt.plot(np.arange(ml, mu, (mu - ml)/5.), frac13, 'go-',  label='SNR>13')
# plt.plot(np.arange(ml, mu, (mu - ml)/5.), frac23, 'mo-', label='SNR>23')

# # Observations made so far. First array of Total Masses of binaries and second that of Chirp Masses.
# Obs = np.array([[65., 36., 21.7, 50.6, 19., 55.8], [30., 15.1, 8.9, 21.1, 7.92, 24.16]])
# # # Observed data
# plt.axvline(x=Obs[Dict[mtyp], 0], ymin=0., ymax = 1., linewidth=2, color='m', alpha=0.5)  #GW150914
# plt.axvline(x=Obs[Dict[mtyp], 1], ymin=0., ymax = 1., linewidth=2, color='c', alpha=0.5)  #LVT151012
# plt.axvline(x=Obs[Dict[mtyp], 2], ymin=0., ymax = 1., linewidth=2, color='g', alpha=0.5)  #GW151226
# plt.axvline(x=Obs[Dict[mtyp], 3], ymin=0., ymax = 1., linewidth=2, color='g', alpha=0.5)  #GW170104
# plt.axvline(x=Obs[Dict[mtyp], 4], ymin=0., ymax = 1., linewidth=2, color='g', alpha=0.5)  #GW170608
# plt.axvline(x=Obs[Dict[mtyp], 5], ymin=0., ymax = 1., linewidth=2, color='g', alpha=0.5)  #GW170814


# # Annotations
# plt.annotate(' GW150914\n (SNR 24) ', (Obs[Dict[mtyp], 0]+0.4, 0.035), fontsize=12, fontweight='bold', color='m')
# plt.annotate(' LVT151012\n (SNR 9.7) ', (Obs[Dict[mtyp], 1]+0.4, 0.035), fontsize=12, fontweight='bold', color='c')
# plt.annotate(' GW151226\n (SNR 13) ', (Obs[Dict[mtyp], 2]+0.4, 0.035), fontsize=12, fontweight='bold', color='g')
# plt.annotate(' GW170104\n (SNR 13) ', (Obs[Dict[mtyp], 3]+0.4, 0.035), fontsize=12, fontweight='bold', color='g')
# plt.annotate(' GW170608\n (SNR 13) ', (Obs[Dict[mtyp], 4]+0.4, 0.035), fontsize=12, fontweight='bold', color='g')
# plt.annotate(' GW170814\n (SNR 18) ', (Obs[Dict[mtyp], 5]+0.4, 0.035), fontsize=12, fontweight='bold', color='g')

# plt.ylabel('Probability of Detection (SNR>9, SNR>13 amd SNR>23)')

# plt.xlabel(mtyp)
# plt.xticks(np.linspace(0., m_lims[Dict[mtyp], -1], 11))
# plt.title('Probability of Detection Vs. %s (Distribution: %s)' %(mtyp, inp_distri))

# plt.grid(True)
# #plt.yticks(np.linspace(0., 0.5, 11))
# #plt.axis([0., 100., 0., 0.7])
# plt.legend(loc='upper right')

# # Plot of histograms against SNR values for different mass bins
# plt.figure(figsize=(10,8))
# Plot(n1, bin1, m_lims[Dict[mtyp], 0], m_lims[Dict[mtyp], 1], 'ko-')
# Plot(n2, bin2, m_lims[Dict[mtyp], 1], m_lims[Dict[mtyp], 2], 'g^-')
# Plot(n3, bin3, m_lims[Dict[mtyp], 2], m_lims[Dict[mtyp], 3], 'rx-')
# Plot(n4, bin4, m_lims[Dict[mtyp], 3], m_lims[Dict[mtyp], 4], 'y*-')
# Plot(n5, bin5, m_lims[Dict[mtyp], 4], m_lims[Dict[mtyp], 5], 'b+-')

# plt.xlabel('SNR')
# plt.ylabel('Probability')
# plt.title('Probability Distribution of SNRs for different bins of %s (Distribution: %s)' %(mtyp, inp_distri))

# plt.grid(True)
# plt.xticks(np.linspace(0., 60., 31))
# # plt.yticks(np.linspace(0, 0.3, 7))
# # plt.axis([0, 60, 0, 0.3])
# plt.legend(loc='upper right')

# plt.show()


###############################################



# Mass of Sun
Msun = 1.989e30
# Get the type of distribution as an input (Exact name to be input in str format)
inp_distri = 'log_flat'
# inp_distri = 'pow_law'

# Dictionary to identify the input mtyp
Dict = {'Total Mass': 2, 'Chirp Mass': 3}
# Reading resp data file into an array
dat = np.genfromtxt('./../Data/Data-for-%s-distri_M1_M2_M_chirpM_rs6_ro1_ro2_rdesign_alpha_delta_iota_SNRs6_SNRo1_SNRo2_SNRdesign.txt' %inp_distri, delimiter=',')

# Plot against Total Mass / Chirp Mass (Exact name to be input in str format)
# mtyp = 'Total Mass'
mtyp = 'Chirp Mass'

# no of bins
nbins = 18
# column no corresponding to run. (S6: 11 .. Design: 14)
run = 13
# SNR > 8.
i8 = np.where(dat[:,run] > 8.)[0]
dat8 = dat[i8, :]
# SNR > 13.
i13 = np.where(dat[:,run] > 13.)[0]
dat13 = dat[i13, :]
# SNR > 23.
i23 = np.where(dat[:,run] > 23.)[0]
dat23 = dat[i23, :]
# hist of chirp/tot mass
plt.figure(figsize=(19, 10))
plt.hist(dat8[:, Dict[mtyp]] / Msun, nbins, color='b', histtype='step', normed=True, label='SNR > 8')
plt.hist(dat13[:, Dict[mtyp]] / Msun, nbins, color='g', histtype='step', normed=True, label='SNR > 13')
plt.hist(dat23[:, Dict[mtyp]] / Msun, nbins, color='m', histtype='step', normed=True, label='SNR > 23')

# # Observed data
Obs = {"GW150914":{"Total Mass": 65., "Chirp Mass": 30., "color": 'm', "SNR": 24}, "LVT151012":{"Total Mass": 36., "Chirp Mass": 15.1, "color": 'c', "SNR": 9.7}, "GW151226":{"Total Mass": 21.7, "Chirp Mass": 8.9, "color": 'g', "SNR": 13}, "GW170104":{"Total Mass": 50.6, "Chirp Mass": 21.1, "color": 'g', "SNR": 13}, "GW170608":{"Total Mass": 19., "Chirp Mass": 7.92, "color": 'g', "SNR": 13}, "GW170814":{"Total Mass": 55.8, "Chirp Mass": 24.16, "color": 'g', "SNR": 18}}

for event in Obs.keys():
    plt.axvline(x=Obs[event[mtyp]], ymin=0., ymax = 1., linewidth=2, color=Obs[event["color"]], alpha=0.5)

# Annotations
annot_y = 0.0045     # Has to be adjusted manually
for event in Obs.keys():
    plt.annotate(' {} (SNR {}) '.format(event, Obs[event["SNR"]]), (Obs[event[mtyp]]+0.4, annot_y), fontsize=12, fontweight='bold', color=Obs[event["color"]], rotation=90)

plt.xlabel(mtyp)
# plt.xlabel('Chirp Mass')
plt.ylabel('Probability of Detection (SNR>8, SNR>13 amd SNR>23)')
plt.title('Probability of Detection Vs. %s (Distribution: %s)' %(mtyp, inp_distri))
plt.legend(loc='upper right')
# plt.axes([0., 100., 0., 1.])

if mtyp == "Chirp Mass":
    plt.xlim([0., 50.])
    plt.ylim([0., 0.06])
else:
    plt.xlim([0., 100.])
    plt.ylim([0., 0.02])

plt.savefig('./../Plots/O2-%s-%s' %(inp_distri, mtyp))
plt.show()
