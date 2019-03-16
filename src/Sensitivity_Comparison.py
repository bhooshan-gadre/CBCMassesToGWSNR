import numpy as np
import pycbc.psd as pp
import matplotlib.pyplot as plt
plt.ion()

# Which run?
# RUN = input('Which run to generate data for? (S6, O1, O2, Design)\n: ')

# Constants
Msun = 1.989e30			# Mass of Sun in 'kg'
G = 6.67259e-11			# Gravitational const
c = 2.99792458e8		# Speed of light in 'm/s'
Mpc = 3.08568025e22		# Megaparsec in 'm'
Md = 1000000			# No of samples in each distribution
Iterations = 10			# No of iterations for which the saved set of templates will be used with shuffled distances
mass_min = 5.			# Mass Limits in Msun
mass_max = 95.

# Array of frequencies
df = 1./4.
# r_extreme = {'Lower Space Cutoff': {'Design': 1587.82130486 * Mpc, 'S6': 82.77513475 * Mpc, 'O1': 454.18511569 * Mpc, 'O2': 576.59669512 * Mpc}, 'Upper Space Cutoff': {'Design': 3279.12309453 * Mpc, 'S6': 153.40926056 * Mpc, 'O1': 901.50600027 * Mpc, 'O2': 1153.03794345 * Mpc}}
r_extreme = {'Lower Space Cutoff': {'Design': 273.38516829 * Mpc, 'S6': 12.80446162 * Mpc, 'O1': 72.82472261 * Mpc, 'O2': 91.44922879 * Mpc}, 'Upper Space Cutoff': {'Design': 19063.31247625 * Mpc, 'S6': 994.86000718 * Mpc, 'O1': 5260.97949502 * Mpc, 'O2': 6597.62537123 * Mpc}}


# Defining a dictionary for detector ASDs and LIGO runs
# These names will be used to pick ASD files
detector = {0:'S6_L1', 1:'S6_H1', 2:'O1_L1', 3:'O1_H1', 4:'O2_L1', 5:'O2_H1', 6:'Design_L', 7:'Design_H', 8:'Design_V'}
# Dictionary to tell which files to pick for a particular run
for_run = {'S6':[0, 1], 'O1':[2, 3], 'O2':[4, 5], 'Design':[6, 7, 8]}

############# Change f_low #############
freq = np.arange(20., 1500., df)
########################################

# Reading the amplitude spectral densities of all detectors into an array 'ASDdata'

ASDdata = {}
# Format : ASDdata = {'S6_H1':<Nx2 ASD data array>, 'O1_L1':<Nx2 ASD data array>}
for i in detector.keys():
	ASDdata[detector[i]] = np.genfromtxt('./../Data/asd_%s.txt' %detector[i], delimiter=",")

# Dict of PSDs
PSD_for_det = {}
for det in detector.keys():
	PSD_for_det[detector[det]] = pp.read.from_numpy_arrays(ASDdata[detector[det]][:,0], ASDdata[detector[det]][:,1]**2., int(1500/df), delta_f=df, low_freq_cutoff=freq[0])


plt.figure(figsize=[12,9])
for RUN, clr in zip(for_run, ['g', 'b', 'm', 'k']):
	alpha = 1.
	for det in for_run[RUN]:
		plt.loglog(PSD_for_det[detector[det]].sample_frequencies.data, PSD_for_det[detector[det]].data**0.5, alpha=alpha, color=clr, label=detector[det])
		alpha -= 0.4

plt.legend()
plt.xlim([20, 1600])
plt.xlabel('Frequency (Hz)')
plt.ylabel('$\mathrm{S}_\mathrm{n}^{1/2}$ $(\mathrm{Hz}^{-1/2})$')
plt.title('Sensitivity Comparison (S6, O1, O2 and Design)')
plt.savefig('Sensitivity_Comparison.png')
plt.show()
