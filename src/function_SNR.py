import numpy as np
from scipy import interpolate
import time

from pycbc.waveform import get_fd_waveform
import pycbc.filter.matchedfilter as mf

# Which run?
# RUN = input('Which run to generate data for? (S6, O1, O2, Design)\n: ')

# Defining a dictionary for detector ASDs and LIGO runs
# These names will be used to pick ASD files
detector = {0:'S6_L1', 1:'S6_H1', 2:'O1_L1', 3:'O1_H1', 4:'O2_L1', 5:'O2_H1', 6:'Design_L', 7:'Design_H', 8:'Design_V'}
# Dictionary to tell which files to pick for a particular run
for_run = {'S6':[0, 1], 'O1':[2, 3], 'O2':[4, 5], 'Design':[6, 7, 8]}

# To check the time taken for 2D interpolations
start_time = time.time()
# Reading the amplitude spectral densities of all detectors into an array 'ASDdata'

ASDdata = {}
# Format : ASDdata = {'S6_H1':<Nx2 ASD data array>, 'O1_L1':<Nx2 ASD data array>}
for i in detector.keys(): 		#for_run[RUN]:
	ASDdata[detector[i]] = np.genfromtxt('./../Data/asd_%s.txt' %detector[i], delimiter=",")

# ASDdata = np.genfromtxt('./../Data/asd_%s.txt' %detector[0], delimiter=",")
# for i in detector.keys()[1:]: 		#for_run[RUN]:
# 	ASDdata = np.append(ASDdata, np.genfromtxt('./../Data/asd_%s.txt' %detector[i], delimiter=","), axis=0)

# # ASDdata = np.array(ASDdata)

# Reading the sky position dependent detector sensitivity data
# Format of data: alpha, delta, F_plus, F_cross
sky_lookup = np.genfromtxt('./../Data/H1_allSky_antenna_patterns.dat')
# 2D Interpolation of data to obtain F+ and Fx functions of right ascension and declination
find_F_plus = interpolate.interp2d(sky_lookup[:,0], sky_lookup[:,1], sky_lookup[:,2])
find_F_cross = interpolate.interp2d(sky_lookup[:,0], sky_lookup[:,1], sky_lookup[:,3])
#clear memory
del sky_lookup
# Printing the time taken for reading and interpolation stuff
print("Time taken to interpolate: %s seconds" % (time.time() - start_time))

# Constants
Msun = 1.989e30			# Mass of Sun in 'kg'
G = 6.67259e-11			# Gravitational const
c = 2.99792458e8		# Speed of light in 'm/s'
Mpc = 3.08568025e22		# Megaparsec in 'm'
# Array of frequencies
df = 0.01
r_extreme = {'Lower Cutoff': {'Design': 1056.*Mpc, 'O1': 179.*Mpc, 'O2': 267.875*Mpc, 'S6': 1.625*Mpc}, 'Upper Cutoff': {'Design': 3452.125*Mpc, 'O1': 964.875*Mpc, 'O2': 1224.375*Mpc, 'S6': 169.*Mpc}}

############# Change f_low #############
freq = np.arange(20., 1500., df)
########################################

# Initializing an array of integrands for all detectors
integrand = np.zeros([len(detector), len(freq)])
# finding integrands for different detectors for all RUNs
for j in detector.keys():		#for_run[RUN]:
	# Interpolating the asd data to obtain an interpolated continuous function of frequency
	asd = interpolate.interp1d(ASDdata[detector[j]][:, 0], ASDdata[detector[j]][:, 1])
	# Calculating ASD values at respective elements of frequency array
	asd_at_freq = asd(freq)
	# Now the frequency dependent integrand of the SNR calculation is calculated for the entire frequency range.
	# Only the values till ISCO frequency will be summed up from this array to get the integration.
	integrand[j] = freq ** (-7. / 3.)
	integrand[j] /= asd_at_freq ** 2.
	integrand[j] *= df
# clean
del ASDdata
# Cumulative sum of rows of integrand matrix along the rows
# So that the integration upto a certain ISCO freq will only pick up the value corresponding to the resp index, instead of summing f_low to f_isco everytime
integrand = np.cumsum(integrand, axis=1)

# Defining a simple function which takes arrays of masses of BHs (in kgs), 
# array of distance (in m) and gives an output array of SNRs
def find_simple_SNR(M1, M2, dist):
	M = M2 + M1
	# Reduced mass array of binary
	mu = (M1 * M2) / M
	n = mu / M
	# Chirp mass array of binary
	chirpM = M * (n ** (3./5.))
	Const = G * chirpM / (c ** 3.)

	# Calculation of ISCO frequency (array)
	f_isco = c ** 3. / (6. ** 1.5 * np.pi * G * M)
	# Now, finding the index of each ISCO freq in freq array
	# the integrand values are to be summed till this index to get integration value
	isco_index = ((f_isco - freq[0]) / df).astype(int)
	
	# Making an SNR dictionary for SNR arrays corresponding to diff RUNs
	SNRs_for_RUN = {}
	for RUN in for_run.keys():
		# Amplitude calculation
		h_tilde_sq = np.sqrt(5. / np.pi / 24.) * (Const * c)
		h_tilde_sq /= dist[RUN]
		h_tilde_sq *= (Const * np.pi) ** (-1./6.)
		h_tilde_sq **= 2.
		# Initialization of SNR array for corresponding RUN
		SNRs_for_RUN[RUN] = np.zeros(len(M))
		for ii in range(len(M)):	# for each CBC injection
			# Adding squares of SNRs in individual detectors
			SNRs_for_RUN[RUN][ii] = (4. * h_tilde_sq[ii] * np.sum(integrand[for_run[RUN], isco_index[ii]]))
		# And taking root to get values of SNR
		SNRs_for_RUN[RUN] = np.sqrt(SNRs_for_RUN[RUN])
	return SNRs_for_RUN


# Defining a more accurate function which takes array of masses of BHs (in kgs), array of distance (in m), 
# the right ascension angle-alpha, declination angle-delta and angle between line of sight and 
# angular momentum vector of binary-iota, to give an array of SNRs
def find_SNR(M1, M2, dist, alpha, delta, iota):
	M = M2 + M1
	# Reduced mass array of binary
	mu = (M1 * M2) / M
	n = mu / M
	# Chirp mass array of binary
	chirpM = M * (n ** (3./5.))
	Const = G * chirpM / (c ** 3.)

	# Creating matrices of required dimensions and substituting the values using interpolation functions for F+ and Fx
	F_plus = np.array([float(find_F_plus(ALPHA, DELTA)) for ALPHA, DELTA in zip(alpha, delta)])
	F_cross = np.array([float(find_F_cross(ALPHA, DELTA)) for ALPHA, DELTA in zip(alpha, delta)])

	# Calculation of ISCO frequency (array)
	f_isco = c ** 3. / (6. ** 1.5 * np.pi * G * M)
	# Now, finding the index of each ISCO freq in freq array
	# the integrand values are to be summed till this index to get integration value
	isco_index = ((f_isco - freq[0]) / df).astype(int)
	
	# Making an SNR dictionary for SNR arrays corresponding to diff RUNs
	SNRs_for_RUN = {}
	for RUN in for_run.keys():
		# Calculating the effective distance to the source
		D_eff = (F_plus * (1. + (np.cos(iota)) ** 2.) / 2.) ** 2.
		D_eff += (F_cross * np.cos(iota)) ** 2.
		D_eff **= -0.5
		D_eff *= dist[RUN]
		#Calculation of amplitude
		h_tilde_sq = np.sqrt(5. / np.pi / 24.) * (Const * c)
		h_tilde_sq /= D_eff
		h_tilde_sq *= (Const * np.pi) ** (-1./6.)
		h_tilde_sq **= 2.
		# Initialization of SNR array for corresponding RUN
		SNRs_for_RUN[RUN] = np.zeros(len(M))
		for ii in range(len(M)):	# for each CBC injection
			# Adding squares of SNRs in individual detectors
			SNRs_for_RUN[RUN][ii] = (4. * h_tilde_sq[ii] * np.sum(integrand[for_run[RUN], isco_index[ii]]))
		# And taking root to get values of SNR
		SNRs_for_RUN[RUN] = np.sqrt(SNRs_for_RUN[RUN])
	return SNRs_for_RUN

# Defining a PyCBC based function which takes array of masses of BHs (in kgs), array of distance (in m), 
# the right ascension angle-alpha, declination angle-delta and angle between line of sight and 
# angular momentum vector of binary-iota, to give an array of SNRs
def find_pycbc_SNR(M1, M2, dist, alpha, delta, iota):
	M1 = M1 / Msun
	M2 = M2 / Msun
	dist = dist / Mpc
	sptilde, sctilde = get_fd_waveform(approximant="TaylorF2", mass1=M1, mass2=M2, distance=dist, inclination=iota, delta_f=df, f_lower=freq[0])
	mf.sigmasq(htilde, psd=None, low_frequency_cutoff=None, high_frequency_cutoff=None)