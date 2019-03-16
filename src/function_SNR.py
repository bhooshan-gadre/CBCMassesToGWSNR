import numpy as np
from scipy import interpolate
import time
import h5py

from pycbc.waveform import get_fd_waveform
import pycbc.filter.matchedfilter as mf
import pycbc.psd as pp
import pycbc.types as pt
from pycbc.population.scale_injections import dlum_to_z

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

# To check the time taken for 2D interpolations
start_time = time.time()
# Reading the amplitude spectral densities of all detectors into an array 'ASDdata'

ASDdata = {}
# Format : ASDdata = {'S6_H1':<Nx2 ASD data array>, 'O1_L1':<Nx2 ASD data array>}
for i in detector.keys(): 		#for_run[RUN]:
	ASDdata[detector[i]] = np.genfromtxt('./../Data/asd_%s.txt' %detector[i], delimiter=",")
# Dict of PSDs
PSD_for_det = {}
for det in detector.keys():
	PSD_for_det[detector[det]] = pp.read.from_numpy_arrays(ASDdata[detector[det]][:,0], ASDdata[detector[det]][:,1]**2., int(1500/df), delta_f=df, low_freq_cutoff=freq[0])


############ Below part is needed only for find_SNR and find_simple_SNR ############

# # Reading the sky position dependent detector sensitivity data
# # Format of data: alpha, delta, F_plus, F_cross
# sky_lookup = np.genfromtxt('./../Data/H1_allSky_antenna_patterns.dat')
# # 2D Interpolation of data to obtain F+ and Fx functions of right ascension and declination
# find_F_plus = interpolate.interp2d(sky_lookup[:,0], sky_lookup[:,1], sky_lookup[:,2])
# find_F_cross = interpolate.interp2d(sky_lookup[:,0], sky_lookup[:,1], sky_lookup[:,3])
# #clear memory
# del sky_lookup
# # Printing the time taken for reading and interpolation stuff
# print("Time taken to interpolate: %s seconds" % (time.time() - start_time))

# # Initializing an array of integrands for all detectors
# integrand = np.zeros([len(detector), len(freq)])
# # finding integrands for different detectors for all RUNs
# for j in detector.keys():		#for_run[RUN]:
# 	# Interpolating the asd data to obtain an interpolated continuous function of frequency
# 	asd = interpolate.interp1d(ASDdata[detector[j]][:, 0], ASDdata[detector[j]][:, 1])
# 	# Calculating ASD values at respective elements of frequency array
# 	asd_at_freq = asd(freq)
# 	# Now the frequency dependent integrand of the SNR calculation is calculated for the entire frequency range.
# 	# Only the values till ISCO frequency will be summed up from this array to get the integration.
# 	integrand[j] = freq ** (-7. / 3.)
# 	integrand[j] /= asd_at_freq ** 2.
# 	integrand[j] *= df

# # clean
# del ASDdata
# # Cumulative sum of rows of integrand matrix along the rows
# # So that the integration upto a certain ISCO freq will only pick up the value corresponding to the resp index, instead of summing f_low to f_isco everytime
# integrand = np.cumsum(integrand, axis=1)

############



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
	SNRs_for_RUN = {}
	M1 = M1 / Msun
	M2 = M2 / Msun
	# print "m1: {}, m2: {}, r: {}".format(M1, M2, dist)
	# Making an SNR dictionary for SNR arrays corresponding to diff RUNs
	SNRs_for_RUN = {}
	# For each RUN like "O1"
	for RUN in for_run.keys():
		Dist = dist[RUN] / Mpc
		# We get an array of SNRs of the size of M1
		SNRs_for_RUN[RUN] = np.zeros(len(M1))
		# Calculation for each binary
		for i in range(len(M1)):
			if i == 2:
				start_t = time.time()
			# Strain calculation only for (+) polarization
			sptilde, _ = get_fd_waveform(approximant="IMRPhenomD", mass1=M1[i], mass2=M2[i], distance=Dist[i], inclination=iota[i], delta_f=df, f_lower=freq[0], f_final=freq[-1])
			sptilde.resize(len(PSD_for_det['O1_L1']))
			# root of sum of sigma squares in individual detectors
			for det in for_run[RUN]:
				# SNR in each detector
				SNRs_for_RUN[RUN][i] += mf.sigmasq(sptilde, psd=PSD_for_det[detector[det]], low_frequency_cutoff=freq[0])
			if i == 2:
				print "For {} run, expected time: {} min".format(RUN, (time.time() - start_t) * len(M1) / 60.)
		# taking root to get values of SNR
		SNRs_for_RUN[RUN] = np.sqrt(SNRs_for_RUN[RUN])
	return SNRs_for_RUN

# Find pycbc based SNR from saved templates
def find_SNR_frm_hdf(path, group, iters):
	Temps = h5py.File(path, "r")
	Alpha, Delta, Iota = (Temps["alpha"].value, Temps["delta"].value, Temps["iota"].value)
	M1 = Temps[group+"/Mass1"].value
	M2 = Temps[group+"/Mass2"].value
	M = M2 + M1
	# Chirp mass array of binary
	chirpM = (M1 * M2)
	chirpM **= (3./5.)
	chirpM /= M ** (1./5.)
	# table of templates
	FD_temps = Temps[group+"/templates"]
	# Dict of distances, and redshifts
	dist = {}
	z = {}
	# Making an SNR dictionary for SNR arrays corresponding to diff RUNs
	SNRs_for_RUN = {}
	# For each RUN like "O1"
	for RUN in for_run:
		print "RUN: {}".format(RUN)
		# We get an array of SNRs of the size of M1
		# Initializing
		SNRs_for_RUN[RUN] = np.zeros(len(FD_temps)*iters)
		# Loading distance array for the specific run
		dist[RUN] = Temps["Distance/{}".format(RUN)]
		# Corresponding redshifts
		z[RUN] = dlum_to_z(dist[RUN].value/Mpc)
		# For 'iters' number of iterations, the template set will be used with shuffled distances
		for itr in np.arange(iters):
			# Calculation for each binary
			for j in range(len(FD_temps)):
				if j % 100 == 0:
					print "{} / {}".format(j+len(FD_temps)*itr, len(FD_temps)*iters)
				# generating frequency series template
				sptilde = pt.frequencyseries.FrequencySeries(FD_temps[j], delta_f=df) * 100. * Mpc / dist[RUN][j+len(FD_temps)*itr]
				# sample freq values
				temp_sample_freqs = sptilde.sample_frequencies.data
				# redshifted freqs will map as f -> f/(1+z)
				f_red = temp_sample_freqs/(1. + z[RUN][j+len(FD_temps)*itr])
				# Appending the redshifted frequencies till the end of prev frequency array i.e. freq[-1]
				f_red_ext = np.append(f_red, np.linspace(f_red[-1], temp_sample_freqs[-1], 1000)[1:])
				# appending the sptilde values with equal no of zeros
				sptilde_ext = np.append(sptilde.data, np.zeros(999))
				# interpolation function
				sp_at_redshifted_f = interpolate.interp1d(f_red_ext, sptilde_ext)
				# interpolated sptilde values 
				sptilde.data = sp_at_redshifted_f(temp_sample_freqs)
				# root of sum of sigma squares in individual detectors
				for det in for_run[RUN]:
					# SNR in each detector
					SNRs_for_RUN[RUN][j+len(FD_temps)*itr] += mf.sigmasq(sptilde, psd=PSD_for_det[detector[det]], low_frequency_cutoff=freq[0])
		# taking root to get values of SNR
		SNRs_for_RUN[RUN] = np.sqrt(SNRs_for_RUN[RUN])
	return np.tile(M1, iters), np.tile(M2, iters), np.tile(M, iters), np.tile(chirpM, iters), dist, np.tile(Alpha, iters), np.tile(Delta, iters), np.tile(Iota, iters), SNRs_for_RUN
