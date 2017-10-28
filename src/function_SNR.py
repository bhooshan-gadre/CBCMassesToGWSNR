import numpy as np
from scipy import interpolate
import time

# Which run?
RUN = input('Which run to generate data for? (S6, O1, O2, Design)\n: ')

# Defining a dictionary for detector ASDs and LIGO runs
# These names will be used to pick ASD files
detector = {0:'S6_L1', 1:'S6_H1', 2:'S6_V1', 3:'O1_L1', 4:'O1_H1', 5:'O2_L1', 6:'O2_H1', 7:'O2_V1', 8:'Design_L', 9:'Design_H', 10:'Design_V'}
# Dictionary to tell which files to pick for a particular run
for_run = {'S6':[0, 1, 2], 'O1':[3, 4], 'O2':[5, 6, 7], 'Design':[8, 9, 10]}

# To check the time taken for 2D interpolations
start_time = time.time()
# Reading the amplitude spectral densities of detectors into an array 'ASDdata'
ASDdata = []
for i in for_run[RUN]:
	ASDdata.append(np.genfromtxt('/home/shreejit/asd_' + detector[i] +'.txt'))

ASDdata = np.array(ASDdata)

# Reading the sky position dependent detector sensitivity data
# Format of data: alpha, delta, F_plus, F_cross
sky_lookup = np.genfromtxt('/home/shreejit/H1_allSky_antenna_patterns.dat')
# 2D Interpolation of data to obtain F+ and Fx functions of right ascension and declination
find_F_plus = interpolate.interp2d(sky_lookup[:,0], sky_lookup[:,1], sky_lookup[:,2])
find_F_cross = interpolate.interp2d(sky_lookup[:,0], sky_lookup[:,1], sky_lookup[:,3])
# clean memory
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
freq = np.arange(9., 8000., df)

# finding integrands for different detectors in the RUN
integrand = np.zeros(len(detector)*len(freq)).reshape(len(detector), len(freq))
for j in for_run[RUN]:
	# Interpolating the asd data (in log-log) to obtain an interpolated continuous function of frequency
	asd = interpolate.interp1d(np.log(ASDdata[j][:, 0]), np.log(ASDdata[j][:, 1]))
	# Calculating ASD values at respective elements of frequency array
	asd_at_freq = asd(np.log(freq))
	asd_at_freq = np.exp(asd_at_freq)
	# Now the frequency dependent integrand of the SNR calculation is calculated for the entire frequency range.
	# Only the values till ISCO frequency will be summed up from this array to get the integration.
	integrand[j] = freq ** (-7. / 3.)
	integrand[j] /= asd_at_freq ** 2.
	integrand[j] *= df
#clean memory
del ASDdata

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

	# Amplitude calculation
	h_tilde = np.sqrt(5. / np.pi / 24.) * (Const * c)
	h_tilde /= dist	
	h_tilde *= (Const * np.pi) ** (-1./6.)
	h_tilde **= 2.

	# Calculation of ISCO frequency (array)
	f_isco = c ** 3. / (6. ** 1.5 * np.pi * G * M)
	# Now, finding the index of each ISCO freq in freq array
	# the integrand values are to be summed till this index to get integration value
	isco_index = ((f_isco - freq[0]) / df).astype(int)
	
	# Making an SNR matrix of required dimension and then substituting values
	SNR = np.zeros(len(M))
	for ii in range(len(M)):
		for jj in for_run[RUN]:
			# Adding squares of SNRs in individual detectors
			SNR[ii] += (4. * h_tilde[ii] * np.sum(integrand[jj][:isco_index[ii]]))
	# And taking root of this quantity
	SNR = np.sqrt(SNR)
	return SNR
	
	
# Defining a more precise function which takes array of masses of BHs (in kgs), array of distance (in m), 
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
	F_plus = np.zeros(len(alpha))
	F_cross = np.zeros(len(alpha))
	for i in range(len(alpha)):
		F_plus[i] = find_F_plus(alpha[i], delta[i])
		F_cross[i] = find_F_cross(alpha[i], delta[i])
	
	# Calculating the effective distance to the source
	D_eff = (F_plus * (1. + (np.cos(iota)) ** 2.) / 2.) ** 2.
	D_eff += (F_cross * np.cos(iota)) ** 2.
	D_eff **= -0.5
	D_eff = D_eff * dist

	#Calculation of amplitude
	h_tilde = np.sqrt(5. / np.pi / 24.) * (Const * c)
	h_tilde /= D_eff
	h_tilde *= (Const * np.pi) ** (-1./6.)
	h_tilde **= 2.

	# Calculation of ISCO frequency (array)
	f_isco = c ** 3. / (6. ** 1.5 * np.pi * G * M)
	# Now, finding the index of each ISCO freq in freq array
	# the integrand values are to be summed till this index to get integration value
	isco_index = ((f_isco - freq[0]) / df).astype(int)

	# Making an SNR matrix of required dimension and then substituting values
	SNR = np.zeros(len(M))
	for ii in range(len(M)):
		for jj in for_run[RUN]:
			# Adding squares of SNRs in individual detectors
			SNR[ii] += (4. * h_tilde[ii] * np.sum(integrand[jj][:isco_index[ii]]))
	# And taking root of this quantity
	SNR = np.sqrt(SNR)
	return SNR
