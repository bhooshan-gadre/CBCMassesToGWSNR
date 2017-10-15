from function_SNR import *
import time

# To check the total time of run
start_time = time.time()
print("Start Time: %s seconds" % (time.time() - start_time))

# Inputs
rl = 648.8 * Mpc    # Lower limit on space
ru = 2030.6 * Mpc   # Upper limit on space
Md = 1000           # No of parts in which volume will be discretized
Iter = 1000         # No of times the code will be run

# Defining the visible volume
vol = 4. / 3. * np.pi * (ru ** 3. - rl ** 3.)

# From uniform discretization of volume of space, we obtain corresponding discretization of radial distance
r = np.linspace(0., vol, Md)
r *= 3. / (4. * np.pi)
r += rl ** 3.
r **= (1. / 3.)

# Open a file for writing the generated data using relevant distribution
rf = open('Data-for-flat-distri_M1_M2_M_chirpM_r_alpha_delta_iota_SNR.csv', 'w')            # Uniform in both masses
# rf = open('Data-for-log-flat-distri_M1_M2_M_chirpM_r_alpha_delta_iota_SNR.csv', 'w')      # Uniform in log of each mass
# rf = open('Data-for-power-law-distri_M1_M2_M_chirpM_r_alpha_delta_iota_SNR.csv', 'w')     # Power Law in primary mass; Uniform in secondary

# We create an empty data array of the required size
data = np.zeros(Iter*Md*9).reshape(Iter*Md, 9)

# We make 'Iter' no of iterations to fill the array with data
# kk = 20
for count in range(Iter): 
# for count in range(Iter/kk):
    
    # Getting time of one iteration, to calculate the estimate of remaining time
    if count==0:
        iter_time = time.time()
    
    # Uncomment relevant patch below for specific model

    # 1. Uniform distribution
    M1 = np.random.uniform(5. * Msun, 50. * Msun, Md)
    M2 = np.random.uniform(5. * Msun, 50. * Msun, Md)

    # 2. Uniform in log
    # M1 = np.random.uniform(np.log(5. * Msun), np.log(50. * Msun), Md)
    
    # M1 = np.exp(M1)
    # M2 = np.random.uniform(np.log(5. * Msun), np.log((100.*Msun - M1)), Md)
    # M2 = np.exp(M2)

    # 3. Power law
    # ALPHA = 2.3
    # A = 5.**(1.-ALPHA) - 95.**(1.-ALPHA)
    # C = 2. + 95.**(1.-ALPHA) / A
    # M1 = np.random.uniform(1., 2., Md)
    # M1 = (A * (C - M1)) ** (1./(1.-ALPHA))
    # M1 *= Msun
    # M_uppr = np.zeros(Md)     # We create an array M_uppr to choose an upper limit
                                # on choice of M2
    
    # for i in range(Md):       # This business is done to ensure that the conditions M1>M2 & M1+M2<=100 Msun are obeyed
    #     if M1[i] < 50.*Msun: M_uppr[i] = M1[i]
    #     else: M_uppr[i] = 100.*Msun - M1[i]

    # for j in range(kk):
    #     M2 = np.random.uniform(5.*Msun, M_uppr, Md)
    
    #### Indent the below section when enabling Power Law Distribution ####

    # Random sky location and orientation wrt line of sight
    alpha = np.random.uniform(0., 2.*np.pi, Md)
    delta = np.random.uniform(-np.pi/2., np.pi/2., Md)
    iota = np.random.uniform(0., np.pi, Md)
    # Total mass array
    M = M2 + M1
    # Reduced mass array
    mu = (M1 * M2) / M
    n = mu / M
    # Chirp mass array
    chirpM = M * (n ** (3./5.))

    # SNR array from the above arrays of parameters
    SNR = find_SNR(M1, M2, r, alpha, delta, iota)

    data[count*Md:(count+1)*Md] = np.array([M1, M2, M, chirpM, r, alpha, delta, iota, SNR]).transpose()
    # data[(count*kk+j)*Md:(count*kk+(j+1))*Md] = np.array([M1, M2, M, chirpM, r, alpha, delta, iota, SNR]).transpose()

    #### Till here ####

    if count==0:
        iter_time2 = time.time()

    print "Remaining time: ", (Iter - count) * (iter_time2 - iter_time) / 60., " Min"
    # print "Remaining time: ", (Iter/kk - count) * (iter_time2 - iter_time) / 60., " Min"

# We sort the data such that the rows are sorted with ascending value of total masses of binaries
# data = data[data[:,2].argsort()]

# We sort the data such that the rows are sorted with ascending value of chirp masses of binaries
data = data[data[:,3].argsort()]

# Top row in the data file written
title = 'BH1 Mass (kg), BH2 Mass (kg), Total Mass (kg), Chirp Mass(kg), Distance (m), alpha, delta, iota, SNR'
np.savetxt(rf, data, delimiter=',', header=title, newline='\n')
rf.close()

print("End Time: %s seconds" % (time.time() - start_time))
