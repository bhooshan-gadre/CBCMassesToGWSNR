from function_SNR import *
# We define a class which contains methods for uniform, log flat and power law distributions
class Distri:
    def __init__(self, r, filenm, Md):
        self.Md = Md            # No of injections per iteration (or vol discretization)
        # We create an empty data array of the required size
        self.data = np.zeros([Md, 15])     # array in which all the data will be stored
        self.r = r              # distance array such that vol is discretized equally
        self.rf = open(('./../Data/Data-for-%s-distri_M1_M2_M_chirpM_rs6_ro1_ro2_rdesign_alpha_delta_iota_SNRs6_SNRo1_SNRo2_SNRdesign.txt' %filenm), 'w')   # Read corresponding file
    # Generating random sky position and orientation (arrays)
    def a_d_i(self):
        alpha = np.random.uniform(0., 2.*np.pi, self.Md)
        delta = np.random.uniform(-np.pi/2., np.pi/2., self.Md)
        iota = np.random.uniform(0., np.pi, self.Md)
        return alpha, delta, iota
    # Calculating total mass and chirp mass arrays
    def M_chirpM(self, m1, m2):
        # Total mass array
        M = m2 + m1
        # Reduced mass array
        mu = (m1 * m2) / M
        n = mu / M
        # Chirp mass array
        chirpM = M * (n ** (3./5.))
        return M, chirpM
    # Calculating data for uniform distribution
    def Uniform(self):
        # Uniform component masses
        m1 = np.random.uniform(5. * Msun, 50. * Msun, self.Md)
        m2 = np.random.uniform(5. * Msun, 50. * Msun, self.Md)
        # Sky position and orientation
        alpha, delta, iota = self.a_d_i()
        # SNRs
        SNR_for_RUN = find_SNR(m1, m2, self.r, alpha, delta, iota)
        M, chirpM = self.M_chirpM(m1, m2)
        # Inputting data into resp 'pockets' of data matrix
        self.data = np.array([m1, m2, M, chirpM, self.r['S6'], self.r['O1'], self.r['O2'], self.r['Design'], alpha, delta, iota, SNR_for_RUN['S6'], SNR_for_RUN['O1'], SNR_for_RUN['O2'], SNR_for_RUN['Design']]).transpose()
        # We sort the data such that the rows are sorted with ascending value of chirp masses of binaries
        self.data = self.data[self.data[:,3].argsort()]
    # Calculating data for log flat distribution
    def log_flat(self):
        # generating component masses uniform in log
        # Deliberately sampling 10 times higher number of samples
        m1 = np.random.uniform(np.log(5. * Msun), np.log(95. * Msun), 10*self.Md)
        m1 = np.exp(m1)
        m2 = np.random.uniform(np.log(5. * Msun), np.log(95. * Msun), 10*self.Md)
        m2 = np.exp(m2)
        # Indices of samples whose total mass is lower than 100 Msun
        i_accpt = np.where((m1+m2) < 100.*Msun)
        # Take only "Md" combinations out of the acceptable ones
        m1 = m1[i_accpt][:self.Md]
        m2 = m2[i_accpt][:self.Md]
        # Sky position and orientation
        alpha, delta, iota = self.a_d_i()
        # SNRs
        SNR_for_RUN = find_SNR(m1, m2, self.r, alpha, delta, iota)
        M, chirpM = self.M_chirpM(m1, m2)
        self.data = np.array([m1, m2, M, chirpM, self.r['S6'], self.r['O1'], self.r['O2'], self.r['Design'], alpha, delta, iota, SNR_for_RUN['S6'], SNR_for_RUN['O1'], SNR_for_RUN['O2'], SNR_for_RUN['Design']]).transpose()
        self.data = self.data[self.data[:,3].argsort()]
    # Calculating data for power law distribution
    def pow_law(self):
        # Power to which 1st component mass is distributed, P(m1) = m1^-ALPHA
        ALPHA = 2.35
        A = 5.**(1.-ALPHA) - 95.**(1.-ALPHA)
        C = 2. + 95.**(1.-ALPHA) / A
        # Deliberately sampling 10 times higher number of samples
        m1 = np.random.uniform(1., 2., 10*self.Md)
        m1 = (A * (C - m1)) ** (1./(1.-ALPHA))
        m1 *= Msun
        m2 = np.random.uniform(5.*Msun, m1, 10*self.Md)
        # Indices of samples whose total mass is lower than 100 Msun
        i_accpt = np.where((m1+m2) < 100.*Msun)
        # Take only "Md" combinations out of the acceptable ones
        m1 = m1[i_accpt][:self.Md]
        m2 = m2[i_accpt][:self.Md]
        # Sky position and orientation
        alpha, delta, iota = self.a_d_i()
        # SNRs
        SNR_for_RUN = find_SNR(m1, m2, self.r, alpha, delta, iota)
        M, chirpM = self.M_chirpM(m1, m2)
        self.data = np.array([m1, m2, M, chirpM, self.r['S6'], self.r['O1'], self.r['O2'], self.r['Design'], alpha, delta, iota, SNR_for_RUN['S6'], SNR_for_RUN['O1'], SNR_for_RUN['O2'], SNR_for_RUN['Design']]).transpose()
        # We sort the data such that the rows are sorted with descending value of chirp masses of binaries
        self.data = self.data[self.data[:,3].argsort()]
    # Write the data into resp file
    def write(self):
        # Top row in the data file written
        title = 'BH1 Mass (kg), BH2 Mass (kg), Total Mass (kg), Chirp Mass(kg), Distance_S6 (m), Distance_O1 (m), Distance_O2 (m), Distance_Design (m), alpha, delta, iota, SNR_S6, SNR_O1, SNR_O2, SNR_Design'
        np.savetxt(self.rf, self.data, delimiter=',', header=title, newline='\n')
        self.rf.close()

# To check the total time of run
start_time = time.time()
print("Start Time: %s seconds" % (time.time() - start_time))

np.random.seed(123)
Md = 1000000        # No of parts in which volume will be discretized
# Empty dict of "r"
r = {}
for RUN in r_extreme['Lower Cutoff'].keys():
    # Limits on space corresponding to respective RUNs
    rl = r_extreme['Lower Cutoff'][RUN]    # Lower limit on space
    ru = r_extreme['Upper Cutoff'][RUN]   # Upper limit on space
    # Defining the visible volume
    vol = 4. / 3. * np.pi * (ru ** 3. - rl ** 3.)
    # From uniform discretization of volume of space, we obtain corresponding discretization of radial distance
    r[RUN] = np.linspace(0., vol, Md//100)
    r[RUN] *= 3. / (4. * np.pi)
    r[RUN] += rl ** 3.
    r[RUN] **= (1. / 3.)
    # We generate array of "r" as Md samples from above distance discretization
    r[RUN] = np.random.choice(r, Md)

# Get the type of distribution as an input (Exact name to be input in str format)
inp_distri = input('Which distribution to generate data for? (Uniform, log_flat, pow_law)\n: ')
run = Distri(r, inp_distri, Md)     # An object of the class
eval('run.%s()' %inp_distri)        # To run the corresponding method for distribution
run.write()                         # Write the file

print("End Time: %s seconds" % (time.time() - start_time))
