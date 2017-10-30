from function_SNR import *
# We define a class which contains methods for uniform, log flat and power law distributions
class Distri:
    def __init__(self, Iter, r, filenm, Md):
        self.Iter = Iter        # No of iterations
        self.Md = Md            # No of injections per iteration (or vol discretization)
        # We create an empty data array of the required size
        self.data = np.zeros(Iter*Md*9).reshape(Iter*Md, 9)     # array in which all the data will be stored
        self.r = r              # distance array such that vol is discretized equally
        self.filenm = filenm    # distribution name used to read file
        self.rf = open(('Data-for-%s-distri_M1_M2_M_chirpM_r_alpha_delta_iota_SNR.csv' %filenm), 'w')   # Read corresponding file
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
        # We make 'Iter' no of iterations to fill the data matrix with values
        for count in range(self.Iter):
            # Getting time of 1st iteration
            if count==0: iter_time = time.time
            # Uniform component masses
            m1 = np.random.uniform(5. * Msun, 50. * Msun, self.Md)
            m2 = np.random.uniform(5. * Msun, 50. * Msun, self.Md)
            # Sky position and orientation
            alpha, delta, iota = self.a_d_i()
            # SNRs
            SNR = find_SNR(m1, m2, self.r, alpha, delta, iota)
            M, chirpM = self.M_chirpM(m1, m2)
            # Inputting data into resp 'pockets' of data matrix
            self.data[count*self.Md:(count+1)*self.Md] = np.array([m1, m2, M, chirpM, self.r, alpha, delta, iota, SNR]).transpose()
            # time for 1st iteration
            if count==0: iter_time2 = time.time()
            print "Remaining time: %s Min" %((self.Iter - count) * (iter_time2 - iter_time) / 60.)
        # We sort the data such that the rows are sorted with ascending value of chirp masses of binaries
        self.data = self.data[self.data[:,3].argsort()]
    # Calculating data for log flat distribution
    def log_flat(self):
        # We make 'Iter' no of iterations to fill the array with data
        for count in range(self.Iter):
            if count==0: iter_time = time.time()
            # generating component masses uniform in log
            m1 = np.random.uniform(np.log(5. * Msun), np.log(50. * Msun), self.Md)
            m1 = np.exp(m1)
            m2 = np.random.uniform(np.log(5. * Msun), np.log((100. * Msun - m1)), self.Md)
            m2 = np.exp(m2)
            alpha, delta, iota = self.a_d_i()
            SNR = find_SNR(m1, m2, self.r, alpha, delta, iota)
            M, chirpM = self.M_chirpM(m1, m2)
            self.data[count*self.Md:(count+1)*self.Md] = np.array([m1, m2, M, chirpM, self.r, alpha, delta, iota, SNR]).transpose()
            if count==0: iter_time2 = time.time()
            print "Remaining time: %s Min" %((self.Iter - count) * (iter_time2 - iter_time) / 60.)
        self.data = self.data[self.data[:,3].argsort()]
    # Calculating data for power law distribution
    def pow_law(self):
        # To ensure uniform distribution of m2 upto m1 (20 values of m2 for each m1)
        # We sub-iterate m2 kk times for each m1
        kk = 20
        for count in range(Iter/kk):
            if count==0: iter_time = time.time()
            # Power to which 1st component mass is distributed, P(m1) = m1^-ALPHA
            ALPHA = 2.3
            A = 5.**(1.-ALPHA) - 95.**(1.-ALPHA)
            C = 2. + 95.**(1.-ALPHA) / A
            m1 = np.random.uniform(1., 2., Md)
            m1 = (A * (C - m1)) ** (1./(1.-ALPHA))
            m1 *= Msun
            m_uppr = np.zeros(Md)     # We create an array m_uppr to choose an upper limit on choice of m2
            for i in range(Md):       # This business is done to ensure that the conditions m1>m2 & m1+m2<=100 Msun are obeyed
                if m1[i] < 50.*Msun: m_uppr[i] = m1[i]
                else: m_uppr[i] = 100.*Msun - m1[i]
            # We sub-iterate kk times for uniform distribution of m2 upto m_uppr, for each m1
            for j in range(kk):
                m2 = np.random.uniform(5.*Msun, m_uppr, Md)
                alpha, delta, iota = self.a_d_i()
                SNR = find_SNR(m1, m2, self.r, alpha, delta, iota)
                M, chirpM = self.M_chirpM(m1, m2)
                self.data[(count*kk+j)*self.Md:(count*kk+(j+1))*self.Md] = np.array([m1, m2, M, chirpM, self.r, alpha, delta, iota, SNR]).transpose()
                if count==0: iter_time2 = time.time()
                print "Remaining time: %s Min" %((self.Iter - count) * (iter_time2 - iter_time) / 60.)
        # We sort the data such that the rows are sorted with descending value of chirp masses of binaries
        self.data = self.data[self.data[:,3].argsort()]
        print self.data
    # Write the data into resp file
    def write(self):
        # Top row in the data file written
        title = 'BH1 Mass (kg), BH2 Mass (kg), Total Mass (kg), Chirp Mass(kg), Distance (m), alpha, delta, iota, SNR'
        np.savetxt(self.rf, self.data, delimiter=',', header=title, newline='\n')
        self.rf.close()

# To check the total time of run
start_time = time.time()
print("Start Time: %s seconds" % (time.time() - start_time))

# Constants
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

# Get the type of distribution as an input (Exact name to be input in str format)
inp_distri = input('Which distribution to generate data for? (Uniform, log_flat, pow_law)\n: ')
run = Distri(Iter, r, inp_distri, Md)       # An object of the class
eval('run.%s()' %inp_distri)                # To run the corresponding method for distribution
run.write()                                 # Write the file

print("End Time: %s seconds" % (time.time() - start_time))
