from function_SNR import *

class Distri:
    def __init__(self, Iter, r, filenm, Md):
        self.Iter = Iter
        # We create an empty data array of the required size
        self.data = np.zeros(Iter*Md*9).reshape(Iter*Md, 9)
        self.r = r
        self.filenm = filenm
        self.rf = open(('Data-for-' + filenm + '-distri_M1_M2_M_chirpM_r_alpha_delta_iota_SNR.csv'), 'w')
        self.Md = Md
    def a_d_i(self):
        alpha = np.random.uniform(0., 2.*np.pi, self.Md)
        delta = np.random.uniform(-np.pi/2., np.pi/2., self.Md)
        iota = np.random.uniform(0., np.pi, self.Md)
        return alpha, delta, iota
    def M_chirpM(self, m1, m2):
        # Total mass array
        M = m2 + m1
        # Reduced mass array
        mu = (m1 * m2) / M
        n = mu / M
        # Chirp mass array
        chirpM = M * (n ** (3./5.))
        return M, chirpM
    def Uniform(self):
        # We make 'Iter' no of iterations to fill the array with data
        for count in range(self.Iter):
            # Getting time of one iteration, to calculate the estimate of remaining time
            if count==0: iter_time = time.time()
            m1 = np.random.uniform(5. * Msun, 50. * Msun, self.Md)
            m2 = np.random.uniform(5. * Msun, 50. * Msun, self.Md)
            alpha, delta, iota = self.a_d_i()
            SNR = find_SNR(m1, m2, self.r, alpha, delta, iota)
            M, chirpM = self.M_chirpM(m1, m2)
            self.data[count*self.Md:(count+1)*self.Md] = np.array([m1, m2, M, chirpM, self.r, alpha, delta, iota, SNR]).transpose()
            if count==0: iter_time2 = time.time()
            print "Remaining time: ", (self.Iter - count) * (iter_time2 - iter_time) / 60., " Min"
        # We sort the data such that the rows are sorted with ascending value of chirp masses of binaries
        self.data = self.data[self.data[:,3].argsort()]
    def log_flat(self):
        # We make 'Iter' no of iterations to fill the array with data
        for count in range(self.Iter):
            # Getting time of one iteration, to calculate the estimate of remaining time
            if count==0: iter_time = time.time()
            m1 = np.random.uniform(np.log(5. * Msun), np.log(50. * Msun), self.Md)
            m1 = np.exp(m1)
            m2 = np.random.uniform(np.log(5. * Msun), np.log((100. * Msun - m1)), self.Md)
            m2 = np.exp(m2)
            alpha, delta, iota = self.a_d_i()
            SNR = find_SNR(m1, m2, self.r, alpha, delta, iota)
            M, chirpM = self.M_chirpM(m1, m2)
            self.data[count*self.Md:(count+1)*self.Md] = np.array([m1, m2, M, chirpM, self.r, alpha, delta, iota, SNR]).transpose()
            if count==0: iter_time2 = time.time()
            print "Remaining time: ", (self.Iter - count) * (iter_time2 - iter_time) / 60., " Min"
        # We sort the data such that the rows are sorted with ascending value of chirp masses of binaries
        self.data = self.data[self.data[:,3].argsort()]
    def pow_law(self):
        kk = 20
        for count in range(Iter/kk):
            # Getting time of one iteration, to calculate the estimate of remaining time
            if count==0: iter_time = time.time()
            ALPHA = 2.3
            A = 5.**(1.-ALPHA) - 95.**(1.-ALPHA)
            C = 2. + 95.**(1.-ALPHA) / A
            m1 = np.random.uniform(1., 2., Md)
            m1 = (A * (C - m1)) ** (1./(1.-ALPHA))
            m1 *= Msun
            m_uppr = np.zeros(Md)     # We create an array M_uppr to choose an upper limit on choice of M2
            for i in range(Md):       # This business is done to ensure that the conditions M1>M2 & M1+M2<=100 Msun are obeyed
                if m1[i] < 50.*Msun: m_uppr[i] = m1[i]
                else: m_uppr[i] = 100.*Msun - m1[i]
            for j in range(kk):
                m2 = np.random.uniform(5.*Msun, m_uppr, Md)
                print m1, m2
                alpha, delta, iota = self.a_d_i()
                SNR = find_SNR(m1, m2, self.r, alpha, delta, iota)
                M, chirpM = self.M_chirpM(m1, m2)
                # print np.array([m1, m2, M, chirpM, self.r, alpha, delta, iota, SNR]).transpose()
                self.data[(count*kk+j)*self.Md:(count*kk+(j+1))*self.Md] = np.array([m1, m2, M, chirpM, self.r, alpha, delta, iota, SNR]).transpose()
                if count==0: iter_time2 = time.time()
                print "Remaining time: ", (self.Iter - count) * (iter_time2 - iter_time) / 60., " Min"
        # We sort the data such that the rows are sorted with ascending value of chirp masses of binaries
        self.data = self.data[self.data[:,3].argsort()]
        print self.data
    def write(self):
        # Top row in the data file written
        title = 'BH1 Mass (kg), BH2 Mass (kg), Total Mass (kg), Chirp Mass(kg), Distance (m), alpha, delta, iota, SNR'
        np.savetxt(self.rf, self.data, delimiter=',', header=title, newline='\n')
        self.rf.close()

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

inp_distri = input('Which distribution to generate data for? (Uniform, log_flat, pow_law)\n: ')
run = Distri(Iter, r, inp_distri, Md)
eval('run.' + str(inp_distri) + '()')
run.write()

print("End Time: %s seconds" % (time.time() - start_time))
