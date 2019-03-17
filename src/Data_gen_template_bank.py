from function_SNR import *
import numpy as np

# Generating random sky position and orientation (arrays)
def a_d_i(Md):
    alpha = np.random.uniform(0., 2.*np.pi, Md)
    delta = np.random.uniform(-np.pi/2., np.pi/2., Md)
    iota = np.random.uniform(0., np.pi, Md)
    spin1 = np.random.uniform(-0.99, 0.99, Md)
    spin2 = np.random.uniform(-0.99, 0.99, Md)
    return alpha, delta, iota, spin1, spin2
# Calculating data for uniform distribution
def Uniform(Md, m_min, m_max):
    # Uniform component masses
    m1 = np.random.uniform(m_min, m_max, 10*Md)
    m2 = np.random.uniform(m_min, m_max, 10*Md)
    # Indices of samples whose total mass is lower than 100 Msun
    i_accpt = np.where((m1+m2) < (m_max+m_min))
    # Take only "Md" combinations out of the acceptable ones
    m1 = m1[i_accpt][:Md]
    m2 = m2[i_accpt][:Md]
    return m1, m2
# Calculating data for log flat distribution
def log_flat(Md, m_min, m_max):
    # generating component masses uniform in log
    # Deliberately sampling 10 times higher number of samples
    m1 = np.random.uniform(np.log(m_min), np.log(m_max), 10*Md)
    m1 = np.exp(m1)
    m2 = np.random.uniform(np.log(m_min), np.log(m_max), 10*Md)
    m2 = np.exp(m2)
    # Indices of samples whose total mass is lower than m_max
    i_accpt = np.where((m1+m2) < (m_max+m_min))
    # Take only "Md" combinations out of the acceptable ones
    m1 = m1[i_accpt][:Md]
    m2 = m2[i_accpt][:Md]
    return m1, m2
# Calculating data for power law distribution
def pow_law(Md, m_min, m_max):
    # Power to which 1st component mass is distributed, P(m1) = m1^-ALPHA
    ALPHA = 2.35
    A = m_min**(1.-ALPHA) - (m_max)**(1.-ALPHA)
    C = 2. + (m_max)**(1.-ALPHA) / A
    # Deliberately sampling 10 times higher number of samples
    m1 = np.random.uniform(1., 2., 10*Md)
    m1 = (A * (C - m1)) ** (1./(1.-ALPHA))
    m2 = np.random.uniform(m_min, m1, 10*Md)
    # Indices of samples whose total mass is lower than 100 Msun
    i_accpt = np.where((m1+m2) < (m_max+m_min))
    # Take only "Md" combinations out of the acceptable ones
    m1 = m1[i_accpt][:Md]
    m2 = m2[i_accpt][:Md]
    return m1, m2


r = {}
for RUN in r_extreme['Lower Space Cutoff'].keys():
    # Limits on space corresponding to respective RUNs
    rl = r_extreme['Lower Space Cutoff'][RUN]		# Lower limit on space
    ru = r_extreme['Upper Space Cutoff'][RUN]		# Upper limit on space
    # Defining the visible volume
    vol = 4. / 3. * np.pi * (ru ** 3. - rl ** 3.)
    # From uniform discretization of volume of space, we obtain corresponding discretization of radial distance
    r[RUN] = np.linspace(0., vol, Md)
    r[RUN] *= 3. / (4. * np.pi)
    r[RUN] += rl ** 3.
    r[RUN] **= (1. / 3.)
    # We generate array of "r" as Md samples from above distance discretization
    r[RUN] = np.random.choice(r[RUN], Md*Iterations)

# Create HDF5 file
temps = h5py.File("./../Data/fd_templates_{}{}.hdf5".format(Md, fd_temps_file_suffix), "a")
distribs = {"Uniform": Uniform, "Log_Flat": log_flat, "Power_Law": pow_law}

Dist = {}
for RUN in r_extreme['Lower Space Cutoff'].keys():
	Dist[RUN] = temps.create_dataset("Distance/{}".format(RUN), data=r[RUN])

# Random sky position and orientation
alpha, delta, iota, spin1, spin2 = a_d_i(Md)
Iota = temps.create_dataset("iota", data=iota)
temps.create_dataset("alpha", data=alpha)
temps.create_dataset("delta", data=delta)
Spin1 = temps.create_dataset("spin1", data=spin1)
Spin2 = temps.create_dataset("spin2", data=spin2)

for distrib in distribs.keys():
    print distrib
    # Group "Distribution"
    in_distrib = temps.create_group(distrib)
    # Masses M1 and M2 generated according to "Distribution"
    m1, m2 = distribs[distrib](Md, mass_min, mass_max)
    Mass1 = in_distrib.create_dataset("Mass1", data=m1)
    Mass2 = in_distrib.create_dataset("Mass2", data=m2)
    # Templates in FD to be saved for dist=100 Mpc
    temps_in_distrib = in_distrib.create_dataset("templates", (Md, len(PSD_for_det['O1_L1'])), dtype='c16')
    # Filling the frequency series values
    for m1, m2, io, i, s1, s2 in zip(Mass1, Mass2, Iota, range(Md), Spin1, Spin2):
        # Get each template
        sptilde, _ = get_fd_waveform(approximant="IMRPhenomD", mass1=m1, mass2=m2, spin1z=s1, spin2z=s2, distance=100., inclination=io, delta_f=df, f_lower=freq[0], f_final=freq[-1])
        sptilde.resize(len(PSD_for_det['O1_L1']))
        temps_in_distrib[i] = np.array(sptilde.data)
        print distrib, i

temps.close()
