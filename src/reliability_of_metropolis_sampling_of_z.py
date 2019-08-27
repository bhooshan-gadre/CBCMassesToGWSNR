# The redshift values are sampled from the Metropolis algorithm used in make_inj.py code
# Corresponding comoving distances (d_c) are calculated
# Test: a PDF of (1+z)*d_c**2 should be flat for uniform distribution in comoving volume

import numpy
import lal
import matplotlib.pyplot as plt

# seed the random number generator to get reproducible results
numpy.random.seed(7)

def pdf(z, omega):
        '''
        This redshift pdf yields a uniform distribution
        in comoving volume divided by (1+z).
        '''
        # FIXME: XLALUniformComovingVolumeDensity() currently implements
        # the factor of 1/(1+z) that converts to source-frame time.
        # If this changes, modify the code below.
        return lal.UniformComovingVolumeDensity(z, omega)

def sample_redshifts_ligo_method(zmax, Md, omega=OMEGA):
    '''
    Yields a random redshift from a cosmologically-correct distribution.
    Uses Metropolis algorithm to draw from the desired pdf.
    '''

    Z = []
    for _ in range(Md):
        z0 = np.random.uniform(0.0, zmax)
        p0 = pdf(z0, omega)
        # acceptance rate is 50% so take every 10th
        # draw from distribution to avoid repeating
        # the same value too often
        for _ in xrange(10):
            z = np.random.uniform(0.0, zmax)
            p = pdf(z, omega)
            if p > p0 or np.random.random() < p / p0:
                z0 = z
                p0 = p
        Z.append(z0)
    return np.array(Z)

max_redshift = 0.7
Uniform_in_Vc = True
# From Planck2015, Table IV
omega = lal.CreateCosmologicalParametersAndRate().omega
lal.SetCosmologicalParametersDefaultValue(omega)
omega.h = 0.679
omega.om = 0.3065
omega.ol = 0.6935
omega.ok = 1.0 - omega.om - omega.ol
omega.w0 = -1.0
omega.w1 = 0.0
omega.w2 = 0.0

zz = sample_redshifts_ligo_method(3, 100000)
dc = []
dvdz = []
for z in zz:
    dc.append(lal.LuminosityDistance(OMEGA, z) / (1 + z))
    dvdz.append(pdf(z, OMEGA))
dc = np.array(dc)
dvdz = np.array(dvdz)


plt.hist(zz, 35)
plt.plot(zz, dvdz/5e7, 'o')
plt.show()
