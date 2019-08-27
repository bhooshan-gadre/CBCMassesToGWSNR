# The redshift values are sampled from the Metropolis algorithm used in make_inj.py code
# Corresponding comoving distances (d_c) are calculated
# Test: a PDF of (1+z)*d_c**2 should be flat for uniform distribution in comoving volume

import numpy
import lal
import matplotlib.pyplot as plt

# seed the random number generator to get reproducible results
numpy.random.seed(7)

# Taken from rates&pop repo (LSC-git)
def draw_redshift(zmax, omega, distr_unif_in_Vc=True):
    '''
    Yields a random redshift from a cosmologically-correct distribution.
    Uses Metropolis algorithm to draw from the desired pdf.
    '''

    def pdf(z):
        '''
        This redshift pdf yields a uniform distribution
        in comoving volume divided by (1+z).
        '''
        # FIXME: XLALUniformComovingVolumeDensity() currently implements
        # the factor of 1/(1+z) that converts to source-frame time.
        # If this changes, modify the code below.
        return lal.UniformComovingVolumeDensity(z, omega)
        #return lal.UniformComovingVolumeDensity(z, omega) / (1.0 + z)


    while True:
        if distr_unif_in_Vc:
            z0 = numpy.random.uniform(0.0, zmax)
            p0 = pdf(z0)
            # acceptance rate is 50% so take every 10th
            # draw from distribution to avoid repeating
            # the same value too often
            for _ in range(10):
                z = numpy.random.uniform(0.0, zmax)
                p = pdf(z)
                if p > p0 or numpy.random.random() < p / p0:
                    z0 = z
                    p0 = p
            yield z0
        else:
            z0 = numpy.random.uniform(0.0, zmax)
            yield z0

# def draw_redshift_metropolis(zmax, omega, distr_unif_in_Vc=True):
#     '''
#     Yields a random redshift from a cosmologically-correct distribution.
#     Uses Metropolis algorithm to draw from the desired pdf.
#     '''

#     def pdf(z):
#         '''
#         This redshift pdf yields a uniform distribution
#         in comoving volume divided by (1+z).
#         '''
#         # FIXME: XLALUniformComovingVolumeDensity() currently implements
#         # the factor of 1/(1+z) that converts to source-frame time.
#         # If this changes, modify the code below.
#         return lal.UniformComovingVolumeDensity(z, omega)
#         #return lal.UniformComovingVolumeDensity(z, omega) / (1.0 + z)


#     while True:
#         if distr_unif_in_Vc:
#             z0 = numpy.random.uniform(0.0, zmax)
#             p0 = pdf(z0)
#             # acceptance rate is 50% so take every 10th
#             # draw from distribution to avoid repeating
#             # the same value too often
#             z = numpy.random.uniform(0.0, zmax)
#             p = pdf(z)
#             if numpy.random.random() < min(1, p / p0):
#                 z0 = z
#                 p0 = p
#             yield z0
#         else:
#             z0 = numpy.random.uniform(0.0, zmax)
#             yield z0

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

random_redshift = iter(draw_redshift(max_redshift, omega, distr_unif_in_Vc=Uniform_in_Vc))
random_redshift_metropolis = iter(draw_redshift(max_redshift, omega, distr_unif_in_Vc=Uniform_in_Vc))

z = []
d_c = []
for i in range(100000):
	zz = next(random_redshift)
	z.append(zz)
	d_c.append(lal.LuminosityDistance(omega, zz) / (1 + zz))

z = numpy.array(z)
d_c = numpy.array(d_c)


# # custom code
# z2 = []
# d_c2 = []
# for i in range(100000):
# 	zz2 = next(random_redshift_metropolis)
# 	z2.append(zz2)
# 	d_c2.append(lal.LuminosityDistance(omega, zz2) / (1 + zz2))

# z2 = numpy.array(z2)
# d_c2 = numpy.array(d_c2)

# Plotting
power = 1
N = 25
plt.hist((1+z)**power*d_c**2, N, histtype='step', density=True, label='lsc-git')
# plt.hist((1+z2)**power*d_c2**2, N, histtype='step', density=True, label='custom')
# plt.legend()
plt.xlabel('(1+z) * $d_c^2$')
plt.ylabel('Probability')
plt.show()
