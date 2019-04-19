from __future__ import division
import logging
from collections import OrderedDict
import numpy as np
import scipy
import scipy.optimize

import lal
import lalsimulation


##################################################

fs = 2048
dt = 1./fs
maxlen = 256
tlen = fs*maxlen
flen = tlen/2 + 1
df = 1./maxlen
flow = df

CANONICAL_SNR = 5.

##################################################
# set cosmological parameters
# From Planck2015, Table IV
OMEGA = lal.CreateCosmologicalParametersAndRate().omega
lal.SetCosmologicalParametersDefaultValue(OMEGA)
OMEGA.h = 0.679
OMEGA.om = 0.3065
OMEGA.ol = 0.6935
OMEGA.ok = 1.0 - OMEGA.om - OMEGA.ol
OMEGA.w0 = -1.0
OMEGA.w1 = 0.0
OMEGA.w2 = 0.0


# default inspiral waveform parameters
# face-on 50-50 Msolar inspiral at 100 Mpc distance
DEFAULT_PARAMS = OrderedDict([('approximant', None), ('distance', 100e6), ('m1', 50.), ('m2', 50.), ('S1x', 0.0), ('S1y', 0.0), ('S1z', 0.0), ('S2x', 0.0), ('S2y', 0.0), ('S2z', 0.0), ('inclination', 0.0), ('f_ref', 0.0), ('phiRef', 0.0), ('longAscNodes', 0.0), ('eccentricity', 0.0), ('meanPerAno', 0.0), ('deltaF', None), ('f_min', None), ('f_max', None), ('LALpars', None)])

DEFAULT_APPROXIMANT_BNS = 'TaylorF2'
DEFAULT_APPROXIMANT_BBH = 'IMRPhenomD'

##################################################

def _get_waveform_params(**kwargs):
    params = OrderedDict(DEFAULT_PARAMS)
    params.update(**kwargs)
    # use waveform approximant appropriate to type
    if not params['approximant']:
        if params['m1'] >= 5 and params['m2'] >= 5:
            params['approximant'] = DEFAULT_APPROXIMANT_BBH
        else:
            params['approximant'] = DEFAULT_APPROXIMANT_BNS
        # print params
    return params

def gen_waveform(freq, z=0, omega=OMEGA, **params):
    """Generate frequency-domain inspiral waveform

    `freq` should be an array of frequency points at which the
    waveform should be interpolated.  Returns a tuple of
    (h_tilde^plus, h_tilde^cross) real-valued (amplitude only) arrays.

    The waveform is generated with
    lalsimulation.SimInspiralChooseFDWaveform().  Keyword arguments
    are used to update the default waveform parameters (1.4/1.4
    Msolar, optimally-oriented, 100 Mpc, (see DEFAULT_PARAMS macro)).
    The mass parameters ('m1' and 'm2') should be specified in solar
    masses and the 'distance' parameter should be specified in
    parsecs**.  Waveform approximants may be given as string names
    (see `lalsimulation` documentation for more info).  If the
    approximant is not specified explicitly, DEFAULT_APPROXIMANT_BNS
    waveform will be used if either mass is less than 5 Msolar and
    DEFAULT_APPROXIMANT_BBH waveform will be used otherwise.

    If a redshift `z` is specified (with optional `omega`), it's
    equivalent distance will be used (ignoring any `distance`
    parameter provided) and the masses will be redshift-corrected
    appropriately.  Otherwise no mass redshift correction will be
    applied.

    For example, to generate a 20/20 Msolar BBH waveform:

    >>> hp,hc = waveform.gen_waveform(freq, 'm1'=20, 'm2'=20)

    **NOTE: The requirement that masses are specified in solar masses
    and distances are specified in parsecs is different than that of
    the underlying lalsimulation method which expects mass and
    distance parameters to be in SI units.

    """
    iparams = _get_waveform_params(**params)

    # if redshift specified use that as distance and correct
    # appropriately, ignoring any distance specified in params.
    if z != 0:
        iparams['distance'] = lal.LuminosityDistance(omega, z) * 1e6
        iparams['m1'] *= 1.0 + z
        iparams['m2'] *= 1.0 + z

    # convert to SI units
    iparams['distance'] *= lal.PC_SI
    iparams['m1'] *= lal.MSUN_SI
    iparams['m2'] *= lal.MSUN_SI
    iparams['approximant'] = lalsimulation.SimInspiralGetApproximantFromString(iparams['approximant'])

    iparams['deltaF'] = freq[1] - freq[0]
    iparams['f_min'] = freq[0]
    # FIXME: the max frequency in the generated waveform is not always
    # greater than f_max, so as a workaround we generate over the full
    # band.  Probably room for speedup here
    # iparams['f_max'] = freq[-1]
    iparams['f_max'] = 10000
    # print iparams

    # logging.debug('waveform params = {}'.format(iparams))

    # generate waveform
    h = lalsimulation.SimInspiralChooseFDWaveform(**iparams)
    # print h

    freq_h = h[0].f0 + np.arange(len(h[0].data.data)) * h[0].deltaF

    def interph(h):
        "interpolate amplitude of h array"
        # FIXME: this only interpolates/returns the amplitude, and not
        # the full complex data (throws out phase), because this was
        # not working:
        # hir = scipy.interpolate.interp1d(freq_h, np.real(h.data.data))(freq)
        # hii = scipy.interpolate.interp1d(freq_h, np.imag(h.data.data))(freq)
        # return hir + 1j * hii
        return scipy.interpolate.interp1d(freq_h, np.absolute(h.data.data))(freq)

    hi = map(interph, h)

    return hi

def SNR(freq, psd, h):
    """SNR for given PSD and waveform

    @returns SNR as a float

    """
    rho = np.sqrt(np.trapz(4*(h**2/psd), freq))
    return rho

def SNR_z(freq, psd, z, omega=OMEGA, **params):
    """SNR of waveform with specified parameters at redshift z

    Remaining keyword arguments are interpreted as waveform generation
    parameters (see waveform.gen_waveform()).

    @returns SNR as a float

    """
    hplus, hcross = gen_waveform(freq, z=z, omega=omega, **params)
    # FIXME: should this be the average of hplus and hcross??
    rho = SNR(freq, psd, hplus)
    print rho
    return rho

def horizon_redshift(freq, psd, omega=OMEGA, **params):
    """Detector horizon redshift

    For the given detector noise PSD and waveform parameters return
    the redshift at which the detection SNR would equal to the
    CANONICAL_SNR.

    Remaining keyword arguments are interpreted as waveform generation
    parameters (see waveform.gen_waveform()).

    @returns redshift as a float

    """
    assert len(freq) == len(psd)

    # we will performa a Brent method optimization to find the root of
    # the following function, thereby returning the z at which
    # SNR = CANONICAL_SNR:
    def opt_SNR_z(z):
        return SNR_z(freq, psd, z, omega, **params) - CANONICAL_SNR

    zmin = 1e-8 # must be less than horizon
    zmax = 10.0 # must be greater than horizon
    print 'opt_SNR_z(zmin): {}'.format(opt_SNR_z(zmin))
    print 'opt_SNR_z(zmax): {}'.format(opt_SNR_z(zmax))

    # A ValueError is returned if the ranges do not cover zero.  This
    # is probably because the zmax is not large enough, so bump the
    # max and try again.
    # FIXME: better checking of this (pre-check?)
    try:
        z_hor = scipy.optimize.brentq(opt_SNR_z, zmin, zmax)
    except ValueError:
        zmax = 100.0
        print "increasing zmax => {}...".format(zmax)
        print 'opt_SNR_z(zmax): {}'.format(opt_SNR_z(zmax))
        z_hor = scipy.optimize.brentq(opt_SNR_z, zmin, zmax)

    return z_hor

def horizon(freq, psd, omega=OMEGA, **params):
    """Detector horizon distance in Mpc

    See horizon_redshift().

    @returns distance in Mpc as a float

    """
    zhor = horizon_redshift(freq, psd, omega=omega, **params)
    return lal.LuminosityDistance(omega, zhor)

##################################################
# RUN="O1"
RUN="O2"

WPATH="/home/shreejit.jadhav/WORK"
# WPATH="/home/shreejit/Dropbox/Academic/WORK"

H1L1_PSD="{0}/CBCMassesToGWSNR/Data/PSDs/{1}/H1L1_{1}_PSD.txt".format(WPATH, RUN)

psd = np.genfromtxt(H1L1_PSD, delimiter=" ")
freqs = psd[:, 0][4000:]
psd = psd[:, 1][4000:]

m=50
z_hor = horizon_redshift(freqs, psd, omega=OMEGA, m1=m, m2=m)
D_hor = lal.LuminosityDistance(OMEGA, z_hor)

print z_hor, D_hor



# # First we create an injection sim file each for upper and lower cutoffs

# # RUN="O1"
# RUN="O2"
# LIMIT="Upper"
# # LIMIT="Lower"

# if [ $LIMIT == "Lower" ]
# then
# MASS=5.
# DMIN=10000
# DMAX=90000
# echo "1"
# elif [ $LIMIT == "Upper" ]
# then
# MASS=50.
# DMIN=1000000
# DMAX=9000000
# echo "2"
# fi

# echo $MASS

# if [ $RUN == "O1" ]
# then
# RUN_START_TIME=1126051217
# elif [ $RUN == "O2" ]
# then
# RUN_START_TIME=1164556817
# fi


# DATA_DURATION=30000
# RUN_END_TIME=$((RUN_START_TIME+DATA_DURATION))
# echo $RUN_END_TIME

# # WPATH="/home/shreejit.jadhav/WORK"
# WPATH="/home/shreejit/Dropbox/Academic/WORK"


# INJ_FLOW=15

# H1L1_PSD="${WPATH}/CBCMassesToGWSNR/Data/PSDs/${RUN}/H1L1_${RUN}_PSD.txt"

# INJ_FILE="${WPATH}/CBCMassesToGWSNR/Data/new_inj_${RUN}_${DATA_DURATION}_${LIMIT}.xml"
# TIME_STEP=10
# TIME_INTERVAL=1


# # upper limit

# lalapps_inspinj -o ${INJ_FILE} --verbose \
#     --d-distr uniform --min-distance ${DMIN} --max-distance ${DMAX} \
#     --m-distr log --min-mass1 ${MASS} --min-mass2 ${MASS} --max-mass1 ${MASS} --max-mass2 ${MASS} \
#     --t-distr uniform --time-step ${TIME_STEP} --time-interval ${TIME_INTERVAL} \
#     --gps-start-time ${RUN_START_TIME} \
#     --gps-end-time ${RUN_END_TIME} \
#     --f-lower 10 \
#     --waveform IMRPhenomPv2threePointFivePN \
#     --l-distr fixed --longitude 4.52389342 --latitude 0.81681409 --i-distr fixed --fixed-inc 0. --disable-spin \
#     --ligo-start-freq ${INJ_FLOW} --taper-injection startend \
#     --ligo-psd ${H1L1_PSD}

# ######################################################

# # Space limits from injection files
# import pycbc.inject as inject
# import pycbc.psd as pp
# import pycbc.filter as pf
# import numpy as np

# flow = 25
# fs = 2048
# dt = 1./fs
# maxlen = 256
# tlen = fs*maxlen
# flen = tlen/2 + 1
# df = 1./maxlen

# r_extreme = {"Lower": {}, "Upper": {}}
# # WPATH="/home/shreejit.jadhav/WORK"
# WPATH="/home/shreejit/Dropbox/Academic/WORK"
# DATA_DURATION=30000

# # for RUN in ["O1", "O2"]:
# # 	for LIMIT in ["Upper", "Lower"]:
# RUN = "O1"
# # RUN = "O2"
# LIMIT = "Lower"
# # LIMIT = "Upper"
# # ASD file
# asd_files = dict(H1_FILE='{0}/CBCMassesToGWSNR/Data/PSDs/{1}/H1_{1}_PSD.txt'.format(WPATH, RUN),
#         L1_FILE='{0}/CBCMassesToGWSNR/Data/PSDs/{1}/L1_{1}_PSD.txt'.format(WPATH, RUN))

# ifos = ['H1', 'L1']

# psds = {}
# for det in ifos:
#     psds['{0}'.format(det)] = pp.from_txt(asd_files['{0}_FILE'.format(det)], flen, df, flow, is_asd_file=False)

# # injection_file = '/home/shreejit.jadhav/WORK/CBCMassesToGWSNR/Data/inj_trial_O1_300000_2.xml'
# injection_file = "{0}/CBCMassesToGWSNR/Data/new_inj_{1}_{2}_{3}.xml".format(WPATH, RUN, DATA_DURATION, LIMIT)
# inj = inject.InjectionSet(injection_file)

# snrs = {}
# m1 = []
# m2 = []
# d = []
# for det in ifos:
#     snrs[det] = []

# for det in ifos:
# 	for inj_id, params in enumerate(inj.table):
# 		if inj_id % 100 == 0:
# 			print "{2}: {0} / {1}".format(inj_id, len(inj.table), det)
# 		h = inj.make_strain_from_inj_object(params, dt, det, flow)
# 		h.resize(tlen)
# 		s = pf.sigma(h, psds[det], flow)
# 		snrs[det].append(s)
# 		if det=="H1":
# 			m1.append(params.mass1)
# 			m2.append(params.mass2)
# 			d.append(params.distance)
# 		# print '\r{0} snr = {1} for injection {2}, m1: {3}, m2: {4}, d:{5}'.format(det, s, inj_id, params.mass1, params.mass2, params.distance)

# m1 = np.array(m1)
# m2 = np.array(m2)
# d = np.array(d)
# for det in ifos:
#     snrs[det] = np.array(snrs[det])

# # print params.inclination, d[ii], m1[ii], m2[ii], snrs['L1'][ii], snrs['H1'][ii] 
# snr = (snrs['H1'] ** 2. + snrs['L1'] ** 2.) ** 0.5
# # sorting
# i_sort = np.argsort(snr)
# m1 = m1[i_sort]
# m2 = m2[i_sort]
# d = d[i_sort]
# snr = snr[i_sort]
# # choose space lim
# if LIMIT=="Upper":
# 	Side, snr_lim = "left", 5
# elif LIMIT=="Lower":
# 	Side, snr_lim = "right", 60
# i_lim = np.searchsorted(snr, snr_lim, side=Side)
# r_extreme[LIMIT][RUN] = d[i_lim]
# message = "{4}: SNR: {0} with mass pair ({1}, {1}) gives {2} limit on space as {3}".format(snr[i_lim], m1[0], LIMIT, r_extreme[LIMIT][RUN], RUN)
# print message
# print r_extreme
# np.save('{}/CBCMassesToGWSNR/Data/new_{}_{}'.format(WPATH, RUN, LIMIT), message)


# # O1: 
# # SNR: 5.00039388468 with mass pair (50.0, 50.0) gives Upper limit on space as 3882.545
# # SNR: 60.0040936381 with mass pair (5.0, 5.0) gives Lower limit on space as 59.46056
# # O2: 
# # SNR: 5.00055571883 with mass pair (50.0, 50.0) gives Upper limit on space as 2828.342



# ######################################################

# from function_SNR import *
# import matplotlib.pyplot as plt

# # Taking r = 100 Mpc
# r = np.array([100.*Mpc])
# r = {'S6':r, 'O1':r, 'O2':r, 'Design':r}
# # limiting masses
# # Change: to cover full SNR window and avoid "leaking", we find space limits such that lightest binary at closest distance will give highest SNR and heaviest binary at farthest distance will give lowest SNR
# mass_lim = {'Lower Space Cutoff': 5., 'Upper Space Cutoff': 50.}
# # SNR cutoffs
# cut = {'Lower Space Cutoff': 60., 'Upper Space Cutoff': 5.}
# # dict of SNRs which will contain SNR dicts against cutoffs. These SNR dicts in turn contain SNR arrays corresponding to diff RUNs
# SNRs = {'Lower Space Cutoff': {}, 'Upper Space Cutoff': {}}
# # limits on space
# r_extreme = {'Lower Space Cutoff': {}, 'Upper Space Cutoff': {}}

# # x_l = 0.
# # x_u = 3000.
# # y_u = 300.
# # tolrnc = 0.2
# # off_x = 50.
# # off_y = 4.
# # RUN = 'S6'

# for cutoff in cut.keys():
# 	# creating injection arrays with cutoff mass pairs
# 	M1 = np.array([mass_lim[cutoff]]) * Msun
# 	M2 = np.array([mass_lim[cutoff]]) * Msun
# 	# finding SNRs for respective cutoffs for diff RUNs
# 	SNRs[cutoff]  = find_pycbc_SNR(M1, M2, r, [0.], [0.], [0.])
# 	# find extremity in space corresponding to the SNR cutoff
# 	for RUN in SNRs[cutoff].keys():
# 		r_extreme[cutoff][RUN] = SNRs[cutoff][RUN] * r[RUN] / cut[cutoff] / Mpc

# print r_extreme
# #########################################
# # 				Output:

# # With find_pycbc_SNR
# # {'Upper Space Cutoff':
# # {'Design': array([3279.12309453]),
# # 'S6': array([153.40926056]),
# # 'O1': array([901.50600027]),
# # 'O2': array([1153.03794345])},
# # 'Lower Space Cutoff':
# # {'Design': array([1587.82130486]),
# # 'S6': array([82.77513475]),
# # 'O1': array([454.18511569]),
# # 'O2': array([576.59669512])}}

# # New
# # {'Upper Space Cutoff': {'Design': array([19063.31247625]), 'S6': array([994.86000718]), 'O1': array([5260.97949502]), 'O2': array([6597.62537123])}, 'Lower Space Cutoff': {'Design': array([273.38516829]), 'S6': array([12.80446162]), 'O1': array([72.82472261]), 'O2': array([91.44922879])}}


# # With find_simple_SNR
# # {'Upper Space Cutoff': 
# # {'Design': array([3452.14263401]),
# #  'O1': array([964.90328158]),
# #  'O2': array([1224.40665554]),
# #  'S6': array([169.03614154])},
# # 'Lower Space Cutoff': 
# # {'Design': array([1056.05949735]),
# #  'O1': array([179.04932493]),
# #  'O2': array([267.89876495]),
# #  'S6': array([1.56752347])}
# #########################################

# # Plotting
# fig = plt.figure(figsize=(15,10))
# # for all RUNs
# for RUN, colr in zip(r_extreme["Lower Space Cutoff"].keys(), ["g", "r", "b", "m"]):
# 	for cutoff in cut.keys():
# 		# SNR vs r curves
# 		plt.plot(r[RUN]/Mpc, SNRs[cutoff][RUN], colr+"-", label=RUN + ": " + cutoff)
# 	plt.fill_between(r[RUN]/Mpc, 0, 100., where=(r[RUN]/Mpc<r_extreme["Upper Space Cutoff"][RUN]/Mpc) & (r[RUN]/Mpc>r_extreme["Lower Space Cutoff"][RUN]/Mpc), edgecolor=colr, facecolor=colr, alpha=0.7, label=RUN)
# plt.xlim(0., 4000.)
# plt.ylim(0., 80.)
# # plt.xscale("log")
# # plt.yscale("log")
# plt.grid(True)
# plt.xlabel('Distance (Mpc)')
# plt.ylabel('SNR')
# plt.title("Limits on observable space for different observational runs of LIGO")
# plt.legend(loc="upper right")
# plt.savefig('./../Plots/Limits_on_space_3March', format='png')
# plt.show()
