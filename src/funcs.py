from __future__ import division
import logging
from collections import OrderedDict

import scipy
import scipy.optimize
import scipy.interpolate
import json
import h5py as hp
import numpy as np
import matplotlib.pyplot as plt

from pycbc.waveform import get_fd_waveform
import pycbc.filter.matchedfilter as mf
import pycbc.types.frequencyseries as pf

from tqdm import tqdm
import datetime

import lal
import lalsimulation


##################################################

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

gps_time = 1129383017 # fixed time for detector orientation

# default inspiral waveform parameters
# face-on 50-50 Msolar inspiral at 100 Mpc distance
DEFAULT_PARAMS = OrderedDict([('approximant', None), ('distance', 100e6), ('m1', 50.), ('m2', 50.), ('S1x', 0.0), ('S1y', 0.0), ('S1z', 0.0), ('S2x', 0.0), ('S2y', 0.0), ('S2z', 0.0), ('inclination', 0.0), ('f_ref', 0.0), ('phiRef', 0.0), ('longAscNodes', 0.0), ('eccentricity', 0.0), ('meanPerAno', 0.0), ('deltaF', None), ('f_min', None), ('f_max', None), ('LALpars', None)])

detectors = {
        'H1' : lal.CachedDetectors[lal.LHO_4K_DETECTOR],
        'L1' : lal.CachedDetectors[lal.LLO_4K_DETECTOR]
        }

DEFAULT_APPROXIMANT_BNS = 'TaylorF2'
DEFAULT_APPROXIMANT_BBH = 'IMRPhenomD'

##################################################

# WPATH="/home/shreejit.jadhav/WORK"
WPATH="/home/shreejit/Dropbox/Academic/WORK"

def _get_waveform_params(**kwargs):
    params = OrderedDict(DEFAULT_PARAMS)
    params.update(**kwargs)
    # use waveform approximant appropriate to type
    if not params['approximant']:
        if params['m1'] >= 5 and params['m2'] >= 5:
            params['approximant'] = DEFAULT_APPROXIMANT_BBH
        else:
            params['approximant'] = DEFAULT_APPROXIMANT_BNS
        # print(params)
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
    # print(iparams)

    # logging.debug('waveform params = {}'.format(iparams))

    # generate waveform
    h = lalsimulation.SimInspiralChooseFDWaveform(**iparams)
    # print(h)

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
    # print(rho)
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
    # print('opt_SNR_z(zmin): {}'.format(opt_SNR_z(zmin)))
    # print('opt_SNR_z(zmax): {}'.format(opt_SNR_z(zmax)))

    # A ValueError is returned if the ranges do not cover zero.  This
    # is probably because the zmax is not large enough, so bump the
    # max and try again.
    # FIXME: better checking of this (pre-check?)
    try:
        z_hor = scipy.optimize.brentq(opt_SNR_z, zmin, zmax)
    except ValueError:
        zmax = 100.0
        print("increasing zmax => {}...".format(zmax))
        print('opt_SNR_z(zmax): {}'.format(opt_SNR_z(zmax)))
        z_hor = scipy.optimize.brentq(opt_SNR_z, zmin, zmax)

    return z_hor

def horizon(freq, psd, omega=OMEGA, **params):
    """Detector horizon distance in Mpc

    See horizon_redshift().

    @returns distance in Mpc as a float

    """
    zhor = horizon_redshift(freq, psd, omega=omega, **params)
    return lal.LuminosityDistance(omega, zhor)

# Generating random sky position and orientation (arrays)
def sample_params(Md):
    alpha = np.random.uniform(0., 2.*np.pi, Md)
    delta = np.arcsin(np.random.uniform(-1.0, 1.0, Md))
    iota = np.arccos(np.random.uniform(-1.0, 1.0, Md))
    spin1 = np.random.uniform(-0.99, 0.99, Md)
    spin2 = np.random.uniform(-0.99, 0.99, Md)
    pol = np.random.uniform(0.0, 2.0 * np.pi, Md)
    # Sample random arrival times for injections
    # t = lal.LIGOTimeGPS(gps_start_time)
    # dt = np.random.uniform(1., gps_end_time-gps_start_time-1., Md)
    # T = t + dt
    # T = lal.GreenwichMeanSiderealTime(T)
    # only one fixed time, so const antenna pattern response to sky location
    T = np.ones(Md) * gps_time
    return alpha, delta, iota, spin1, spin2, pol, T

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

# Redshift sampling with uniform density in comoving vol (based on the one used by rates&pop group of LIGO)
def sample_redshifts(zmax, Md, omega=OMEGA):
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

    Z = []
    for _ in range(Md):
        z0 = np.random.uniform(0.0, zmax)
        p0 = pdf(z0)
        # acceptance rate is 50% so take every 10th
        # draw from distribution to avoid repeating
        # the same value too often
        for _ in xrange(10):
            z = np.random.uniform(0.0, zmax)
            p = pdf(z)
            if p > p0 or np.random.random() < p / p0:
                z0 = z
                p0 = p
        Z.append(z0)
    return np.array(Z)

# Convert redshift array to luminosity distance array
def convert_to_dlum(z, omega=OMEGA):
    Dlum = []
    for rs in z:
        Dlum.append(lal.LuminosityDistance(omega, rs))
    return np.array(Dlum)

# Find the effective distance of a merger from its sky location, inclination etc
def find_eff_d(dlum, alpha, delta, pol, iota, T_inj):
    d_eff_dict = {}
    # calculate and set detector-specific columns
    for det_site, det in detectors.iteritems():
        d_eff_dict[det_site] = []
        for i in range(len(dlum)):
            fp, fc = lal.ComputeDetAMResponse(det.response, alpha[i], delta[i], pol[i], T_inj[i])
            cosi = np.cos(iota[i])
            deff = dlum[i] * ((0.5 * (1.0 + cosi**2) * fp)**2 + (cosi * fc)**2)**-0.5
            d_eff_dict[det_site].append(deff)
        d_eff_dict[det_site] = np.array(d_eff_dict[det_site])
    return d_eff_dict

def ChirpMass(Mass1, Mass2):
    Mtot = Mass1 + Mass2
    MChirp = Mass1 * Mass2
    MChirp **= 0.6
    MChirp /= Mtot**0.2
    return MChirp

# plots the histogram of vals for those with corresponding snrs falling in SNRmin and SNRmax
def Plot_Hist(vals, Run, SNRs, SNRmax, Thresh_Color, nbins=50, Lwid=1):
    for SNRmin in Thresh_Color:
        i_accpt = np.where((SNRs > SNRmin) & (SNRs < SNRmax))[0]
        param_vals = vals[i_accpt]
        # histogram of the parameter
        plt.hist(param_vals, nbins, color=Thresh_Color[SNRmin], linewidth=Lwid, \
                         histtype='step', density=True, label='SNR > {}'.format(SNRmin))
    return param_vals

def Plot_events(Obs, Run, param, Distrib, repopath, dtTime, Mass_min, Mass_max, N, Thresh_Color, Lwid=1, fsiz=13, Savefig=False):
    yr_of_run = {'O1': '15', 'O2': '17'}
    # Observed events so far
    for event in Obs['data'].keys():
        # Selecting only the observations done in resp runs
        if event[2:4] == yr_of_run[Run]:
            # Calculate the resp param
            if param == 'Chirp Mass':
                val = ChirpMass(Obs['data'][event]['mass1']['best'], \
                                Obs['data'][event]['mass2']['best'])
            elif param == 'Total Mass':
                val = Obs['data'][event]['mass1']['best']
                val += Obs['data'][event]['mass2']['best']
            elif param == 'q':
                val = Obs['data'][event]['mass1']['best']
                val /= Obs['data'][event]['mass2']['best']
                if val>1.: val = 1./val
            for thresh in np.sort(Thresh_Color.keys()):
                if Obs['data'][event]['snr_pycbc']['best'] > thresh:
                    colr = Thresh_Color[thresh]
            # Vertical line
            plt.axvline(x=val, ymin=0., ymax = 1., \
                        color=colr, linewidth=Lwid, alpha=0.5)
            plt.tick_params(axis='both', which='major', labelsize=fsiz)
    plt.xlabel(param, fontsize=fsiz, fontweight='bold')
    plt.ylabel('Probability of Detection', fontsize=fsiz, fontweight='bold')
    plt.title('Run: %s, Distribution: %s' %(Run, Distrib), fontsize=fsiz, fontweight='bold')
    plt.legend(fontsize=fsiz)
    if Savefig:
        plt.savefig('{}/Plots/{}_PDF_{}-{}-{}_{}_{}_{}.png'.\
                    format(repopath, dtTime, Run, param.replace(' ', '_'), Distrib, int(Mass_min), int(Mass_max), N))