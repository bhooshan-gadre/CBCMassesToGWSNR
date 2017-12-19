Obs = 'S5'
eventname = 'S5_data'

# Obs = 'S6'
# eventname = 'S6_data'

# Obs = 'O1'
# eventname = 'GW150914' 
#eventname = 'GW151226' 
#eventname = 'LVT151012'

# Obs = 'O2'
# eventname = 'GW170104'

# Obs = 'Design'


# want plots?
make_plots = 1
make_psds = 1
plottype = "png"

# In[21]:
# Standard python numerical analysis imports:
import numpy as np
from scipy import signal
from scipy.interpolate import interp1d
from scipy.signal import butter, filtfilt, iirdesign, zpk2tf, freqz
import h5py
import json
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
# LIGO-specific readligo.py 
import readligo as rl

# In[22]:
# Read the event properties from a local json file
fnjson = "BBH_events_v3.json"
try:
    events = json.load(open(fnjson,"r"))
except IOError:
    print("Cannot find resource file "+fnjson)
    print("You can download it from https://losc.ligo.org/s/events/"+fnjson)
    print("Quitting.")
    quit()
# did the user select the eventname ?
try: 
    events[eventname]
except:
    print('You must select an eventname that is in '+fnjson+'! Quitting.')
    quit()

# In[4]:
# Extract the parameters for the desired event:
event = events[eventname]
fn_H1 = event['fn_H1']              # File name for H1 data
fn_L1 = event['fn_L1']              # File name for L1 data
fs = event['fs']                    # Set sampling rate
print("Reading in parameters for event " + event["name"])
print(event)

# ## Read in the data
# We will make use of the data, and waveform template, defined above.

# In[5]:
#----------------------------------------------------------------
# Load LIGO data from a single file.
# FIRST, define the filenames fn_H1 and fn_L1, above.
#----------------------------------------------------------------
try:
    # read in data from H1 and L1, if available:
    strain_H1, time_H1, chan_dict_H1 = rl.loaddata(fn_H1, 'H1')
    strain_L1, time_L1, chan_dict_L1 = rl.loaddata(fn_L1, 'L1')
except:
    print("Cannot find data files!")
    print("You can download them from https://losc.ligo.org/s/events/"+eventname)
    print("Quitting.")
    quit()

# ## Data Gaps
# ## First look at the data from H1 and L1

# In[6]:
# both H1 and L1 will have the same time vector, so:
time = time_H1
# the time sample interval (uniformly sampled!)
dt = time[1] - time[0]

# Let's look at the data and print out some stuff:
print('time_H1: len, min, mean, max = ', len(time_H1), time_H1.min(), time_H1.mean(), time_H1.max())
print('strain_H1: len, min, mean, max = ', len(strain_H1), strain_H1.min(),strain_H1.mean(),strain_H1.max())
print('strain_L1: len, min, mean, max = ', len(strain_L1), strain_L1.min(),strain_L1.mean(),strain_L1.max())

#What's in chan_dict?  (See also https://losc.ligo.org/tutorials/)
bits = chan_dict_H1['DATA']
print("For H1, {0} out of {1} seconds contain usable DATA".format(bits.sum(), len(bits)))
bits = chan_dict_L1['DATA']
print("For L1, {0} out of {1} seconds contain usable DATA".format(bits.sum(), len(bits)))

# In[8]:
if make_psds:
    # number of sample for the fast fourier transform:
    NFFT = 4*fs
    Pxx_H1, freqs = mlab.psd(strain_H1, Fs = fs, NFFT = NFFT)
    Pxx_L1, freqs = mlab.psd(strain_L1, Fs = fs, NFFT = NFFT)

    # We will use interpolations of the ASDs computed above for whitening:
    psd_H1 = interp1d(freqs, Pxx_H1)
    psd_L1 = interp1d(freqs, Pxx_L1)

    # Here is an approximate, smoothed PSD for H1 during O1, with no lines. We'll use it later.    
    # Pxx = (1.e-22*(18./(0.1+freqs))**2)**2+0.7e-23**2+((freqs/2000.)*4.e-23)**2
    # psd_smooth = interp1d(freqs, Pxx)

if make_plots:
    # plot the ASDs, with the template overlaid:
    f_min = 20.
    f_max = 2000.
    plt.figure(figsize=(10,8))
    plt.loglog(freqs, np.sqrt(Pxx_L1),'g',label='L1 strain')
    plt.loglog(freqs, np.sqrt(Pxx_H1),'r',label='H1 strain')
    # plt.loglog(freqs, np.sqrt(Pxx),'k',label='H1 strain, O1 smooth model')
    plt.axis([f_min, f_max, 1e-24, 1e-19])
    plt.grid('on')
    plt.ylabel('ASD (strain/rtHz)')
    plt.xlabel('Freq (Hz)')
    plt.legend(loc='upper center')
    plt.title('Advanced LIGO strain data near '+eventname)
    # plt.savefig(eventname+'_ASDs.'+plottype)

plt.show()

# Save the PSD File
File_H1 = open(('asd_%s_H1.txt' %Obs), 'w')
File_L1 = open(('asd_%s_L1.txt' %Obs), 'w')
np.savetxt(File_H1, np.array([freqs, np.sqrt(Pxx_H1)]).transpose(), delimiter=',', newline='\n')
np.savetxt(File_L1, np.array([freqs, np.sqrt(Pxx_L1)]).transpose(), delimiter=',', newline='\n')
File_H1.close()
File_L1.close()
