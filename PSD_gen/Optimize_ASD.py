import numpy as np
from scipy import interpolate

def opt(ASDdata):
    asd = interpolate.interp1d(np.log(ASDdata[:, 0]), np.log(ASDdata[:, 1]))
    # Calculating ASD values at respective elements of frequency array
    asd_at_freq = asd(np.log(freq))
    asd_at_freq = np.exp(asd_at_freq)
    return asd_at_freq

# Which run?
# RUN = 'S5'
RUN = 'S6'
# RUN = 'O1'
# RUN = 'O2'
# RUN = 'Design'

# Defining a dictionary for detector ASDs and LIGO runs
# These names will be used to pick ASD files
detector = {0:'S6_L1', 1:'S6_H1', 3:'O1_L1', 4:'O1_H1', 5:'O2_L1', 6:'O2_H1', 8:'Design_L', 9:'Design_H', 11:'S5_L1', 12:'S5_H1'}
# Dictionary to tell which files to pick for a particular run
# not considering Virgo detector in O2
for_run = {'S5':[11, 12], 'S6':[0, 1], 'O1':[3, 4], 'O2':[5, 6], 'Design':[8, 9]}

# Necessary freq range
df = 0.125
freq = np.arange(20., 1500.+df, df)

##### Currently we are considering only one PSD for Design specification #####

if RUN == 'Design':
    # Read the PSD file (heavy)
    dat = np.genfromtxt('/home/shreejit/Design_PSD.txt')
    # Convert to ASD file (still heavy)
    dat[:, 1] = np.sqrt(dat[:, 1])
    # Picking out indices of concerned freq range
    # indices = []
    ASD = [[], []]
    for i in range(len(dat[:,0])):
        if float(dat[i,0]/df).is_integer() and float(dat[i,0])>=20. and float(dat[i,0])<=1500.:
            # indices.append(i)
            ASD[0].append(dat[i,0])
            ASD[1].append(dat[i,1])

    # numpy array of freq and ASD
    ASD = np.array(ASD).transpose()
    # Saving the new file
    f = open('/home/shreejit/ASD/asd_Design.txt', 'w')
    np.savetxt(f, ASD, delimiter=',', newline='\n')
    f.close()

else:
    for i in for_run[RUN]:
        ASDdata = np.genfromtxt('/home/shreejit/asd_%s.txt' %detector[i], delimiter=',')
        opt_asd = np.array([freq, opt(ASDdata)]).transpose()
        f2 = open('/home/shreejit/ASD/asd_%s.txt' %detector[i], 'w')
        np.savetxt(f2, opt_asd, delimiter=',', newline='\n')
        f2.close()