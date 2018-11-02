import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# getting the values of a specific parameter of all BBH mergers which 
# give SNR in the specified range
def param_in_SNR_range(data, Run, SNRmin, SNRmax):
	i_accpt = np.where((data[:,run_dict[Run]] > SNRmin) & (data[:,run_dict[Run]] < SNRmax))[0]
	param_vals = data[i_accpt, :]
	# for mass ratio
	if mtyp == 'q':
		param_vals = param_vals[:, 0] / param_vals[:, 1]
	# For total/chirp mass
	else: param_vals = param_vals[:, Dict[mtyp]]
	return param_vals, data[i_accpt, :][:, run_dict[Run]]


Md = 100000
# Mass of Sun
Msun = 1.989e30
# Dictionary to identify the input mtyp
Dict = {'Total Mass': 2, 'Chirp Mass': 3}
# SNR Thresholds
Thresh_color = {8: 'b', 13: 'g', 23: 'm'}
# Parameter dict
param = {'Total Mass': {}, 'Chirp Mass': {}, 'q': {}}
# Observed data
Obs = {"GW150914":{"Total Mass": 65., "Chirp Mass": 30., "color": 'm', "SNR": 24}, "LVT151012":{"Total Mass": 36., "Chirp Mass": 15.1, "color": 'c', "SNR": 9.7}, "GW151226":{"Total Mass": 21.7, "Chirp Mass": 8.9, "color": 'g', "SNR": 13}, "GW170104":{"Total Mass": 50.6, "Chirp Mass": 21.1, "color": 'g', "SNR": 13}, "GW170608":{"Total Mass": 19., "Chirp Mass": 7.92, "color": 'g', "SNR": 13}, "GW170814":{"Total Mass": 55.8, "Chirp Mass": 24.16, "color": 'g', "SNR": 18}}
# column no corresponding to run. (S6: 11 .. Design: 14)
run_dict = {'S6': 11, 'O1': 12, 'O2': 13, 'Design': 14}

# Mass limits in Msun
mass_min = 5.
mass_max = 95.

# no of bins
nbins = 30
# Run: 'S6', 'O1', 'O2', 'Design'
RUN = 'O2'
SNR_max = 60.
# Annotation y coordinate placing
annot_y = 0.026     # Has to be adjusted manually
event_col = 'r'

# Parameter histogram
# Type of distribution
for inp_distri in ["Uniform", "Log_Flat", "Power_Law"]:
	# Reading resp data file into an array
	dat = np.genfromtxt('./../Data/Data-for-%s-distri_%s_%s_%s.txt' %(inp_distri, int(mass_min), int(mass_max), Md), delimiter=',')
	# Plot against Total Mass / Chirp Mass / q
	for mtyp in ["Total Mass", "Chirp Mass"]:
		plt.figure(figsize=(19, 10))
		for thresh in Thresh_color:
			param[mtyp][thresh], _ = param_in_SNR_range(dat, RUN, thresh, SNR_max)
			# histogram of the parameter
			plt.hist(param[mtyp][thresh] , nbins, color=Thresh_color[thresh], histtype='step', density=True, label='{}; SNR > {}'.format(mtyp, thresh))
		# Observed events so far
		for event in Obs.keys():
			# Vertical line
		    plt.axvline(x=Obs[event][mtyp], ymin=0., ymax = 1., linewidth=2, color=Obs[event]["color"], alpha=0.5)
		    # Annotations
		    plt.annotate(' {} (SNR {}) '.format(event, Obs[event]["SNR"]), (Obs[event][mtyp]+0.4, annot_y), fontsize=12, fontweight='bold', color=Obs[event]["color"], rotation=90)
		plt.xlabel(mtyp)
		plt.ylabel('Probability of Detection')
		plt.title('Probability of Detection Vs. %s (Distribution: %s)' %(mtyp, inp_distri))
		plt.legend(loc='upper right')
		# plt.axes([0., 100., 0., 1.])
		# if mtyp == "Chirp Mass":
		#     plt.xlim([0., 50.])
		#     plt.ylim([0., 0.06])
		# else:
		#     plt.xlim([0., 100.])
		#     plt.ylim([0., 0.02])
		plt.savefig('./../Plots/%s-%s-%s_%s_%s_%s.png' %(RUN, inp_distri, mtyp, int(mass_min), int(mass_max), Md))
		plt.show()
