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
Obs = {"GW150914":{"RUN": "O1", "Chirp Mass": 28.6, "Total Mass": 63.1, "color": 'm', "SNR": 23.6}, \
	"GW151012":{"RUN": "O1", "Chirp Mass": 15.2, "Total Mass": 35.7, "color": 'm', "SNR": 9.5}, \
	"GW151226":{"RUN": "O1", "Chirp Mass": 8.9, "Total Mass": 20.5, "color": 'm', "SNR": 13.1}, \
	"GW170104":{"RUN": "O2", "Chirp Mass": 21.5, "Total Mass": 49.1, "color": 'm', "SNR": 13}, \
	"GW170608":{"RUN": "O2", "Chirp Mass": 7.9, "Total Mass": 17.8, "color": 'm', "SNR": 15.4}, \
	"GW170729":{"RUN": "O2", "Chirp Mass": 35.7, "Total Mass": 80.3, "color": 'm', "SNR": 9.8}, \
	"GW170809":{"RUN": "O2", "Chirp Mass": 25, "Total Mass": 56.4, "color": 'm', "SNR": 12.2}, \
	"GW170814":{"RUN": "O2", "Chirp Mass": 24.2, "Total Mass": 53.4, "color": 'm', "SNR": 16.3}, \
	"GW170818":{"RUN": "O2", "Chirp Mass": 26.7, "Total Mass": 59.8, "color": 'm', "SNR": 11.3}, \
	"GW170823":{"RUN": "O2", "Chirp Mass": 29.3, "Total Mass": 65.6, "color": 'm', "SNR": 11.1}}
# Color coding according to SNR thresholds
for thr in Thresh_color.keys():
	for obs in Obs.keys():
		if Obs[obs]["SNR"] >= thr: Obs[obs]["color"] = Thresh_color[thr]
# column no corresponding to run. (S6: 11 .. Design: 14)
run_dict = {'S6': 11, 'O1': 12, 'O2': 13, 'Design': 14}

# Mass limits in Msun
mass_min = 5.
mass_max = 95.

# no of bins
nbins = 30
# Runs: 'S6', 'O1', 'O2', 'Design'
SNR_max = 60.
# Annotation y coordinate placing
annot_y = 0.026     # Has to be adjusted manually
event_col = 'r'

# line width
lwid = 3
# fontsize
fsiz = 13

# Parameter histogram
for RUN in ["O1", "O2"]:
	# Type of distribution
	for inp_distri in ["Uniform", "Log_Flat", "Power_Law"]:
		# Reading resp data file into an array
		dat = np.genfromtxt('./../Data/Data-for-%s-distri_%s_%s_%s_3march.txt' %(inp_distri, int(mass_min), int(mass_max), Md), delimiter=',')
		# Plot against Total Mass / Chirp Mass / q
		for mtyp in ["Total Mass", "Chirp Mass"]:
			plt.figure(figsize=(7, 7))
			for thresh in Thresh_color:
				param[mtyp][thresh], _ = param_in_SNR_range(dat, RUN, thresh, SNR_max)
				# histogram of the parameter
				plt.hist(param[mtyp][thresh] , nbins, color=Thresh_color[thresh], linewidth=lwid, histtype='step', density=True, label='SNR > {}'.format(thresh))
			# Observed events so far
			for event in Obs.keys():
				# Selecting only the observations done in resp runs
				if Obs[event]["RUN"] == RUN:
					# Vertical line
				    plt.axvline(x=Obs[event][mtyp], ymin=0., ymax = 1., color=Obs[event]["color"], linewidth=lwid, alpha=0.5)
				    # Annotations
				    # plt.annotate(' {} (SNR {}) '.format(event, Obs[event]["SNR"]), (Obs[event][mtyp]+0.4, annot_y), fontsize=fsiz, fontweight='bold', color=Obs[event]["color"], rotation=90)
				    plt.tick_params(axis='both', which='major', labelsize=fsiz)
			plt.xlabel(mtyp, fontsize=fsiz, fontweight='bold')
			plt.ylabel('Probability of Detection', fontsize=fsiz, fontweight='bold')
			# plt.title('Probability of Detection Vs. %s \n(Distribution: %s)' %(mtyp, inp_distri), fontsize=fsiz, fontweight='bold')
			plt.title('Run: %s, Distribution: %s' %(RUN, inp_distri), fontsize=fsiz, fontweight='bold')
			plt.legend(loc='upper right', fontsize=fsiz)
			# plt.axes([0., 100., 0., 1.])
			# if mtyp == "Chirp Mass":
			#     plt.xlim([0., 50.])
			#     plt.ylim([0., 0.06])
			# else:
			#     plt.xlim([0., 100.])
			#     plt.ylim([0., 0.02])
			plt.savefig('./../Plots/bold_PDF_%s-%s-%s_%s_%s_%s.png' %(RUN, mtyp, inp_distri, int(mass_min), int(mass_max), Md))
			plt.show()













# Cumulative 2D Histogram
# Type of distribution
# for inp_distri in ["Power_Law"]:	#"Uniform", "Log_Flat", "Power_Law"
# 	# Reading resp data file into an array
# 	dat = np.genfromtxt('./../Data/Data-for-%s-distri_%s_%s_%s.txt' %(inp_distri, int(mass_min), int(mass_max), Md), delimiter=',')
# 	# Plot against Total Mass / Chirp Mass / q
# 	for mtyp in ["Total Mass"]:		#, "Chirp Mass"
# 		plt.figure(figsize=(19, 10))
# 		thresh = 5.
# 		param[mtyp][thresh], SNRs = param_in_SNR_range(dat, RUN, thresh, SNR_max)
# 		# 2D histogram of SNRs vs parameter
# 		H, xedges, yedges = np.histogram2d(param[mtyp][thresh], SNRs, bins=[nbins, nbins])
# 		H = H.T  # Let each row list bins with common y range.
# 		H = H[::-1].cumsum(axis=0)[::-1]
# 		b = H.sum(axis=1)
# 		H = H/b[:, np.newaxis]
# 		# Discrete
# 		X, Y = np.meshgrid(xedges, yedges)
# 		plt.pcolormesh(X, Y, H, cmap='RdBu')
# 		plt.colorbar()
# 		# Smooth
# 		# im = mpl.image.NonUniformImage(ax, interpolation='bilinear')
# 		# xcenters = (xedges[:-1] + xedges[1:]) / 2
# 		# ycenters = (yedges[:-1] + yedges[1:]) / 2
# 		# im.set_data(xcenters, ycenters, H)
# 		# ax.images.append(im)
# 		# Observed events so far
# 		for event in Obs.keys():
# 			# Points
# 		    plt.plot(Obs[event][mtyp], Obs[event]["SNR"], marker='*', markersize=10, color=event_col)
# 		    # Annotations
# 		    # plt.annotate(' {} (SNR {}) '.format(event, Obs[event]["SNR"]), (Obs[event][mtyp]+0.2, Obs[event]["SNR"]+3.), fontsize=12, fontweight='bold', color=event_col, rotation=90)
# 		plt.xlabel(mtyp)
# 		plt.ylabel('SNR Threshold')
# 		plt.title('Probability of Detection : SNR Threshold Vs. %s (Distribution: %s)' %(mtyp, inp_distri))
# 		plt.legend(loc='upper right')
# 		# plt.savefig('./../Plots/O2-%s-%s-%s_%s_%s_%s.png' %('2D-cum', inp_distri, mtyp, int(mass_min), int(mass_max), Md))
# 		plt.show()




# H, xedges, yedges = np.histogram2d(M, SNRs, bins=[30, 30])
# H = H.T  # Let each row list bins with common y range.
# H = H.cumsum(axis=0)[::-1]

# # :func:`imshow <matplotlib.pyplot.imshow>` can only display square bins:

# fig = plt.figure(figsize=(7, 3))
# ax = fig.add_subplot(131, title='imshow: square bins')
# plt.imshow(H, interpolation='nearest', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])

# # :func:`pcolormesh <matplotlib.pyplot.pcolormesh>` can display actual edges:

# ax = fig.add_subplot(132, title='pcolormesh: actual edges', aspect='equal')
# X, Y = np.meshgrid(xedges, yedges)
# ax.pcolormesh(X, Y, H)

# # :class:`NonUniformImage <matplotlib.image.NonUniformImage>` can be used to
# # display actual bin edges with interpolation:

# ax = fig.add_subplot(133, title='NonUniformImage: interpolated', aspect='equal', xlim=xedges[[0, -1]], ylim=yedges[[0, -1]])
# im = mpl.image.NonUniformImage(ax, interpolation='bilinear')
# xcenters = (xedges[:-1] + xedges[1:]) / 2
# ycenters = (yedges[:-1] + yedges[1:]) / 2
# im.set_data(xcenters, ycenters, H)
# ax.images.append(im)
# plt.show()
