from function_SNR import *
import matplotlib.pyplot as plt

	
# This function finds out the value of distance at which we get a specific value of SNR (5 or 60 in our case)
# Gets input in the form (tolerance in SNR, value of SNR at which r is to be found out, array of r, array of SNR)
def r_for_SNR(tol, val, x, y):
	yy = y[np.where(y>(val-tol))]
	# print np.where(yy<(val+tol))
	return x[np.where(yy<(val+tol))]


# Main function which finds SNR, cutoff on distance and plots the graphs
# Takes input as a matrix containing info: [SNRcutoff, M_BH, low_lim_x_axis, uppr_lim_x_axis, uppr_lim_y_axis, tolerance, offset for annotation in x, in y, annotation]
def calc_plot(Mat):
	cut = Mat[0]
	ML = Mat[1]
	x_lL = Mat[2]
	x_L = Mat[3]
	y_L = Mat[4]
	tolrnc = Mat[5]

	M1 = M2 = np.ones(len(r)) * ML * Msun
	alpha = np.array([Mat[9]])
	delta = np.array([Mat[10]])
	iota = np.array([Mat[11]])
	
	SNR  = find_SNR(M1, M2, r, alpha, delta, iota)

	rL = r_for_SNR(tolrnc, cut, r, SNR)/Mpc
	print '(', rL, ', ', cut, ')'

	# Plotting
	fig = plt.figure()
	sb = fig.add_subplot(111)
	sb2 = fig.add_subplot(111)

	sb.plot(r/Mpc, SNR, '-')
	sb2.plot([rL], [cut], 'o')
	sb2.annotate('  (%1.1f Mpc,  %1.1f)' %(rL, cut), (rL + Mat[6], cut + Mat[7]))
	sb.grid(True)
	plt.xlabel('Distance (Mpc)')
	plt.ylabel('SNR')
	plt.title(Mat[8])
	plt.axis([x_lL, x_L, 0., y_L])


# Defining numpy array of r
r = np.linspace(100, 5000, 1001)
r *= Mpc

Mats = [[60., 50., 0., 3000., 300., 0.7, 50., 4., 'Lower Cutoff', 0., 0., np.pi], [5., 5., 0., 5000., 100., 0.01, 70., 4., 'Higher Cutoff', 0., 0., np.pi]]
		# [[60., 50., 0., 3000., 300., 0.1, 50., 4., 'Lower Cutoff', sky_lookup[4844, 0], sky_lookup[4844, 1], np.pi], [5., 5., 0., 5000., 100., 0.2, 70., 4., 'Higher Cutoff', sky_lookup[17039, 0], sky_lookup[17039, 1], np.pi/2.]]

calc_plot(Mats[0])
calc_plot(Mats[1])

plt.show()

