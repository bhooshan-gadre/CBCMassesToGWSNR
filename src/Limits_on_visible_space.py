from function_SNR import *
import matplotlib.pyplot as plt


# This function finds out the value of distance at which we get a specific value of SNR (5 or 60 in our case)
# Gets input in the form (tolerance in SNR, value of SNR at which r is to be found out, array of r, array of SNR)
def r_for_SNR(tol, val, x, y):
	yy = y[np.where(y>(val-tol))]
	# print np.where(yy<(val+tol))
	return x[np.where(yy<(val+tol))]
	
# Main function which finds SNR, cutoff on distance and plots the graphs
# Takes input : (SNRcutoff, M_BH, low_lim_x_axis, uppr_lim_x_axis, uppr_lim_y_axis, tolerance, offset for annotation in x, in y, annotation)
def calc_plot(cut, ML, x_l, x_u, y_u, tolrnc, off_x, off_y, annot):
	M1 = np.ones(len(r)) * ML * Msun
	M2 = np.ones(len(r)) * ML * Msun
	SNR  = find_simple_SNR(M1, M2, r)
	rL = r_for_SNR(tolrnc, cut, r, SNR)/Mpc
	print '(', rL, ', ', cut, ')'

	# Plotting
	fig = plt.figure()
	# For the plot SNR vs. r
	sb = fig.add_subplot(111)
	# For the limiting point
	sb2 = fig.add_subplot(111)
	# Plot SNR vs. r
	sb.plot(r/Mpc, SNR, '-')
	# Plot the limiting point
	sb2.plot([rL], [cut], 'o')
	sb2.annotate('  (%1.1f Mpc,  %1.1f)' %(rL, cut), (rL + off_x, cut + off_y))
	sb.grid(True)
	plt.xlabel('Distance (Mpc)')
	plt.ylabel('SNR')
	plt.title(annot)
	plt.axis([x_l, x_u, 0., y_u])

# Defining numpy array of r
r = np.linspace(100, 5000, 1001)
r *= Mpc

# Input : (SNRcutoff, M_BH, low_lim_x_axis, uppr_lim_x_axis, uppr_lim_y_axis, tolerance, offset for annotation in x, in y, annotation)
calc_plot(60., 50., 0., 3000., 300., 0.2, 50., 4., 'Lower Cutoff')
calc_plot(5., 5., 0., 5000., 100., 0.005, -700., 4., 'Higher Cutoff')

plt.show()
