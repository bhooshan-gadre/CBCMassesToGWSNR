from function_SNR import *
import matplotlib.pyplot as plt

# Defining numpy array of r
r = np.arange(1., 5000., 0.125)
r *= Mpc
r = {'S6':r, 'O1':r, 'O2':r, 'Design':r}
# limiting masses
mass_lim = {'Lower Cutoff': 50., 'Upper Cutoff': 5.}
# SNR cutoffs
cut = {'Lower Cutoff': 60., 'Upper Cutoff': 5.}
# dict of SNRs which will contain SNR dicts against cutoffs. These SNR dicts in turn contain SNR arrays corresponding to diff RUNs
SNRs = {'Lower Cutoff': {}, 'Upper Cutoff': {}}
# limits on space
r_extreme = {'Lower Cutoff': {}, 'Upper Cutoff': {}}

# x_l = 0.
# x_u = 3000.
# y_u = 300.
# tolrnc = 0.2
# off_x = 50.
# off_y = 4.
# RUN = 'S6'

for cutoff in cut.keys():
	# creating injection arrays with cutoff mass pairs
	M1 = np.ones(len(r['S6'])) * mass_lim[cutoff] * Msun
	M2 = np.ones(len(r['S6'])) * mass_lim[cutoff] * Msun
	# finding SNRs for respective cutoffs for diff RUNs
	SNRs[cutoff]  = find_simple_SNR(M1, M2, r)
	# find extremity in space corresponding to the SNR cutoff
	for RUN in SNRs[cutoff].keys():
		r_extreme[cutoff][RUN] = r[RUN][np.where(min((SNRs[cutoff][RUN] - cut[cutoff])**2.) == ((SNRs[cutoff][RUN] - cut[cutoff])**2.))]

print r_extreme
#########################################
# 				Output:

# {'Lower Cutoff': {'Design': array([3.25847834e+25]), 'S6': array([5.01423041e+22]), 'O1': array([5.52336765e+24]), 'O2': array([8.26576597e+24])}, 'Upper Cutoff': {'Design': array([1.06521539e+26]), 'S6': array([5.21479962e+24]), 'O1': array([2.97729573e+25]), 'O2': array([3.77802976e+25])}}

# OR

# {'Lower Cutoff': 
# {'Design': array([ 1056.]),
#   'O1': array([ 179.]),
#   'O2': array([ 267.875]),
#   'S6': array([ 1.625])},

# 'Upper Cutoff': 
# {'Design': array([ 3452.125]),
#   'O1': array([ 964.875]),
#   'O2': array([ 1224.375]),
#   'S6': array([ 169.])}}
#########################################

# Plotting
fig = plt.figure(figsize=(15,10))
# for all RUNs
for RUN, colr in zip(r_extreme["Lower Cutoff"].keys(), ["g", "r", "b", "m"]):
	for cutoff in cut.keys():
		# SNR vs r curves
		plt.plot(r[RUN]/Mpc, SNRs[cutoff][RUN], colr+"-", label=RUN + ": " + cutoff)
	plt.fill_between(r[RUN]/Mpc, 0, 100., where=(r[RUN]/Mpc<r_extreme["Upper Cutoff"][RUN]/Mpc) & (r[RUN]/Mpc>r_extreme["Lower Cutoff"][RUN]/Mpc), edgecolor=colr, facecolor=colr, alpha=0.7, label=RUN)
plt.xlim(0., 4000.)
plt.ylim(0., 80.)
# plt.xscale("log")
# plt.yscale("log")
plt.grid(True)
plt.xlabel('Distance (Mpc)')
plt.ylabel('SNR')
plt.title("Limits on observable space for different observational runs of LIGO")
plt.legend(loc="upper right")
# plt.savefig('Limits_on_space', format='png')
plt.show()
