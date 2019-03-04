from function_SNR import *
import matplotlib.pyplot as plt

# Taking r = 100 Mpc
r = np.array([100.*Mpc])
r = {'S6':r, 'O1':r, 'O2':r, 'Design':r}
# limiting masses
# Change: to cover full SNR window and avoid "leaking", we find space limits such that lightest binary at closest distance will give highest SNR and heaviest binary at farthest distance will give lowest SNR
mass_lim = {'Lower Space Cutoff': 5., 'Upper Space Cutoff': 50.}
# SNR cutoffs
cut = {'Lower Space Cutoff': 60., 'Upper Space Cutoff': 5.}
# dict of SNRs which will contain SNR dicts against cutoffs. These SNR dicts in turn contain SNR arrays corresponding to diff RUNs
SNRs = {'Lower Space Cutoff': {}, 'Upper Space Cutoff': {}}
# limits on space
r_extreme = {'Lower Space Cutoff': {}, 'Upper Space Cutoff': {}}

# x_l = 0.
# x_u = 3000.
# y_u = 300.
# tolrnc = 0.2
# off_x = 50.
# off_y = 4.
# RUN = 'S6'

for cutoff in cut.keys():
	# creating injection arrays with cutoff mass pairs
	M1 = np.array([mass_lim[cutoff]]) * Msun
	M2 = np.array([mass_lim[cutoff]]) * Msun
	# finding SNRs for respective cutoffs for diff RUNs
	SNRs[cutoff]  = find_pycbc_SNR(M1, M2, r, [0.], [0.], [0.])
	# find extremity in space corresponding to the SNR cutoff
	for RUN in SNRs[cutoff].keys():
		r_extreme[cutoff][RUN] = SNRs[cutoff][RUN] * r[RUN] / cut[cutoff] / Mpc

print r_extreme
#########################################
# 				Output:

# With find_pycbc_SNR
# {'Upper Space Cutoff':
# {'Design': array([3279.12309453]),
# 'S6': array([153.40926056]),
# 'O1': array([901.50600027]),
# 'O2': array([1153.03794345])},
# 'Lower Space Cutoff':
# {'Design': array([1587.82130486]),
# 'S6': array([82.77513475]),
# 'O1': array([454.18511569]),
# 'O2': array([576.59669512])}}


# With find_simple_SNR
# {'Upper Space Cutoff': 
# {'Design': array([3452.14263401]),
#  'O1': array([964.90328158]),
#  'O2': array([1224.40665554]),
#  'S6': array([169.03614154])},
# 'Lower Space Cutoff': 
# {'Design': array([1056.05949735]),
#  'O1': array([179.04932493]),
#  'O2': array([267.89876495]),
#  'S6': array([1.56752347])}
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
# plt.savefig('./../Plots/Limits_on_space', format='png')
plt.show()
