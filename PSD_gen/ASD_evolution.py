import numpy as np
import matplotlib.pyplot as plt

# TO read and plot ASD for a specific run
def ASD(Run, colr, Alpha):
    # Read ASD data
    dat = np.genfromtxt('./../Data/asd_%s_L1.txt' %Run, delimiter=',')
    # Plot ASD
    plt.loglog(dat[:, 0], dat[:, 1], colr, alpha=Alpha, label=Run)

plt.figure(figsize=(10, 8))
# reading and plotting ASD data for all the runs along with Design ASD
for run, color, A in zip(['S5', 'S6', 'O1', 'O2', 'Design'], ['g', 'g', 'b', 'b', 'k'], [0.4, 0.8, 0.5, 1.0, 1.0]):
    ASD(run, color, A)

# Detailing of graph
fnt = 17
f_min = 10.
f_max = 2048.
plt.axis([f_min, f_max, 1e-24, 1e-19])
plt.grid('on')
plt.xticks(fontsize=fnt)
plt.yticks(fontsize=fnt)
plt.xlabel('Frequency (Hz)', fontsize=fnt)
plt.ylabel('ASD (strain/rtHz)', fontsize=fnt)
plt.title('Evolution of ASD curve', fontsize=fnt)
legend = plt.legend(title='Runs', fontsize=fnt, loc='upper right')
legend.get_title().set_fontsize(fnt) #legend 'Title' fontsize
plt.setp(plt.gca().get_legend().get_texts(), fontsize='15') #legend 'list' fontsize

plt.show()
