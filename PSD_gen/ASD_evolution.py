import numpy as np
import matplotlib.pyplot as plt

# TO read and plot ASD for a specific run
def ASD(Run, colr, Alpha):
    # Read ASD data
    dat = np.genfromtxt('/home/shreejit/asd_%s_L1.txt' %Run, delimiter=',')
    # Plot ASD
    plt.loglog(dat[:, 0], dat[:, 1], colr, alpha=Alpha, label=Run)

plt.figure(figsize=(10, 8))
# reading and plotting ASD data for all the runs along with Design ASD
for run, color, A in zip(['S5', 'S6', 'O1', 'O2', 'Design'], ['g', 'g', 'b', 'b', 'k'], [0.4, 0.8, 0.5, 1.0, 1.0]):
    ASD(run, color, A)

# Detailing of graph
plt.grid('on')
plt.xlabel('Frequency (Hz)')
plt.ylabel('ASD (strain/rtHz)')
plt.title('Evolution of ASD curves')
plt.legend(title='Runs', loc='upper right')

plt.show()
