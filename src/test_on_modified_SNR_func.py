from function_SNR import *

# Input
rl = 648.8 * Mpc
ru = 2030.6 * Mpc
Md = 10
Iter = 100

vol = 4. / 3. * np.pi * (ru ** 3. - rl ** 3.)

# Arrays
r = np.linspace(0., vol, Md)
r *= 3. / (4. * np.pi)
r += rl ** 3.
r **= (1. / 3.)

M1 = np.random.uniform(5. * Msun, 50. * Msun, Md)
M2 = np.random.uniform(5. * Msun, 50. * Msun, Md)
alpha = np.random.uniform(0., 2. * np.pi, Md)
delta = np.random.uniform(-np.pi/2., np.pi/2., Md)
iota = np.random.uniform(0., np.pi, Md)

SNR = find_SNR(M1, M2, r, alpha, delta, iota)

print SNR
