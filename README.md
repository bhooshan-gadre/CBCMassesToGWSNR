# CBCMassesToGWSNR
To study effect of CBC mass distribution on GW signal SNR


## Codes

1. **_function_SNR.py_**
	- Module containing two functions to calculate SNRs (SIMPLE: dependent only on component masses and distance from source; DETAILED: containing all the parameters like sky location and orientation etc.)
	- Reads the ASD data and Antenna sensitivity to sky location data
	- Prforms 1D and 2D interpolation functions of F+, Fx and ASD
	- Contains all the necessary constants
	
2. **_Limits_on_visible_space.py_**
	- We are looking at a SNR window of 5 to 60.
	- Lowest value of SNR will be obtained from a (5 Msun, 5 Msun) binary at farthest distance and highest value of SNR from (50 Msun, 50 Msun) binary at the closest distance.
	- Includes two functions:
		- one to calculate the distance at which specific SNR value is obtained
		- second to plot the graphs of SNR vs Distance with the limiting values

3. **_Data_gen-M1_M2_M_chirpM_r_alpha_delta_iota_SNR_**
	- Here we inject 10^6 binaries uniformly distributed in space, having random positions in sky, with random orientations with respect to the line of sight.
	- Data for three distributions, viz., uniform (5Msun to 50Msun), log flat and power law (ALPHA = 2.3), is generated
	- For all the above distributions, m1, m2 >= 5Msun and M_tot <= 100Msun
	- We record the data in a .csv file with format 'BH1 Mass (Msun), BH2 Mass (Msun), Total Mass (Msun), Chirp Mass (Msun), Distance (Mpc), alpha, delta, iota, SNR'
	- (Note the difference between (small) alpha (right inclination angle) and (big) ALPHA (parameter determining weightage: M^-ALPHA))
	- The data stored, are sorted according to chirp mass

4. **_Histograms_of_SNR_and_likelihood_of_detection_for_diff_mass_bins_**
	- We saggregate the above recorded data in different mass bins (total mass of binaries).
	- We plot histograms of SNRs for these bins.
	- Mean, median and variance are calculated for each bin
	- From the above histograms we check for the number of observations above specific SNR values (>8, >13 and >23), i.e., area under the curves ahead of respective SNR values
	- These areas above a specific SNR cutoff, for different mass bins give us the PDFs of masses of detectable sources above that perticular SNR threshold
	- Thus we obtain the PDFs of chirp/total mass of detectable binaries




