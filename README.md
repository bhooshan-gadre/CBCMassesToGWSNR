# CBCMassesToGWSNR
To study effect of CBC mass distribution on GW signal SNR


## Codes

1. **_function_SNR.py_**
	- Module containing definitions of functions to calculate SNRs (simple-dependent only on component masses and distance from source; detailed- containing all the info like sky location and orintation etc.)
	- Reads the ASD data and Antenna sensitivity to sky location data
	- Prforms 1D and 2D interpolation functions of F+, Fx and ASD
	- Contains all the necessary constants
	
2. **_Limits_on_visible_space.py_**
	- We are looking at a SNR window of 5 to 60.
	- Lowest value of SNR will be obtained from a (5 Msun, 5 Msun) binary at farthest distance and highest value of SNR from (50 Msun, 50 Msun) binary at the closest distance.
	- Includes two functions:
		- one to calculate the distance at which specific SNR value is obtained
		- second to plot the graphs of SNR vs Distance with the limiting values
		
3. **_Data_gen-Gamma_N.py_**
	- We take different values of merger rates and calculate the no of observations (> SNR13) taking different values of gamma
	- Here gamma is the combined parameter obtained as follows.
	- We assume mass dependence of source distribution coming from merger time of binaries (M^alpha) and number density of binaries (M^-beta). Thus, combined mass dependence becomes M^(alpha-beta) => M^gamma (say).
	- Also, we create a separate file for the N-gamma cobinations giving 6 detections per year (2 detections in 4 months)
	
4. **_Plot_gamma_vs_N_(2D_histo).py_**
	- 2D histogram with N and gamma as its axes, shows colour coded density of bins of N-gamma combinations giving 6 detections per year.

5. **_Data_gen-M1_M2_r_alpha_delta_iota_SNR_**
	- Here we inject 10^6 binaries uniformly distributed in space, having random positions in sky, with random orientations wrt line of sight.
	- We record the data in a .csv file with format 'BH1 Mass (Msun), BH2 Mass (Msun), Total Mass (Msun), Distance (Mpc), alpha, delta, iota, SNR'
	
6. **_Histograms_of_SNR_for_diff_mass_bins_**
	- We saggregate the above recorded data in different mass bins (total mass of binaries).
	- We weigh the data with weighing factor M^gamma
	- We plot histograms of SNRs for these bins.
	- Mean, median and variance are calculated for each bin
	
7. **_Fraction_of_mergers_above_SNR_8_13_**
	- From the above histograms we check for the number of observations above specific SNR value (8 and 13 here)
	- These numbers, plotted against the mean mass of the bins, give the likelihood of the mass bins to give a detection with SNR above a particular value.
	- Thus we obtain the relation between total mass of binaries and their likelihood to be detected.




























	
