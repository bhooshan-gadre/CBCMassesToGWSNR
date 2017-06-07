# CBCMassesToGWSNR
To study effect of CBC mass distribution on GW signal SNR

## Codes
- **_function_SNR.py_**
	- Module containing definitions of functions to calculate SNRs (simple-dependent only on component masses and distance from source; detailed- containing all the info like sky location and orintation etc.)
	- Reads the ASD data and Antenna sensitivity to sky location data
	- Prforms 1D and 2D interpolation functions of F+, Fx and ASD
	- Contains all the necessary constants
	
- **_Limits_on_visible_space.py_**
	- We are looking at a SNR window of 5 to 60.
	- Lowest value of SNR will be obtained from a (5 Msun, 5 Msun) binary at farthest distance and highest value of SNR from (50 Msun, 50 Msun) binary at the closest distance.
	- Includes two functions:
		- one to calculate the distance at which specific SNR value is obtained
		- second to plot the graphs of SNR vs Distance with the limiting values
		
- **_Alpha_vs_N_data_gen_modified.py_**
	- 
