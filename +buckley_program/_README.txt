README for Hogan Lab 2D Inversion Code
(last updated 9/7/2017)

-------------------------------------------------------------------------
Files necessary (in same directory as Main.m):

APMTransfer.m
Cc.m
DMATransfer.m
eps1.m
f_charge.m
f_charge_W.m
InputData.csv
InversionParameters.m
Main.m
readAndPlotResults.m
simpleSmooth.m
smooth.m
Twomey.m
zp2dpCc.m
Results (dummy folder)

IMPORTANT: These scripts are written to operate with MATLAB 2016b or 2017a. Using a previous version will most likely result in bugs, and newer versions have not been tested.

-------------------------------------------------------------------------
Actions necessary:

APMTransfer.m (APM or 2nd dimension of separation transfer function):
	- Confirm APM settings, including: outer and inner radii (R2, R1), APM Length (L), rotational speed (RPM), and aerosol flowrate (Qae). All units in SI.	
	

DMATransfer.m (DMA or 1st dimension of separation transfer function):
	- Confirm DMA settings, including: sheath flowrate (Qsh), aerosol flowrate (Qae), DMA length (L), inner and outer radii (R1 and R2), and gas parameters (T,e)
	- If accounting for cpc counting efficiency and transmission losses, confirm tubing length (Ltube) and uncomment line 65 (G = DMATfer.*Pen.*CPCEff;) and comment line 67 (G = DMATfer;).
	- If inlet and outlet sheath flow are not equal and/or inlet and outlet aerosol flow are not equal, please reference Stolzenburg and McMurry for slight modifications to the transfer function calculation.
	
	
InputData.csv:
	Input measurement data into this file in the format:
		Column 1/A: singly charged diameters classified by DMA (in meters). zp2dpCc.m can be used to convert mobilities (which can be simply calculated from the voltage applied to the DMA) to mobility diameters.
		Column 2/B: singly charged masses classified by APM (in kilograms)
		Columns 3/C on: measurement data. Number of columns = number of APM measurements 
	
	
InversionParameters.m
	- Confirm settings to be used during inversion, including: number of points describing diameters and masses (ll and mm), factor to be used in Twomey routine except in special cases (factor), factor to be used for first few iterations of Twomey loop (factorInit), number of Twomey steps employing factorInit (numInitTwomeySteps), weighting factor for smoothing (smoothFactor), and the maximum number of Twomey steps to be taken in one Twomey loop (maxTwomeySteps).

	
Main.m
	- Confirm aerosol flowrate (Qae), measurement time for each channel (t_ms), CPC error value (errorCPC), and the name which should be used to describe the results when saving (resultsFile)
	- Consider and confirm the values for source error (errorSource) and the buffer error (errorBuffer). Descriptions for these terms and how to estimate them are given in the comments of Main (lines 54-75)
	- Run or submit this file when all of the above have been completed.
	- comment out line 52 ('errorCounting = 0;') to account for measurement error
	
	
-------------------------------------------------------------------------

Tips:

Achieving convergence and modifying InversionParameters.m:
	- If the chi-squared error is increasing significantly with each Twomey step, the most likely culprits are either poor data or incorrect instrument settings (e.g. incorrect length or radii of the APM).
	- If data are good and instrument settings are correct, but chi-squared is slowly decreasing because it is jumping up and down, changing factor or factorInit can be effective in speeding convergence.
	- If chi-squared < 1 is difficult to achieve, but chi-squared converges to a slightly larger value (1-10), it is likely that unknown or unaccounted for error is impacting the inversion. It is reasonable (but should be noted) to increase the overall error by either adding a term to the overall error or increasing the errorCPC value in this case.
	- If ll and mm are changed, smoothFactor may also need to be changed to get same smoothing effect, since smoothing is done on a nearest neighbor basis.
		 
	
Using instruments other than the DMA or APM:
	- generate *.m file similar to DMATransfer.m or APMTransfer.m to represent the transfer function of correct instrument. For ease of use, this file should take similar inputs as DMATransfer.m and APMTransfer.m
	- use find and replace function (MATLAB) in Main.m file to change DMATransfer(...) or APMTransfer(...) to the name of your new *.m file.
	
	
	
	
	
