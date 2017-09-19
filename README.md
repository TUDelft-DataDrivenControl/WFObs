
 ##   WindFarmObserver (WFObs)
Developed by Doekemeijer et al., Delft University of Technology, 2017
            


## Summary:
WindFarmObserver (WFObs) is a collection of state estimators for the nonlinear dynamic flow model WindFarmSimulator (WFSim). The states are the time-varying two-dimensional flow fields in a wind farm. This state information can be used in real-time for state-feedback control algorithms, see work by Vali et al.

## Quick use:
Make sure you clone the repository recursively, by
	
	git clone https://github.com/TUDelft-DataDrivenControl/WFObs --recursive

Once downloaded, open WFObs_demo.m with any recent version of MATLAB. Follow the instructions therein to perform offline estimation simulations of various wind farm scenarios. If necessary, high-fidelity data will be automatically downloaded within WFSim.
	
## Folder hierarchy:

	/bin/:            contains all the functions and scripts used by WFObs.
	/cluster_jobs/:   contains a folder with a demo code to run WFObs on a Linux cluster
	/configurations/: contains the settings files for various simulation scenarios
	/dev_tools/:      contains scripts to analyse results, import raw LES data, and other tools useful for debugging and development
	/setup_sensors/:  contains scripts that will allow you to set up the location of measurements, and generate the relevant input file
	
	/WFSim/:          contains all the model information, identical to the original repository (https://github.com/TUDelft-DataDrivenControl/WFSim)
	
## Debugging:
For any serious issues, reach out to us on the Github page. 

All credit goes to the Delft University of Technology. WFObs was written by ir. Bart Doekemeijer with support from ir Sjoerd Boersma, and under the supervision of dr.ir. Jan-Willem van Wingerden and prof. Lucy Pao.             
