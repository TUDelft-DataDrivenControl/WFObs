clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%             WIND FARM OBSERVER (WFObs) by B.M. Doekemeijer
%                 Delft University of Technology, 2017
%              Repo: https://github.com/Bartdoekemeijer/WFObs
%
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%
%%   Quick use:
%     1. Specify your script preferences in lines 49 - 59
%     3. Specify your simulation settings file in line  62
%      a. To create a simulation case manually, have a look at the ./configurations folder
%          i.  Use  '/Setup_sensors/WFObs_s_sensors.m' to specify measurement locations
%          ii. All WFSim model/site information is contained in '../WFSim/bin/core/meshing.m'
%     3. Press start.
%
%%   Relevant input/output variables

%%   Debugging and contributing:
%     - First, try to locate any errors by turning all possible outputs
%       on (printProgress, printConvergence, Animate, plotMesh).
%     - If you cannot solve your problems, reach out on the Github.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define script settings
% Command window reporting settings
scriptOptions.printProgress     = 1;  % Print progress every timestep
scriptOptions.printConvergence  = 0;  % Print convergence parameters every timestep

% Visualization settings
scriptOptions.plotMesh          = 0;  % Show meshing and turbine locations
scriptOptions.Animate           = 0;  % Show results every x iterations (0: no plots)
   scriptOptions.plotContour    = 1;  % Show flow fields
   scriptOptions.plotPower      = 0;  % Plot true and predicted power capture vs. time
   scriptOptions.plotError      = 0;  % plot RMS and maximum error vs. time
   scriptOptions.plotCenterline = 1;  % Plot centerline speed of the wake (m/s)

% Saving settings
scriptOptions.savePlots         = 0;  % Save all plots in external files at each time step
scriptOptions.saveEst           = 1;  % Save estimated flow fields & powers in an external file at each time step
scriptOptions.saveWorkspace     = 1;  % Save complete workspace at the end of simulation
scriptOptions.savePath          = ['results/YawCase3_enkf']; % Destination folder of saved files

% Model and observer configuration file
configName = 'YawCase3.m'; % configuration filename. See './configurations' for options: 'NoPrecursor', 'YawCase3'



%% Execute the WFObs core code
run('WFObs_addpaths.m'); % Import libraries for WFObs & WFSim
outputData = WFObs_core(scriptOptions,configName);