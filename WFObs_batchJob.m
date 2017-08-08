clear all; close all; clc;

%% Batch job code
%  Summary:
%     This code performs a batch of simulations. It iterates over all the
%     configuration files defined in /cluster_jobs/queue/*.
%

% Command window reporting settings
scriptOptions.printProgress     = 1;  % Print progress every timestep
scriptOptions.printConvergence  = 0;  % Print convergence parameters every timestep

% Visualization settings
scriptOptions.plotMesh          = 0;  % Show meshing and turbine locations
scriptOptions.Animate           = 0;  % Show results every x iterations (0: no plots)
   scriptOptions.plotContour    = 0;  % Show flow fields
   scriptOptions.plotPower      = 0;  % Plot true and predicted power capture vs. time
   scriptOptions.plotError      = 0;  % plot RMS and maximum error vs. time
   scriptOptions.plotCenterline = 0;  % Plot centerline speed of the wake (m/s)

% Saving settings
scriptOptions.savePlots         = 0;  % Save all plots in external files at each time step
scriptOptions.saveEst           = 1;  % Save estimated flow fields & powers in an external file at each time step
scriptOptions.saveWorkspace     = 1;  % Save complete workspace at the end of simulation


%% Execute the WFObs core with all configuration files
run('WFObs_addpaths.m');       % Import libraries for WFObs & WFSim
rmpath configurations;         % Remove configurations dir from paths
addpath('cluster_jobs/queue'); % Add cluster queue folder as path

% Collect all configuration files
filenames = dir('cluster_jobs/queue');
filenames = filenames(3:end); % Skip '.' and '..' (navigation folders)

% Perform simulation for each config file
for j = 1:size(filenames,1)
    % Set configuration filename
    filename_tmp = filenames(j).name;
    configName   = filename_tmp;
    
    % Set destination folder for results
    str_i_tmp    = strfind(filename_tmp,'.m');
    filename_tmp = filename_tmp(1:str_i_tmp-1);
    scriptOptions.savePath = ['results/' filename_tmp];
    
    % Perform simulation
    outputData = WFObs_core(scriptOptions,configName);
end