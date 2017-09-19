clear all; close all; clc;
%
% exportFigsFromWorkspace.m
%
% This code will create a folder inside your 'path' called 'figures', where
% all the contour plot and centerline figures will be saved to, using the
% WFObs_s_animations(..) script. It will plot a figure every [x] instances,
% with [x] defined as the variable 'expFreq'.
%
%
time_vec = [400 800 1200];
plotError   = false; % Export final contour, centerline, power-time and error-time plots
plotPower   = true;
powerFC     = 500; % Plot x-seconds open-loop forecast
plotContour = true;
plotCline   = true;

% Visualization settings
scriptOptions.plotMesh          = 0;  % Show meshing and turbine locations
scriptOptions.Animate           = 50;  % Show results every x iterations (0: no plots)
   scriptOptions.plotContour    = 1;  % Show flow fields
   scriptOptions.plotPower      = 1;  % Plot true and predicted power capture vs. time
    scriptOptions.powerForecast = 0;  % Plot power forecast (0 = disabled, x = number of steps) (only if plotPower = 1)
   scriptOptions.plotError      = 0;  % plot RMS and maximum error vs. time
   scriptOptions.plotCenterline = 1;  % Plot centerline speed of the wake (m/s)
   
% Data sources
data = {};
data{end+1} = struct(...
    'name','sim',...
    'path','../../results/WE2017/axi_2turb_alm_turb_sim/workspace.mat');
data{end+1} = struct(...
    'name','ExKF',...
    'path','../../results/WE2017/axi_2turb_alm_turb_measFlow_exkf_stateEst/workspace.mat');
data{end+1} = struct(...
    'name','EnKF',...
    'path','../../results/WE2017/axi_2turb_alm_turb_measFlow_enkf_stateEst/workspace.mat');
data{end+1} = struct(...
    'name','UKF',...
    'path','../../results/WE2017/axi_2turb_alm_turb_measFlow_ukf_stateEst/workspace.mat');

% Add libraries
addpath('../../bin'); % Add binary files from WFObs for plotting
addpath('../../WFSim/libraries/export_fig'); % Add export_fig library
addpath('../../WFSim/bin/core');
addpath('libraries'); % Add libraries for VAF calculations
addpath('../LES_import/bin'); % Used for 'secs2timestr.m'

for di = 1:length(data)
    WS_tmp = load(data{di}.path);
    Wp     = WS_tmp.Wp; 
    sys    = WS_tmp.sys;
    strucObs = WS_tmp.strucObs;
    scriptOptions = WS_tmp.scriptOptions;
    
    
    % Load measurements from LES simulation (*.mat file)
    Wp_tmp = meshing(WS_tmp.Wp.name,false,false);
    LESData    = load(Wp_tmp.sim.measurementFile); % Load measurements
    LESData.ud = LESData.u + strucObs.noise_obs*randn(size(LESData.u)); % Add noise
    LESData.vd = LESData.v + strucObs.noise_obs*randn(size(LESData.v)); % Add noise
    clear Wp_tmp
    
    % Format figures
    scriptOptions.Animate        = 1;
    scriptOptions.plotContour    = plotContour;
    scriptOptions.plotCenterline = plotCline;
    scriptOptions.plotPower      = 0;
    scriptOptions.powerForecast  = powerFC;
    scriptOptions.plotError      = 0;
    scriptOptions.savePlots      = 1;
    
    % Determine and create output directory
    scriptOptions.savePath = [fileparts(data{di}.path) '/figures'];
    mkdir(scriptOptions.savePath);
    
    % Export figures 
    ticLoop = tic;
    NN       = length(time_vec);
    for k = 1:NN
        t = time_vec(k);
        sol_array = (WS_tmp.sol_array(1:t));
        [ ~,scriptOptions ] = WFObs_s_animations( Wp,sol_array,sys,LESData,scriptOptions,strucObs,{} );
        elapsedTime = toc(ticLoop);
        ETA = ceil((NN-k)*(elapsedTime/k));
        disp([data{di}.name ': k = ' num2str(k) '/' num2str(NN) '.  ETA: ' secs2timestr(ETA) '.']);
    end
    if plotPower + plotError > 0
        scriptOptions.plotPower = plotPower;
        scriptOptions.plotError = plotError;    
        [ ~,scriptOptions ] = WFObs_s_animations( Wp,sol_array,sys,LESData,scriptOptions,strucObs,{} );
    end
end