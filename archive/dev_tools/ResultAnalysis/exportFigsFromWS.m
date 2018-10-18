clear all; close all; clc;
%
% exportFigsFromWorkspace.m
%
time_vec = [1:50:3500];

saveFigs = true;
savePath = ['figures/6turb_varyingTurb/sim/']; % Save path

% Visualization settings
expOptions.plotContour    = 1;  % Show flow fields
expOptions.plotPower      = 0;  % Plot true and predicted power capture vs. time
 expOptions.powerForecast = 0;  % Plot power forecast (0 = disabled, x = number of steps) (only if plotPower = 1)
expOptions.plotError      = 0;  % plot RMS and maximum error vs. time
expOptions.plotCenterline = 0;  % Plot centerline speed of the wake (m/s)
   
% Data sources
data = {};
data{end+1} = struct(...
    'name','Sim',...
    'path','D:\bmdoekemeijer\My Documents\MATLAB\WFObs\results\TORQUE2018\2turb_varyingTurb_sim\workspace.mat');

% data{end+1} = struct(...
%     'name','EnKF (old tuning)',...
%     'path','D:\bmdoekemeijer\My Documents\SurfDrive\PhD\Dissemination\2017 Wind Energy\MATLAB_out\APC\apc_9turb_alm_turb_n50\workspace.mat');
% data{end+1} = struct(...
%     'name','EnKF (new Q, n=50)',...
%     'path','D:\bmdoekemeijer\My Documents\SurfDrive\PhD\Dissemination\2017 Wind Energy\MATLAB_out\APC\APC_enkf_newQ_n50\workspace.mat');
% data{end+1} = struct(...
%     'name','EnKF (new Q, n=70)',...
%     'path','D:\bmdoekemeijer\My Documents\SurfDrive\PhD\Dissemination\2017 Wind Energy\MATLAB_out\APC\APC_enkf_newQ_n70\workspace.mat');
% data{end+1} = struct(...
%     'name','EnKF (Power meas.)',...
%     'path','D:\bmdoekemeijer\My Documents\SurfDrive\PhD\Dissemination\2017 Wind Energy\MATLAB_out\MEAScomparison_2turb_alm_turb_enkf_n30_Pwr\workspace.mat');
% data{end+1} = struct(...
%     'name','EnKF (Upw. LiDAR)',...
%     'path','D:\bmdoekemeijer\My Documents\SurfDrive\PhD\Dissemination\2017 Wind Energy\MATLAB_out\MEAScomparison_2turb_alm_turb_enkf_n30_upwLidar\workspace.mat');

% data{end+1} = struct(...
%     'name','ExKF',...
%     'path','../../results/WE2017/axi_2turb_alm_turb_measFlow_exkf_stateEst/workspace.mat');
% data{end+1} = struct(...
%     'name','EnKF',...
%     'path','../../results/WE2017/axi_2turb_alm_turb_measFlow_enkf_stateEst/workspace.mat');
% data{end+1} = struct(...
%     'name','UKF',...
%     'path','../../results/WE2017/axi_2turb_alm_turb_measFlow_ukf_stateEst/workspace.mat');

% Add libraries
addpath('../../bin'); % Add binary files from WFObs for plotting
addpath('../../WFSim/libraries/export_fig'); % Add export_fig library
addpath('../../WFSim/bin/core');
addpath('libraries'); % Add libraries for VAF calculations
addpath('../LES_import/bin'); % Used for 'secs2timestr.m'

for di = 1:length(data)
    WS_tmp     = load(data{di}.path);
    Wp         = WS_tmp.Wp; 
    sys        = WS_tmp.sys;
    options    = WS_tmp.scriptOptions;
    strucObs   = WS_tmp.strucObs;
    sol_array  = WS_tmp.sol_array;
    clear WS_tmp
    
    % Load measurements from LES simulation (*.mat file)
    Wp_tmp = meshing(Wp.name,false,false);
    LESData = load(Wp_tmp.sim.measurementFile); % Load measurements
    clear Wp_tmp
    
    % Format figures
    options           = expOptions;
    options.Animate   = 1;  % Show results every x iterations (0: no plots)
    options.savePlots = saveFigs;
    options.savePath  = [savePath '/' data{di}.name];
    
    % Determine and create output directory
    mkdir(options.savePath);
    
    % Export figures 
    ticLoop = tic;
    NN       = length(time_vec);
    
    for k = 1:NN
        if k == 1 || NN < 6; hFigs = {}; end
        t = time_vec(k);
        sol_array_short = (sol_array(1:t));
        [ hFigs,options ] = WFObs_s_animations( Wp,sol_array_short,sys,LESData,options,strucObs,hFigs );
        elapsedTime = toc(ticLoop);
        ETA = ceil((NN-k)*(elapsedTime/k));
        disp([data{di}.name ': k = ' num2str(k) '/' num2str(NN) '.  ETA: ' secs2timestr(ETA) '.']);
    end
%     if plotPower + plotError > 0
%         options.plotPower = plotPower;
%         options.plotError = plotError;    
%         [ ~,options ] = WFObs_s_animations( Wp,sol_array,sys,LESData,options,strucObs,{} );
%     end
end