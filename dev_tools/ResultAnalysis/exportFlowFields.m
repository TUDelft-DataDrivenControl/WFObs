clear all; close all; clc;
%
% exportFlowFields.m
%
% This code will create a folder inside your 'path' called 'figures', where
% all the contour plot and centerline figures will be saved to, using the
% WFObs_s_animations(..) script. It will plot a figure every [x] instances,
% with [x] defined as the variable 'expFreq'.
%
%
expFreq = 10; % Export figure for every * steps

% Data sources
data = {};
data{end+1} = struct(...
    'name','Sim (nominal)',...
    'path','../../results/YawCase3_sim_nominal/workspace.mat');
data{end+1} = struct(...
    'name','Sim (tuned)',...
    'path','../../results/YawCase3_sim_tuned/workspace.mat');
% Periodic parameter tuning
data{end+1} = struct(...
    'name','UKF (nominal)',...
    'path','../../results/YawCase3_UKF_nominal/workspace.mat');
data{end+1} = struct(...
    'name','UKF (tuned)',...
    'path','../../results/YawCase3_UKF_tuned/workspace.mat');
data{end+1} = struct(...
    'name','EnKF (nominal)',...
    'path','../../results/YawCase3_EnKF_nominal/workspace.mat');
data{end+1} = struct(...
    'name','EnKF (tuned)',...
    'path','../../results/YawCase3_EnKF_tuned/workspace.mat');


% Add libraries
addpath('../../bin'); % Add binary files from WFObs for plotting
addpath('../../WFSim/libraries/export_fig'); % Add export_fig library
addpath('libraries'); % Add libraries for VAF calculations
addpath('../LES_import/bin'); % Used for 'secs2timestr.m'

hFigs = {};
for di = 1:length(data)
    WS_tmp = load(data{di}.path);
    Wp_import  = WS_tmp.Wp;
    Wp         = Wp_import; 
    Wp.turbine = rmfield(Wp.turbine,'input'); 
    for j = 1:1794
        Wp.turbine.input(j) = Wp_import.turbine.input{j};
    end
    strucObs      = WS_tmp.strucObs;
    scriptOptions = WS_tmp.scriptOptions;
    
    % Format figures
    scriptOptions.Animate        = 1;
    scriptOptions.plotContour    = 1;
    scriptOptions.plotCenterline = 1;
    scriptOptions.plotPower      = 0;
    scriptOptions.plotError      = 0;
    scriptOptions.savePlots      = 1;
    
    % Determine and create output directory
    scriptOptions.savePath = [fileparts(data{di}.path) '/figures'];
    mkdir(scriptOptions.savePath);
    
    % Export figures 
    ticLoop = tic;
    time_vec = 1:expFreq:length(WS_tmp.sol_array);
    NN       = length(time_vec);
    for k = 1:NN
        t = time_vec(k);
        sol_array = {WS_tmp.sol_array{1:t}};
        [ hFigs,scriptOptions ] = WFObs_s_animations( Wp,sol_array,scriptOptions,strucObs,hFigs );
        elapsedTime = toc(ticLoop);
        ETA = ceil((NN-k)*(elapsedTime/k));
        disp([data{di}.name ': k = ' num2str(k) '/' num2str(NN) '.  ETA: ' secs2timestr(ETA) '.']);
    end
    scriptOptions.plotPower = 1;
    scriptOptions.plotError = 1;    
    [ hFigs,scriptOptions ] = WFObs_s_animations( Wp,WS_tmp.sol_array,scriptOptions,strucObs,hFigs );
end