%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      WindFarmObserver (WFObs)                        
%                                                                                                   
%  Description:  
%   WindFarmObserver (WFObs) is a state observer that gives accurate
%   estimations of the time-varying two-dimensional flow field    
%   in a wind farm while requiring few measurements (<= 1% of the           
%   variables of interest) and a low computational cost (<= 1 s). This      
%   state information can be used in real-time for state-feedback con-
%   trol algorithms to determine a temporally optimal control policy,
%   accounting for time-varying atmospheric conditions and unmod-
%   eled (and stochastic) flow dynamics, both typically neglected in
%   existing wind farm control methods.
%                                                       
%  Publications on WFObs (sorted on from newest to oldest):
%   1. Doekemeijer, B.M.; Boersma, S.; Pao, L.Y.; van Wingerden, J.W. 'Ensemble Kalman filtering for wind field estimation in wind farms', in proceedings of the American Control Conference (ACC), 2017 
%   2. Doekemeijer, B.M.; van Wingerden, J.W.; Boersma, S.; Pao, L.Y. 'Enhanced Kalman filtering for a 2D CFD NS wind farm flow model', Journal of Physics: Conference Series, Volume 753, Issue 5, 2016, Pages 052015.
%  
%  Recommended usage:
%   1. Make sure your WFSim and WFObs versions are up-to-date
%   2. Specify your script     preferences in lines 49 - 59
%   3. Specify your simulation settings file in line  62
%     a. To create a simulation case manually, have a look at the ./configurations folder
%        i.  Use  '/Setup_sensors/WFObs_s_sensors.m' to specify measurement locations 
%        ii. All WFSim model/site information is contained in '../WFSim/bin/core/meshing.m'
%
%  Filter options:
%   1. Ensemble Kalman Filter (EnKF) including localization, inflation,
%      and parallelization. Best results in both accuracy and speed.
%   2. Extended Kalman Filter (ExKF) including sparsification options
%      to reduce computational cost. Good reconstruction results, but slow.
%   3. Do Nothing (sim), which simply simulates the WFSim model without
%      using any kind of filtering/reconstruction algorithm.
%
%
%  Last updated: 04/07/2017
%  Author:       Bart Doekemeijer, TUDelft (B.M.Doekemeijer@tudelft.nl)
%                                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Import libraries
clear all; close all; clc; tic; 
addpath bin                           % Observer directory
addpath WFSim/bin/core                % WFSim model directory
addpath WFSim/libraries/sparse_null   % Model supplementary library
addpath WFSim/libraries/export_fig    % Graphics library (get here: http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig)

%% Script settings
strucScript = struct(...
    'Animation'       , 10, ...  % Plot figures every # iteration (no animation: 0)
       'plotcontour'  , 1, ...  % plot flow fields (contourf)
       'plotpower'    , 0, ...  % Plot true and predicted power capture vs. time
       'ploterror'    , 0, ...  % plot RMS and maximum error vs. time
    'plotMesh'        , 0, ...  % Plot meshing layout (grid)
    'saveplots'       , 1, ...  % Save all plots in external files at each time step
    'saveest'         , 1, ...  % Save estimated flow fields & powers in an external file at each time step
    'saveworkspace'   , 0, ...  % Save complete workspace at the end of simulation
    'savepath'        , ['Results/yawCase3_UKF/'] ... % Destination folder of saved files
    );  

%% Model and observer configuration file
configName = 'YawCase3.m'; % configuration filename. See './configurations' for options: 'NoPrecursor', 'YawCase3'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%                    Internal code                      %% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WFObs_s_initialize;     % Initialize model, observer, and script
for k = [1 (2+strucObs.obsv_delay):1:Wp.sim.NN ]; 
    tic;                          % Start stopwatch for iteration k
    it        = 0;                % Convergence parameter
    eps       = 1e19;             % Convergence parameter
    epss      = 1e20;             % Convergence parameter
    sol.k     = k;                % Write timestep to solution
    timeindex = Wp.sim.time(k+1); % Setup index number corresponding to time at k
    
    % Load measurement data
    measured = WFObs_s_loadmeasurements( sourcepath, datanroffset, timeindex, strucObs.noise_obs ); 
     
    while ( eps>conv_eps && it<max_it && eps<epss ); % Convergence to a solution
        it   = it+1; epss = eps;        
        if k>1; max_it = max_it_dyn; end;
        % Pre-processing: update freestream conditions from SCADA data
        [sol,Wp,B1,B2,bc,strucObs] = WFObs_s_determineFreestream(Wp,input{k},measured,sol,strucObs);
        
        % Calculate optimal solution according to filter of choice
        [sol,strucObs]      = WFObs_o(strucObs,Wp,sys,B1,B2,bc,input{timeindex},measured,sol,k,it,options);  % Perform the observer update
        [sol,eps]           = MapSolution(Wp.mesh.Nx,Wp.mesh.Ny,sol,k,it,options);   % Map solution to flowfields
        [~,~,~,sol.power,~] = Actuator(Wp,input{timeindex},sol,options);             % Calculate power after analysis update
        
    end;        
    WFObs_s_reporting;  % Display progress in the command window
    WFObs_s_animations; % Display animations on screen
end;

if strucScript.saveworkspace; save([strucScript.savepath 'workspace.mat']); end; % save workspace to folder
disp([datestr(rem(now,1)) ' __  Completed simulations']);