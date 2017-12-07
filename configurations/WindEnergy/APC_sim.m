% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%      WFOBS CONFIGURATION FILE  [TEMPLATE LAST UPDATED: 14 Nov. 2017]  %%
%                                                                         %
%        this file contains the complete set of inputs required for the   %
%        observer to do its work. The following settings are defined:     %
%           - Wind farm simulation case from WFSim (Wp.name)              %
%           - Loading of a consistent random seed for 1:1 comparisons     %
%           - Estimation settings of the freestream conditions            %
%           - Measurement definitions and artificial noise levels         %
%           - Kalman filter settings                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% General settings
% Wind farm case to simulate
Wp.name = 'apc_9turb_alm_turb_sim'; % Name of WFSim meshing (from '/WFSim/bin/core/meshing.m')

% Estimate freestream conditions
strucObs.U_Inf.estimate  = false; % Estimate freestream (inflow) u_Inf and v_Inf
   
% Measurement definitions
strucObs.measPw = true;    % Use power measurements (SCADA) from turbines in estimates
strucObs.measFlow = false; % Use flow measurements (LIDAR) in estimates
    strucObs.sensorsPath = ''; % Flow measurement setup filename (see '/setup_sensors/sensors_layouts')     

%% Kalman filter settings
strucObs.se.enabled = true;   % Estimate the model states (flow fields)?
strucObs.pe.enabled = true;   % Estimate model parameters?
strucObs.filtertype = 'sim';  % Choice of observer

% Load universal settings
run('WE_universal_settings.m');  
