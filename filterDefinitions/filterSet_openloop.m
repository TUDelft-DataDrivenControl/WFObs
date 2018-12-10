function [strucObs] = filterSet_TORQUE_axi_2turb_alm_turb()
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
% Load a predefined random seed (for one-to-one comparisons between simulation cases)
strucObs.loadRandomSeed = true;

% Estimate freestream conditions
strucObs.U_Inf.estimate  = false;  % Estimate freestream (inflow) u_Inf and v_Inf
strucObs.U_Inf.intFactor = 0.99;  % LPF gain (1: do not change, 0: instant change)
strucObs.U_Inf.pwrScaling = 0.95; % pwrScaling = P_meas/P_ADM, roughly. Scaling factor for difference between ADM and measurement data, including efficiency losses etc.

%% Kalman filter settings
% State estimation settings
strucObs.se.enabled = false; % Estimate the model states (flow fields)?
% Initial state error covariance matrix (diagonal values on P_0 matrix)
strucObs.se.P0.u = 0.10; % Initial state covariance for long. flow
strucObs.se.P0.v = 0.10; % Initial state covariance for lat. flow

% Process noise covariance matrix (diagonal values on Q matrix)
strucObs.se.Qk.u = 1e-2; % Autocovariance for long. flow process noise
strucObs.se.Qk.v = 1e-4; % Autocovariance for lat.  flow process noise

% Measurement noise covariance matrix (diagonal values on R matrix)
strucObs.se.Rk.P = 4e8;  % Autocovariance for turbine power measurements [if strucObs.measPw == 1]
strucObs.se.Rk.u = 1e-2; % Autocovariance for long. flow measurements  [if strucObs.measFlow == 1]
strucObs.se.Rk.v = 1e-2; % Autocovariance for lat.  flow measurements  [if strucObs.measFlow == 1]

% Covariance entries for pressure (only applicable if modelOptions.exportPressures = true)
strucObs.se.P0.p = 0.00; % Initial state covariance for pressure terms
strucObs.se.Qk.p = 0.00; % Autocovariance for pressure process noise

% Parameter estimation settings
strucObs.pe.enabled = false; % Estimate model parameters?
strucObs.pe.vars = {'site.lm_slope'}; % If empty {} then no estimation
strucObs.pe.P0   = [0.5]; % Initial state covariance(s) for model variable(s)
strucObs.pe.Qk   = [1e-4]; % Autocovariance(s) process noise for model variable(s)
strucObs.pe.lb   = [0.001]; % Lower bound(s) for model variable(s)
strucObs.pe.ub   = [4.000]; % Upper bound(s) for model variable(s)

% Observer-specific settings
strucObs.filtertype = 'sim';
switch lower(strucObs.filtertype)
    case {'exkf'} % Extended Kalman filter (ExKF)
        % ... The ExKF does not have any specific model settings
        disp('WARNING: You selected to perform an ExKF simulation.')
        disp('Please ensure that modelOptions.Linearversion = true.')
        disp(' ')
        
    case {'ukf'}  % Unscented Kalman filter (UKF)
        % Sigma-point generation settings
        strucObs.alpha = 1.0; % Tuning parameter
        strucObs.beta  = 2.0; % Tuning parameter (2 optimal for Gaussian distributions)
        strucObs.kappa = 0.0; % Tuning parameter (state est.: "0" or param est.: "3-L")
        
    case {'enkf'} % Ensemble Kalman filter (EnKF)
        strucObs.nrens  = 50;        % Ensemble size
        strucObs.r_infl = 1.025;     % Covariance inflation factor (typically 1.00-1.20, no inflation: 1)
        strucObs.f_locl = 'gaspari'; % Localization method: 'off', 'gaspari' (Gaspari-Cohn 1999) or 'heaviside' (Heaviside step function: 0s or 1s)
        strucObs.l_locl = 131;       % Gaspari-Cohn: typically sqrt(10/3)*L with L the cut-off length. Heaviside: cut-off length (m).
        
    case {'sim'} % No filter, just open-loop simulation
        disp('WARNING: You selected to perform an open-loop simulation.')
        disp('Please ensure that modelOptions.exportPressures = true.')
        disp(' ')
        
    otherwise
        error('not a valid filter/simulation specified.');
end
end