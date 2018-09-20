%% SOWFA source directories, meshing and measurement options
Wp.name = '6turb_adm_turb';  % Name of meshing (from '/WFSim/bin/core/meshing.m')

%% WFSim model settings
scriptOptions.startUniform = false; % Start from a uniform flow field (1) or from a fully developed waked flow field (0).
scriptOptions.conv_eps     = 1e-6;  % Convergence parameter
scriptOptions.max_it_dyn   = 1;     % Convergence parameter

if scriptOptions.startUniform
    scriptOptions.max_it = 1;   % Iteration limit for simulation start-up
else
    scriptOptions.max_it = 50;  % Iteration limit for simulation start-up
end


%% Filter settings
% General settings
strucObs.loadRandomSeed = true; % Load a predefined random seed (for one-to-one comparisons between simulation cases)
strucObs.noise_obs      = 0.1;  % Disturbance amplitude (m/s) in output data by randn*noiseampl ('0' for no noise)
strucObs.noise_init     = 0.0;  % Disturbance amplitude (m/s) in initial flow field by randn*noiseinit ('0' recommended)

% Estimate freestream conditions
strucObs.U_Inf.estimate  = false;  % Estimate freestream (inflow) u_Inf and v_Inf
strucObs.U_Inf.intFactor = 0.99;  % LPF gain (1: do not change, 0: instant change)

% Measurement definitions
strucObs.measPw      = false; % Use power measurements (SCADA) from turbines in estimates
strucObs.measFlow    = true;  % Use flow measurements (LIDAR) in estimates
strucObs.sensorsPath = 'sensors_apc_6turb_alm'; % measurement setup filename (see '/setup_sensors/sensors_layouts')
    
% Kalman filter settings
strucObs.filtertype = 'sim'; % Observer types are outlined next
switch lower(strucObs.filtertype)
    
    % Extended Kalman filter (ExKF)
    case {'exkf'}
        % Covariances
        strucObs.R_k = 1.0; % Measurement   covariance matrix
        strucObs.Q_k = 1.0; % Process noise covariance matrix
        strucObs.P_0 = 0.5; % Initial state covariance matrix

        % Other model settings
        scriptOptions.exportPressures = false; % Model/predict/filter pressure terms
        scriptOptions.Linearversion   = true;  % Calculate linearized system matrices: necessary for ExKF
        
    % Sliding mode observer (SMO)
    case {'smo'}
        % tuning parameters
        strucObs.alpha = 0.01; 
        
        % Other model settings
        scriptOptions.exportPressures = false; % Model/predict/filter pressure terms
        scriptOptions.Linearversion   = true;  % Calculate linearized system matrices: necessary for SMO
        
        
    % Unscented Kalman filter (UKF)
    case {'ukf'}
        % General settings
        strucObs.stateEst             = true;  % Do state estimation: true/false
        scriptOptions.exportPressures = false; % Model, predict and filter pressure terms
        
        % Covariances
        strucObs.R_k   = 0.10;  % Measurement   covariance matrix
        strucObs.R_ePW = 1e-3;  % Measurement noise for turbine power measurements        
        strucObs.Q_k.u = 0.10;  % Process noise covariance matrix
        strucObs.Q_k.v = 0.01;  % Process noise covariance matrix
        strucObs.Q_k.p = 0.0;   % Process noise covariance matrix        
        strucObs.P_0.u = 0.10;  % Initial state covariance matrix
        strucObs.P_0.v = 0.10;  % Initial state covariance matrix
        strucObs.P_0.p = 0.0;   % Initial state covariance matrix      

        % Online model parameter adaption/estimation/tuning
        strucObs.tune.est  = false; % Do estimation
        strucObs.tune.vars = {'turbine.forcescale','site.lmu'}; % If empty {} then no estimation
        strucObs.tune.Q_k  = [1e-3,1e-2]; % Standard dev. for process noise 'u' in m/s
        strucObs.tune.P_0  = [1e-3,1e-2]; % Width of uniform dist. around opt. estimate for initial ensemble
        strucObs.tune.lb   = [0.00,0.00]; % Lower bound
        strucObs.tune.ub   = [6.00,6.00]; % Upper bound
        
        % Sigma-point generation settings
        strucObs.alpha = 1e0;
        strucObs.beta  = 2; % 2 is optimal for Gaussian distributions
        strucObs.kappa = 0;% "0" or "3-L"
        
        % Other model settings
        scriptOptions.Linearversion = false;   % Calculate linearized system matrices
        

    % Ensemble Kalman filter (EnKF)    
    case {'enkf'}
        % General settings
        strucObs.nrens = 50; % Ensemble size
        scriptOptions.exportPressures = false; % Include pressure terms in ensemble members (default: false)
        
        % Model state covariances
        strucObs.stateEst = true;  % Estimate model states
        strucObs.R_e   = 0.10; % Standard dev. for measurement noise ensemble
        strucObs.R_ePw = 5e3;  % Measurement noise for turbine power measurements        
        strucObs.Q_e.u = 0.10; % Standard dev. for process noise 'u' in m/s
        strucObs.Q_e.v = 0.01; % Standard dev. for process noise 'v' in m/s
        strucObs.Q_e.p = 0.00; % Standard dev. for process noise 'p' in m/s        
        strucObs.W_0.u = 0.90; % Width (in m/s) of uniform dist. around opt. estimate for initial ensemble
        strucObs.W_0.v = 0.30; % Width (in m/s) of uniform dist. around opt. estimate for initial ensemble
        strucObs.W_0.p = 0.00; % Only used for case Projection = 0
        
        % Inflation and localization
        strucObs.r_infl = 1.025;  % Covariance inflation factor (typically 1.00-1.20, no inflation: 1)
        strucObs.f_locl = 'gaspari'; % Localization method: 'off', 'gaspari' (Gaspari-Cohn 1999) or 'heaviside' (Heaviside step function: 0s or 1s)
        strucObs.l_locl = 131;    % Gaspari-Cohn: typically sqrt(10/3)*L with L the cut-off length. Heaviside: cut-off length (m).
        
        % Parameter estimation settings
        strucObs.tune.est  = false; % Estimate model parameters
        strucObs.tune.vars = {'turbine.forcescale','site.lmu'};
        strucObs.tune.Q_e  = [0.01,0.01]; % Standard dev. for process noise 'u' in m/s
        strucObs.tune.W_0  = [0.15,0.10]; % Width of uniform dist. around opt. estimate for initial ensemble
        strucObs.tune.lb   = [0.20,0.50]; % Lower bounds
        strucObs.tune.ub   = [2.50,2.50]; % Upper bounds
                
        % Other settings
        scriptOptions.Linearversion = false; % Disable unnecessary calculations in model
        
        
    % No filter, just open-loop simulation
    case {'sim'}
        scriptOptions.exportPressures = true;  % Must be 'true' for sim case.
        scriptOptions.Linearversion   = false; % Calculate linearized system matrices
        strucObs.measFlow             = false; % Use flow measurements (LIDAR) in estimates
        strucObs.measPw               = false; % Use power measurements (SCADA) from turbines in estimates
        
    otherwise
        error('not a valid filter/simulation specified.');
end