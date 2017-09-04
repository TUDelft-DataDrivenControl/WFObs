%% SOWFA source directories, meshing and measurement options
Wp.name = '2turb_adm_noturb';  % Name of meshing (from '/WFSim/bin/core/meshing.m')

%% WFSim model settings
scriptOptions.startUniform    = 0;    % Start from a uniform flow field (1) or from a fully developed waked flow field (0).
scriptOptions.conv_eps        = 1e-6; % Convergence parameter
scriptOptions.max_it_dyn      = 1;    % Convergence parameter

if scriptOptions.startUniform==1
    scriptOptions.max_it = 1;   % Iteration limit for simulation start-up
else
    scriptOptions.max_it = 50;  % Iteration limit for simulation start-up
end


%% Filter settings
% General settings
strucObs.loadRandomSeed  = 1;      % Load a predefined random seed (for one-to-one comparisons between simulation cases)
strucObs.noise_obs       = 0.1;    % Disturbance amplitude (m/s) in output data by randn*noiseampl ('0' for no noise)
strucObs.noise_init      = 0.0;    % Disturbance amplitude (m/s) in initial flow field by randn*noiseinit ('0' recommended)

% DISABLED FOR NOW:
% % Estimate freestream conditions
% strucObs.U_Inf.estimate  = false;  % Estimate freestream (inflow) u_Inf and v_Inf
% strucObs.U_Inf.intFactor = 0.99;  % LPF gain (1: do not change, 0: instant change)

% Measurement definitions
strucObs.measPw      = true;  % Use power measurements (SCADA) from turbines in estimates
strucObs.measFlow    = false; % Use flow measurements (LIDAR) in estimates
strucObs.sensorsPath = 'sensors_2turb_adm_wakecenter'; % measurement setup filename (see '/setup_sensors/sensors_layouts')
        
% Kalman filter settings
strucObs.filtertype      = 'enkf'; % Observer types are outlined next
switch lower(strucObs.filtertype)
    
    % Extended Kalman filter (ExKF)
    case {'exkf'}
        % Covariances
        strucObs.R_k = 1.0; % Measurement   covariance matrix
        strucObs.Q_k = 1.0; % Process noise covariance matrix
        strucObs.P_0 = 0.5; % Initial state covariance matrix
        
        % Simplifications/speed-ups
        strucObs.diagP   = true;    % Neglect all off-diagonal elements in P
        strucObs.sparseF = true;    % Sparsify F matrix to reduce number of operations
        strucObs.Fthresh = 0.01;    % Neglect values smaller than [*] in F (if above is set to true)
        
        % Other model settings
        scriptOptions.exportPressures = 0; % Model/predict/filter pressure terms
        scriptOptions.Linearversion   = 1; % Calculate linearized system matrices: necessary for ExKF
        
        
    % Unscented Kalman filter (UKF)
    case {'ukf'}
        % General settings
        strucObs.stateEst             = false;  % Do state estimation: true/false
        scriptOptions.exportPressures = 0;  % Model, predict and filter pressure terms
        
        % Covariances
        strucObs.R_k   = 0.10;  % Measurement   covariance matrix
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
        
        % Measurement definitions
        strucObs.R_ePw = 5e3;  % INCOMPATIBLE [Measurement noise for turbine power measurements]
        strucObs.R_e   = 0.10; % Standard dev. for measurement noise ensemble

        % Other model settings
        scriptOptions.Linearversion   = 0;   % Calculate linearized system matrices
        

    % Ensemble Kalman filter (EnKF)    
    case {'enkf'}
        % General settings
        strucObs.nrens      = 50; % Ensemble size
        strucObs.resampling = 0;  % Redistribute particles every timestep. false = classical EnKF.
        scriptOptions.exportPressures = 0; % Include pressure terms in ensemble members (default: false)
        
        % Model state covariances
        strucObs.stateEst = false;  % Estimate model states
        strucObs.Q_e.u    = 0.10; % 1e-3;% 0.10 % Standard dev. for process noise 'u' in m/s
        strucObs.Q_e.v    = 0.01; % 1e-4;% 0.01 % Standard dev. for process noise 'v' in m/s
        strucObs.Q_e.p    = 0.00;  % Standard dev. for process noise 'p' in m/s        
        strucObs.W_0.u    = 0.90; % Width (in m/s) of uniform dist. around opt. estimate for initial ensemble
        strucObs.W_0.v    = 0.30; % Width (in m/s) of uniform dist. around opt. estimate for initial ensemble
        strucObs.W_0.p    = 0.00;  % Only used for case Projection = 0
        
        % Inflation and localization
        strucObs.r_infl         = 1.0;%1.025;  % Covariance inflation factor (typically 1.00-1.20, no inflation: 1)
        strucObs.f_locl         = 'gaspari'; % Localization method: 'off', 'gaspari' (Gaspari-Cohn 1999) or 'heaviside' (Heaviside step function: 0s or 1s)
        strucObs.l_locl         = 131;    % Gaspari-Cohn: typically sqrt(10/3)*L with L the cut-off length. Heaviside: cut-off length (m).
        
        % Parameter estimation settings
        strucObs.tune.est  = true; % Estimate model parameters
        strucObs.tune.vars = {'site.lmu'}; %{'turbine.forcescale','site.lmu'};
        strucObs.tune.Q_e  = [0.15]; %[1e-3,1e-3]; % Standard dev. for process noise 'u' in m/s
        strucObs.tune.W_0  = [0.18]; %[0.15,0.10]; % Width of uniform dist. around opt. estimate for initial ensemble
        strucObs.tune.lb   = [0.05]; %[0.40,0.05]; % Lower bounds
        strucObs.tune.ub   = [15.00]; %[3.00,2.50]; % Upper bounds
        
        % Measurement definitions
        strucObs.R_ePw       = 5e3;   % Measurement noise for turbine power measurements
        strucObs.R_e         = 0.10; % Standard dev. for measurement noise ensemble
        
        % Other settings
        scriptOptions.Linearversion   = 0; % Disable unnecessary calculations in model
        

    % No filter, just open-loop simulation
    case {'sim'}
        scriptOptions.exportPressures = 1; % Must be 'true' for sim case.
        scriptOptions.Linearversion   = 0; % Calculate linearized system matrices
        strucObs.measFlow             = false; % Use flow measurements (LIDAR) in estimates
        strucObs.measPw               = false; % Use power measurements (SCADA) from turbines in estimates
        
    otherwise
        error('not a valid filter/simulation specified.');
end