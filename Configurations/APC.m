% SOWFA source directories and meshing options
sourcepath      = 'WFSim\Data_SOWFA\APC\9turb_100x50_lin';        % Specify location of SOWFA data (excluding backslash at the end)
datanroffset    = 22000;                                          % Numbering offset (first filenumber is datanroffset+1)
Wp.name         = 'APC_3x3turb_noyaw_9turb_100x50_lin';           % Name of meshing (from "meshing.m")
sensors_path    = ['Setup_sensors\sensor_layouts\' ...            % Specify file with sensor (measurement) locations
                   'sensors_downwind_APC_3x3turb_noyaw_9turb_100x50_lin.mat']; 

% Model settings
options.startUniform    = 1;    % Start from a uniform flow field (1) or from a fully developed waked flow field (0).
conv_eps                = 1e-6; % Convergence parameter
max_it_dyn              = 1;    % Convergence parameter
max_it                  = 1;    % Convergence parameter

% Filter settings
strucObs.filtertype      = 'enkf';  % Observer types are outlined below in "Filter settings"
strucObs.obsv_delay      = 000;    % Number of time steps after which the observer is enabled (between 0 and NN-1)
strucObs.loadrandomseed  = 1;      % Load a predefined random seed (for one-to-one comparisons between simulation cases)
strucObs.noise_obs       = 0.1;    % Disturbance amplitude (m/s) in output data by randn*noiseampl ('0' for no noise)
strucObs.noise_init      = 0.0;    % Disturbance amplitude (m/s) in initial flow field by randn*noiseinit ('0' recommended)
strucObs.noise_input     = 0.0;    % Noise on input vector beta, enforced by the command "randn*beta"

switch lower(strucObs.filtertype)
    case {'exkf'}
        strucObs.R_k            = 1.0; % Measurement   covariance matrix   
        strucObs.Q_k            = 1.0; % Process noise covariance matrix 
        strucObs.P_0            = 0.5; % Initial state covariance matrix 
        options.exportPressures = 0;   % Model/predict/filter pressure terms
        
        strucObs.diagP        = true;   % Neglect all off-diagonal elements in P
        strucObs.sparseF      = true;   % Sparsify F matrix to reduce number of operations
        strucObs.Fthresh      = 0.01;   % Neglect values smaller than [*] in F (if above is set to true)
        options.Linearversion = 1;      % Calculate linearized system matrices    
        
    case {'enkf'}
        strucObs.nrens   =   50;         % Ensemble size
        strucObs.R_e     =   0.10;       % Standard dev. for measurement noise ensemble
        strucObs.Q_e.u   =   0.08;       % Standard dev. for process noise 'u' in m/s
        strucObs.Q_e.v   =   0.02;       % Standard dev. for process noise 'v' in m/s
        strucObs.Q_e.p   =   0.00;       % Standard dev. for process noise 'p' in m/s
        strucObs.W_0.u   =   0.90;       % Width (in m/s) of uniform dist. around opt. estimate for initial ensemble
        strucObs.W_0.v   =   0.30;       % Width (in m/s) of uniform dist. around opt. estimate for initial ensemble
        strucObs.W_0.p   =   0.00;       % Only used for case Projection = 0 
        
        strucObs.measPw       = false;    % Use power measurements from turbines in estimates
        strucObs.usePwforFlow = false;    % Have direct correlation between states and measured Pw (recommended: off)
          strucObs.pwLocFactor.auto  = 1; % Correction factor between 0 (uncorrelated) and 1 (no correction): covar. entries are multiplied with this to discouple power and states
          strucObs.pwLocFactor.cross = 1; % Correction factor between 0 (uncorrelated) and 1 (no correction): covar. entries are multiplied with this to discouple power and states
        strucObs.R_ePW        = 1e-3;     % Measurement noise for turbine power measurements
        
        options.exportPressures = false;      % Include pressure terms in ensemble members (default: false)
        strucObs.r_infl         = 1;          % Covariance inflation factor (typically 1.00-1.20, no inflation: 1)
        strucObs.f_locl         = 'gaspari';  % Localization method: 'off', 'gaspari' (Gaspari-Cohn 1999) or 'heaviside' (Heaviside step function: 0s or 1s)
        strucObs.l_locl         = 50;         % Gaspari-Cohn: typically sqrt(10/3)*L with L the cut-off length. Heaviside: cut-off length (m).
        options.Linearversion   = 0;          % Calculate linearized system matrices. Do not change, keep this '0'.      
        
    case {'sim'}
        options.exportPressures = 1; % Do not change for sim case.
        options.Linearversion = 0; % Calculate linearized system matrices     
        
    otherwise
        error('not a valid filter/simulation specified.');
end;