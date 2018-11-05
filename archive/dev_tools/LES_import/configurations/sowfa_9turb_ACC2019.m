% Source files
scriptOptions.outputFilename = 'sowfa_9turb_ACC2019.mat';
scriptOptions.plotFrequency  = 10;         % Plot mapping every * instances (will always plot k == 1, set to high value for no plots after k == 1)
scriptOptions.sourcePath     = 'W:\OpenFOAM\bmdoekemeijer-2.4.0\simulationCases\ACC2019\NREL5MW_9turb_SSC_aligned_CL_deterministic\postProcessing'

% Turbine properties directly from PALM or SOWFA. The reference frame is 
%   x (vertical, upwards pos.) - y (horizontal, rightwards pos.).
rawTurbData           = struct('Crx',[868 868 868 1500 1500 1500 2132 2132 2132],...
                               'Cry',repmat([1120.8 1500.0 1879.2],1,3));
rawTurbData.Drotor    = 126.4*ones(1,9); % Rotor diameter in (m)
rawTurbData.hubHeight = 90.0             % Hub height in (m)

% Filtering
filterSettings.turbData.MM = true; % Apply moving average to turbine data
filterSettings.turbData.tL = 1;    % Window width to the left (seconds)
filterSettings.turbData.tR = 1;    % Window width to the right (seconds)

filterSettings.Ur.nPts = 50;    % Number of rotor points (only if CT not available)
filterSettings.Ur.MM   = true;  % Moving-average for rotor velocity (Only when CT not given)
filterSettings.Ur.tL   = 5;     % Moving-average for rotor velocity (Only when CT not given)
filterSettings.Ur.tR   = 5;     % Moving-average for rotor velocity (Only when CT not given)

filterSettings.CTp.MM   = true;  % Additional moving-mean average for Ct_prime
filterSettings.CTp.tL   = 3;     % Additional moving-mean average for Ct_prime
filterSettings.CTp.tR   = 3      % Additional moving-mean average for Ct_prime

% Desired output settings
meshSetup.dt          = 5.0 ; % Timestep in seconds
meshSetup.rho         = 1.20; % Air density (kg m^-3)
meshSetup.distance_S  = 600 ; % distance (m) upwind   first  turbine to export
meshSetup.distance_N  = 650;  % distance (m) downwind  last  turbine to export
meshSetup.distance_W  = 600 ; % distance (m) west most left  turbine (from hub) to export
meshSetup.distance_E  = 600 ; % distance (m) east most right turbine (from hub) to export
meshSetup.Nx          = 50;   % Number of grid points in x-direction (-)
meshSetup.Ny          = 25    % Number of grid points in y-direction (-)