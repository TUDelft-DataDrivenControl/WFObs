% Source files
scriptOptions.outputFilename = 'APC';
scriptOptions.plotFrequency  = 1e9;         % Plot mapping every * instances (will always plot k == 1, set to high value for no plots after k == 1)
scriptOptions.saveMemory     = false;       % turn on if you are having memory issues (SOWFA data only)
scriptOptions.sourcePath     = '/tudelft/ls/staff-group/3me/dcsc/DataDriven/Data/SOWFA/APC_12_wake_noyaw_pitch_agc_50/sliceDataInstant';

% Turbine properties directly from PALM or SOWFA. The reference frame is 
%   x (vertical, upwards pos.) - y (horizontal, rightwards pos.).

tLocs = [1086.7, 1554.9 ; ... % turbine 1 % APC 9 turbine case
         1634.0, 1870.9 ; ... % turbine 2
         1465.9, 898.1  ; ... % turbine 3
         1276.3, 1226.5 ; ... % turbine 4
         2013.2, 1214.1 ; ... % turbine 5
         1823.6, 1542.5 ; ... % turbine 6
         2560.5, 1530.1 ; ... % turbine 7
         2370.9, 1858.5 ; ... % turbine 8
         2181.3, 2186.9 ];    % turbine 9
rawTurbData           = struct('Crx',tLocs(:,2),'Cry',tLocs(:,2));			   
rawTurbData.Drotor    = 126.4*ones(1,9); % Rotor diameter in (m)
rawTurbData.hubHeight = 90.0;            % Hub height in (m)

% Desired output settings
meshSetup.dt          = 1.0 ; % Timestep in seconds
meshSetup.rho         = 1.20; % Air density (kg m^-3)
meshSetup.distance_S  = 300 ; % distance (m) upwind   first  turbine to export
meshSetup.distance_N  = 1100; % distance (m) downwind  last  turbine to export
meshSetup.distance_W  = 250 ; % distance (m) west most left  turbine (from hub) to export
meshSetup.distance_E  = 250 ; % distance (m) east most right turbine (from hub) to export
meshSetup.Nx          = 100 ; % Number of grid points in x-direction (-)
meshSetup.Ny          = 40  ; % Number of grid points in y-direction (-)