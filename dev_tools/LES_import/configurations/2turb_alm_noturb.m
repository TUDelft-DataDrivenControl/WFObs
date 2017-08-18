% Source files
scriptOptions.outputFilename = '2turb_alm_noturb';
scriptOptions.plotFrequency  = 1e9;         % Plot mapping every * instances (will always plot k == 1, set to high value for no plots after k == 1)
scriptOptions.saveMemory     = false;       % turn on if you are having memory issues (SOWFA data only)
scriptOptions.sourcePath     = '/tudelft/ls/staff-group/3me/dcsc/DataDriven/Data/SOWFA/sampleForPODexcitationPRBSpositiveNoPrecursor/horizontal_slices';

% Turbine properties directly from PALM or SOWFA. The reference frame is 
%   x (vertical, upwards pos.) - y (horizontal, rightwards pos.).
rawTurbData           = struct('Crx',[1342.0, 1658.0],'Cry',[1226.3, 1773.7]);
rawTurbData.Drotor    = [126.4, 126.4]; % Rotor diameter in (m)
rawTurbData.hubHeight = 90.0;           % Hub height in (m)
rawTurbData.tau       = 3;              % Time constant tau in low pass filter 1/(tau*s+1)

% Desired output settings
meshSetup.dt          = 1.0 ; % Timestep in seconds
meshSetup.rho         = 1.20; % Air density (kg m^-3)
meshSetup.distance_S  = 400 ; % distance (m) upwind   first  turbine to export
meshSetup.distance_N  = 850;  % distance (m) downwind  last  turbine to export
meshSetup.distance_W  = 400 ; % distance (m) west most left  turbine (from hub) to export
meshSetup.distance_E  = 400 ; % distance (m) east most right turbine (from hub) to export
meshSetup.Nx          = 50;   % Number of grid points in x-direction (-)
meshSetup.Ny          = 25;   % Number of grid points in y-direction (-)