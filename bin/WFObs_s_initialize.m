function [ Wp,sol,sys,strucObs,scriptOptions ] = WFObs_s_initialize( scriptOptions,configName )
% WFOBS_S_INITIALIZE  Initialize the WFSim model and the estimator settings
%
%   SUMMARY
%    This code does the necessary initializations for the WFSim model, for
%    the estimator, and for the relevant script settings.
%
%   RELEVANT INPUT/OUTPUT VARIABLES
%     - configName: name of the simulation case that is to be simulated.
%     All simulation scenarios can be found in the '/configurations/'
%     folder as seperate files. The default case is 'YawCase3.m'.
%
%     - scriptOptions: this struct contains all simulation settings, not
%     related to the wind farm itself (solution methodology, outputs, etc.)
%
%     - Wp: this struct contains all the simulation settings related to the
%           wind farm, the turbine inputs, the atmospheric properties, etc.
%
%     - sol: this struct contains the system states at a certain timestep.
%         sol.u:     Instantaneous longitudinal flow field over the mesh (in m/s)
%         sol.v:     Instantaneous longitudinal flow field over the mesh (in m/s)
%         sol.uu:    Same as sol.u, used for convergence
%         sol.vv:    Same as sol.v, used for convergence
%
%     - sys: this struct contains the system matrices at a certain timestep.
%         sys.pRCM:  Reverse Cuthill-McKee algorithm for solving A*x=b faster.
%         sys.B1:    Important matrix in the boundary conditions.
%         sys.B2:    Important matrix in the boundary conditions.
%         sys.bc:    Important vector in the boundary conditions.
%
%

% Load configuration file from the 'configurations' folder
run(configName);    

% Decide whether we need to linearize WFSim
if strcmp(lower(strucObs.filtertype),'exkf') || strcmp(lower(strucObs.filtertype),'smo')
    % Calculate linearized system matrices: necessary for ExKF & SMO
    scriptOptions.Linearversion = true; 
else
    % Disable calculation of lin. matrices: not needed for EnKF & UKF
    scriptOptions.Linearversion = false; 
end

% Check KF settings compatibility
if strcmp(lower(strucObs.filtertype),'sim') == false
    if strucObs.se.enabled == 0 && strucObs.pe.enabled == 0
        error(['Please turn on state and/or parameter estimation. '...
            'Alternatively, select "sim" for open-loop simulations.']);
    end
end
    
if scriptOptions.printProgress
    disp(' WindFarmObserver (WFObs)');
    disp([' Case:  ' configName ]);
    disp(' ');
end

% load a default random seed for consistency
if strucObs.loadRandomSeed
    load('randomseed')
    rng(randomseed)
    clear randomseed;
end

% Default settings: following WFSim options are never used in WFObs
scriptOptions.Projection      = 0;    % Use projection
scriptOptions.exportLinearSol = 0;    % Export linear solution
scriptOptions.Derivatives     = 0;    % Calculate derivatives/gradients for system

if scriptOptions.printProgress
    disp([datestr(rem(now,1)) ' __  Initializing simulation model.']);
end;

[Wp,sol,sys] = InitWFSim(Wp,scriptOptions); % Initialize model

% Set default model convergence settings
scriptOptions.conv_eps     = 1e-6; % Convergence parameter
scriptOptions.max_it_dyn   = 1;    % Convergence parameter

if Wp.sim.startUniform
    scriptOptions.max_it = 1;   % Iteration limit for simulation start-up
else
    scriptOptions.max_it = 50;  % Iteration limit for simulation start-up
end

% Define what the system should predict (with or without pressures)
strucObs.size_state = Wp.Nu + Wp.Nv + Wp.Np;
if scriptOptions.exportPressures == 0
    strucObs.size_output = Wp.Nu + Wp.Nv;
else
    strucObs.size_output = Wp.Nu + Wp.Nv + Wp.Np;
end;

% Create global RCM vector
[~, sysRCM] = WFSim_timestepping( sol, sys, Wp, scriptOptions );
sys.pRCM    = sysRCM.pRCM;

scriptOptions.klen = length(num2str(Wp.sim.NN));        % used for proper spacing in cmd output window
scriptOptions.tlen = length(num2str(Wp.sim.time(end))); % length

if scriptOptions.printProgress
    disp([datestr(rem(now,1)) ' __  Finished initialization sequence.']);
end;
end