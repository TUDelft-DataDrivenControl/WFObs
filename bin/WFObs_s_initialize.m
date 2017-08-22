function [ Wp,sol,sys,strucObs,scriptOptions,LESData,hFigs ] = WFObs_s_initialize( scriptOptions,configName )
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
%     - LESData: this struct contains all the flow fields and turbine data
%                from the LES data.
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
%     - hFigs: cell array of Figures to (re)plot figures into.
%

% Load configuration file from the 'configurations' folder
run(configName);    

if scriptOptions.printProgress
    disp(' WindFarmObserver (WFObs)');
    disp([' Case:  ' configName ]);
    disp(' ');
end;

% load a default random seed for consistency
if strucObs.loadRandomSeed; load('randomseed'); rng(randomseed); clear randomseed; end;

% Create destination folder for output files
if (scriptOptions.savePlots + scriptOptions.saveWorkspace > 0)
    mkdir(scriptOptions.savePath);
end;

% Default settings: following WFSim options are never used in WFObs
scriptOptions.Projection      = 0;    % Use projection
scriptOptions.exportLinearSol = 0;    % Export linear solution
scriptOptions.Derivatives     = 0;    % Calculate derivatives/gradients for system

if scriptOptions.printProgress
    disp([datestr(rem(now,1)) ' __  Initializing simulation model.']);
end;

[Wp,sol,sys] = InitWFSim(Wp,scriptOptions); % Initialize model

% Add noise to initial conditions
[sol.u,sol.uu]  = deal(sol.u + randn(Wp.mesh.Nx,Wp.mesh.Ny)*strucObs.noise_init);
[sol.v,sol.vv]  = deal(sol.v + randn(Wp.mesh.Nx,Wp.mesh.Ny)*strucObs.noise_init);

% Define what the system should predict (with or without pressures)
strucObs.size_state = Wp.Nu + Wp.Nv + Wp.Np;
if scriptOptions.exportPressures == 0
    strucObs.size_output = Wp.Nu + Wp.Nv;
else
    strucObs.size_output = Wp.Nu + Wp.Nv + Wp.Np;
end;

% Define measurement locations
sensorsfile        = load(strucObs.sensorsPath);
strucObs.obs_array = unique([sensorsfile.sensors{1}.obsid; sensorsfile.sensors{2}.obsid]);

% Load measurements from LES simulation (*.mat file)
LESData    = load(Wp.sim.measurementFile); % Load measurements
LESData.ud = LESData.u + strucObs.noise_obs*randn(size(LESData.u)); % Add noise
LESData.vd = LESData.v + strucObs.noise_obs*randn(size(LESData.v)); % Add noise

% Setup blank figure windows
hFigs = {};

% Create global RCM vector
% soltemp     = sol; soltemp.k = 1;
% [~, sysRCM] = Make_Ax_b(Wp,sys,soltemp,scriptOptions);
% sys.pRCM    = sysRCM.pRCM;  clear sysRCM soltemp;

scriptOptions.klen = length(num2str(Wp.sim.NN));        % used for proper spacing in cmd output window
scriptOptions.tlen = length(num2str(Wp.sim.time(end))); % length

% Save simulation & filter settings
if (scriptOptions.savePlots + scriptOptions.saveWorkspace > 0)
    save([scriptOptions.savePath '/' strucObs.filtertype '_settings.mat']);
end;

if scriptOptions.printProgress
    disp([datestr(rem(now,1)) ' __  Finished initialization sequence.']);
end;
end