function [ Wp,sol,sys,strucObs ] = WFObs_s_initialize( Wp, strucObs, modelOptions, verboseOptions )
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

% Check KF settings compatibility
if strcmp(lower(strucObs.filtertype),'sim') == false
    if strucObs.se.enabled == 0 && strucObs.pe.enabled == 0
        error(['Please turn on state and/or parameter estimation. '...
            'Alternatively, select "sim" for open-loop simulations.']);
    end
end
    
if verboseOptions.printProgress
    disp(' WindFarmObserver (WFObs)');
    disp(' ');
end

% load a default random seed for consistency
if strucObs.loadRandomSeed
    load('randomseed')
    rng(randomseed)
    clear randomseed;
end

if verboseOptions.printProgress
    disp([datestr(rem(now,1)) ' __  Initializing simulation model.']);
end;

[Wp,sol,sys] = InitWFSim(Wp,modelOptions,verboseOptions.plotMesh); % Initialize model

% Define what the system should predict (with or without pressures)
strucObs.size_state = Wp.Nu + Wp.Nv + Wp.Np;
if modelOptions.exportPressures == 0
    strucObs.size_output = Wp.Nu + Wp.Nv;
else
    strucObs.size_output = Wp.Nu + Wp.Nv + Wp.Np;
end;

% Create global RCM vector
turbInput = struct('t',0,'CT_prime',2*ones(Wp.turbine.N,1),'phi',zeros(Wp.turbine.N,1));
[~, sysRCM] = WFSim_timestepping( sol, sys, Wp, turbInput, modelOptions );
sys.pRCM    = sysRCM.pRCM;

if verboseOptions.printProgress
    disp([datestr(rem(now,1)) ' __  Finished initialization sequence.']);
end;
end