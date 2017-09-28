function [ outputData ] = WFObs_core( scriptOptions, configName, WpOverwrite )
% WFOBS_CORE  Perform a complete time simulation with state estimation
%
%   SUMMARY
%    This code will complete a full wind farm simulation including state
%    estimation, as set up in configurations/*configName*.m.  It will use
%    measurements from data_SOWFA/* or data_PALM/* to improve the flow
%    estimations compared to the WFSim model.
%
%   RELEVANT INPUT/OUTPUT VARIABLES
%     - configName: name of the simulation case that is to be simulated.
%     All simulation scenarios can be found in the '/configurations/'
%     folder as seperate files. The default case is 'YawCase3.m'.
%
%     - scriptOptions: this struct contains all simulation settings, not
%     related to the wind farm itself (solution methodology, outputs, etc.)
%
%     - OutputData (*): this struct contains all interesting outputs. It packs
%       several substructs, namely:
%       - *.Wp: this struct contains all the simulation settings related to
%         the wind farm, the turbine inputs, the atmospheric properties, etc.
%         See WFSim.m for a more elaborate description of 'Wp'.
%
%       - *.sol_array: this cell array contains the system states at every
%         simulated time instant. Each cell entry contains a sol struct.
%         See WFSim.m for a more elaborate description of 'sol'. In
%         addition to the entries described in WFSim.m, each 'sol' struct
%         contains in addition:
%           *.sol.score: a struct containing estimation performance scores
%           such as the maximum estimation error, the RMS error, and the
%           computational cost (CPU time).
%           *.sol.measuredData: a struct containing the true (to be
%           estimated) values, and the measurement data fed into the
%           estimation algorithm.
%
%       - *.strucObs: this struct contains all the observer settings and
%       (temporary) files used for updates, such as covariance matrices,
%       ensemble/sigma point sets, measurement noise, etc.
%

%% Pre-processing
timerScript = tic;       % Start script timer

% Initialize model and observer variables
[Wp,sol,sys,strucObs,scriptOptions,LESData,hFigs] = ...
    WFObs_s_initialize(scriptOptions,configName);

max_it   = scriptOptions.max_it;    % Convergence constraints
conv_eps = scriptOptions.conv_eps;  % Convergence constraints

% Overwrite variables if WpOverwrite is specified
if nargin > 2
    if scriptOptions.printProgress
        disp('Overwriting variables in Wp...');
    end
    Wp = mergeStruct(Wp,WpOverwrite);
    [sys.B1,sys.B2,sys.bc] = Compute_B1_B2_bc(Wp); % Update boundary conditions
end

%% Core: time domain simulations
while sol.k < Wp.sim.NN
    timerCPU = tic;                 % Start iteration timer
    sol.k    = sol.k + 1;           % Timestep forward
    sol.time = Wp.sim.time(sol.k+1);% Timestep forward
    
    % Load measurement data
    sol.measuredData = WFObs_s_loadmeasurements(LESData,sol.k);
    
    % Determine freestream inflow properties from SCADA data
    [ Wp,sol,sys,strucObs ] = WFObs_s_freestream(Wp,sol,sys,strucObs);
    
    % Calculate optimal solution according to filter of choice
    [Wp,sol,strucObs] = WFObs_o(strucObs,Wp,sys,sol,scriptOptions);
    
    % Display progress in the command window
    sol = WFObs_s_reporting(timerCPU,Wp,sol,strucObs,scriptOptions);
              
    % Save reduced-size solution to an array
    sol.measuredData = rmfield(sol.measuredData,{'u','v','sol'});
    if nnz(strcmp(fieldnames(sol),'uk')) >= 1
        sol_array(sol.k) = rmfield(sol,{'uu','vv','pp','uk','vk'});
    else
        sol_array(sol.k) = rmfield(sol,{'uu','vv','pp'});
    end
    sol_array(sol.k).site = Wp.site; % Save site info too
    
    % Display animations on screen
    [hFigs,scriptOptions] = WFObs_s_animations(Wp,sol_array,sys,LESData,scriptOptions,strucObs,hFigs);
end


%% Post-processing
% save workspace variables, if necessary
if scriptOptions.saveWorkspace
    save([scriptOptions.savePath '/workspace.mat'],'configName',...
        'Wp','sys','sol_array','scriptOptions','strucObs');
end

% Put all relevant outputs in one structure
if nargout >= 1
    outputData.sol_array     = sol_array;
    outputData.Wp            = Wp;
    outputData.strucObs      = strucObs;
    outputData.scriptOptions = scriptOptions;
    outputData.configName    = configName;
end;

% Print end of simulation statement
if scriptOptions.printProgress
    disp([datestr(rem(now,1)) ' __  Completed simulations. ' ...
        'Total CPU time: ' num2str(toc(timerScript)) ' s.']);
end
end