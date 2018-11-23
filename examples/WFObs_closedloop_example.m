clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%             WIND FARM OBSERVER (WFObs) by B.M. Doekemeijer
%                 Delft University of Technology, 2018
%          Repo: https://github.com/TUDelft-DataDrivenControl/WFObs
%
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%
%%   Quick use:
%     1. Specify your script preferences in lines 58 - 73
%     3. Specify your simulation settings file in line  76
%      a. if you want (OPTIONAL!), you can create a simulation case
%         manually, then have a look at the ./configurations folder
%          i.  Use  '/Setup_sensors/WFObs_s_sensors.m' to specify measurement locations
%          ii. All WFSim model/site information is contained in '../WFSim/bin/core/meshing.m'
%     3. Press start.
%
%%   Relevant input/output variables
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
%%   Debugging and contributing:
%     - First, try to locate any errors by turning all possible outputs
%       on (printProgress, printConvergence, Animate, plotMesh).
%     - If you cannot solve your problems, reach out on the Github.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set-up WFSim model
addpath('../WFSim/layoutDefinitions') % Folder with predefined wind farm layouts
Wp = layoutSet_clwindcon_3turb(); % Choose which scenario to simulate. See 'layoutDefinitions' folder for the full list.
Wp.turbine.Cry(3) = 480; % Move turbine a little closer to aligned
Wp.sim.h = 5;
addpath('../WFSim/solverDefinitions'); % Folder with model options, solver settings, etc.
modelOptions = solverSet_default(Wp); % Choose model solver options. Default for EnKF/UKF: solverSet_minimal. For ExKF: solverSet_linearmatrices.

%% Setup KF settings
addpath('../filterDefinitions') % Folder with predefined KF settings
strucObs = filterSet_openloop(); % Observer/KF settings

%% Setup sensors
addpath('../sensorDefinitions')
measOptions = sensorSet_power_only(Wp);

%% Setup WFSim as simulation model
addpath('../WFSim/bin/core');
simModel.Wp = Wp; % Can differ from Wp to mimic model discrepancies
% simModel.Wp.site.lm_slope = 0.8; % Purposely introduce error with controller model
% simModel.Wp.site.u_Inf = Wp.site.u_Inf - 1.0; % Purposely introduce error with controller model
simModel.modelOptions = solverSet_default(Wp); % Choose model solver options. Default for EnKF/UKF: solverSet_minimal. For ExKF: solverSet_linearmatrices.
[simModel.Wp,simModel.sol,simModel.sys] = InitWFSim(simModel.Wp,simModel.modelOptions,0); % Initialize WFSim model

% External animations for this script specifically
scriptOptions.Animate     = 12; % Animation frequency
scriptOptions.plotContour = 1;  % Show flow fields
scriptOptions.plotPower   = 0;  % Plot true and predicted power capture vs. time
scriptOptions.plotError   = 0;  % plot RMS and maximum error vs. time
scriptOptions.savePlots   = 0;  % Save all plots in external files at each time step
scriptOptions.savePath    = ['results/tmp']; % Destination folder of saved files


%% Initialize WFObs object
addpath('../bin','../bin_supplementary/online_wfsim'); % Add the main 'bin' folder and the postProcessing folder
addpath('../bin_supplementary/offline_tools');
WFObj = WFObs_obj( Wp,modelOptions,strucObs ); % Initialize WFObj object

%% Execute the WFObs core code
hFigs = [];
optimizationFrequency = 60; % Every [x] seconds

% Turbine input
inputData = struct('phi',zeros(simModel.Wp.turbine.N,1),...
    'CT_prime',2*ones(simModel.Wp.turbine.N,1));
                   
while 1
    % Random event after 5 minutes: sudden change in TI
    if simModel.sol.time == 300
        disp('DISCRETE EVENT AT T=300 s: Simulating a change in TI...')
        simModel.Wp.site.lm_slope = 0.21;   
        WFObj.model.Wp.site.lm_slope = 0.21;
    end
    
    % Grab turbine inputs and measurements from simulation model
    [simModel.sol,simModel.sys] = WFSim_timestepping(simModel.sol,simModel.sys,simModel.Wp,inputData,simModel.modelOptions); % forward timestep: x_k+1 = f(x_k)
    [measuredData] = getDataWFSim(simModel,measOptions);
    
    % Perform estimation
    sol = WFObj.timestepping(inputData,measuredData);
    
    if rem(WFObj.model.sol.time,optimizationFrequency) == 0
        tic
        Np = Inf; % Prediction horizon in seconds
        [optYawAngles] = optimizeStaticYaw(WFObj.model,Np)
        inputData.phi = optYawAngles';
        toc
    end
    
    % Post-processing
    solTrue = struct('u',simModel.sol.u,'v',simModel.sol.v,'P',simModel.sol.turbine.power);
    sol_array(sol.k) = formatSol(WFObj.model,solTrue); % Save reduced-size solution to an array
    [ hFigs,scriptOptions ] = WFObs_p_animations( WFObj.strucObs,WFObj.model.Wp,sol_array,scriptOptions,hFigs ); % Create figures
end


% Optimization function
function [optYawAngles] = optimizeStaticYaw(model,Np)
    
    % Define cost function
    J = @(x) power_predhorizon(x,model,Np);

    % Grid search optimization
    costOpt = -1; % Initial optimal cost
    for gamma1 = max([-40,model.sol.turbInput.phi(1)-10]):5:min([40, model.sol.turbInput.phi(1)+10])
        gamma2 = 0;
%         for gamma2 = -20:10:20
            x = [gamma1 gamma2 0];
            cost = power_predhorizon(x',model,Np);
            if cost > costOpt
                costOpt = cost;
                xOpt = x;
            end
%         end
    end
    optYawAngles = xOpt;
    
    
    function cost = power_predhorizon(xYaw,model,Np)
        Wp = model.Wp;
        options = model.modelOptions;
        nTurbs = length(Wp.turbine.Crx);
        
        % Create input set
        turbInput = struct('phi',xYaw,'CT_prime',2*ones(nTurbs,1));
        
        % Evaluate result
        if Np == Inf % Infinite horizon: use steady-state
            Wp.sim.h = Inf;
            [Wp,sol,sys] = InitWFSim(Wp,options,0);
            [options.max_it_dyn,options.max_it] = deal(100);
            [sol,sys] = WFSim_timestepping(sol,sys,Wp,turbInput,options);
            [sol,~]   = WFSim_timestepping(sol,sys,Wp,turbInput,options);
            cost = sum(sol.turbine.power); % Time-averaged power
        else
            cost = 0;
            sol = model.sol;
            sys = model.sys;
            initTime = sol.time;
            while sol.time - initTime < Np % Predict over finite horizon
                sol = WFSim_timestepping(sol,sys,Wp,turbInput,options);
                cost = cost + sum(sol.turbine.power)*Wp.sim.h; % Integrated energy
            end
        end
    end
end
