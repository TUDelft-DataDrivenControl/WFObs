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
addpath('WFSim/layoutDefinitions') % Folder with predefined wind farm layouts
Wp = layoutSet_sowfa_9turb_apc_alm_turbl(); % Choose which scenario to simulate. See 'layoutDefinitions' folder for the full list.
addpath('WFSim/solverDefinitions'); % Folder with model options, solver settings, etc.
modelOptions = solverSet_minimal(Wp); % Choose model solver options. Default for EnKF/UKF: solverSet_minimal. For ExKF: solverSet_linearmatrices.

%% Setup KF settings
addpath('filterDefinitions') % Folder with predefined KF settings
strucObs = filterSet_WES2018(); % Observer/KF settings

%% Setup sensors
addpath('sensorDefinitions')
measOptions = sensorSet_power_only(Wp);

%% Setup WFSim as simulation model
simModel.Wp = Wp; % Can differ from Wp to mimic model discrepancies
simModel.modelOptions = solverSet_default(Wp); % Choose model solver options. Default for EnKF/UKF: solverSet_minimal. For ExKF: solverSet_linearmatrices.
[simModel.Wp,simModel.sol,simModel.sys] = InitWFSim(simModel.Wp,simModel.modelOptions,0); % Initialize WFSim model

% External animations for this script specifically
scriptOptions.Animate     = 5; % Animation frequency
scriptOptions.plotContour = 1;  % Show flow fields
scriptOptions.plotPower   = 1;  % Plot true and predicted power capture vs. time
scriptOptions.plotError   = 1;  % plot RMS and maximum error vs. time
scriptOptions.savePlots   = 0;  % Save all plots in external files at each time step
scriptOptions.savePath    = ['results/tmp']; % Destination folder of saved files


%% Initialize WFObs object
addpath('bin','bin_supplementary/online_wfsim'); % Add the main 'bin' folder and the postProcessing folder
WFObj = WFObs_obj( Wp,modelOptions,strucObs ); % Initialize WFObj object

%% Execute the WFObs core code
hFigs = [];
while WFObj.model.sol.time < 1000 % While measurements available

    % Turbine input
    inputData = struct('phi',zeros(simModel.Wp.turbine.N,1),...
                       'CT_prime',2*ones(simModel.Wp.turbine.N,1));
    
    % Grab turbine inputs and measurements from simulation model
    [simModel.sol,simModel.sys] = WFSim_timestepping(simModel.sol,simModel.sys,simModel.Wp,inputData,simModel.modelOptions); % forward timestep: x_k+1 = f(x_k)
    [measuredData] = getDataWFSim(simModel,measOptions);
    
    % Perform estimation
    sol = WFObj.timestepping(inputData,measuredData);
    
    % Post-processing
    solTrue = struct('u',simModel.sol.u,'v',simModel.sol.v,'P',simModel.sol.turbine.power);
    sol_array(sol.k) = formatSol(WFObj.model,solTrue); % Save reduced-size solution to an array
    [ hFigs,scriptOptions ] = WFObs_p_animations( WFObj.strucObs,WFObj.model.Wp,sol_array,scriptOptions,hFigs ); % Create figures
end