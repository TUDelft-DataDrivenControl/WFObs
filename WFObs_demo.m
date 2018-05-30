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

%% Initialize object
WFObj=WFObs_obj('apc_9turb_alm_turb'); % See './configurations' for options

%% Preload LES data
LESData = load(WFObj.Wp.sim.measurementFile); % Load offline measurements database

%% Execute the WFObs core code
hFigs = [];
scriptOptions.Animate = 10;
    scriptOptions.plotContour    = 1;  % Show flow fields
    scriptOptions.plotPower      = 1;  % Plot true and predicted power capture vs. time
    scriptOptions.powerForecast  = 0;  % Plot power forecast (0 = disabled, x = number of steps) (only if plotPower = 1)
    scriptOptions.plotError      = 0;  % plot RMS and maximum error vs. time
    scriptOptions.savePlots      = 0;  % Save all plots in external files at each time step
    scriptOptions.savePath       = ['results/tmp']; % Destination folder of saved files

while WFObj.sol.k < WFObj.Wp.sim.NN
    % Load and format measurements from offline database
    measuredData = struct();
    for i = 1:length(WFObj.Wp.turbine.Crx)
        measuredData(i).idx  = i; % Turbine number
        measuredData(i).type = 'P'; % Power measurement
        measuredData(i).value = LESData.turbData.power(WFObj.sol.k+1,i) + ...
                                1e4*randn();
        measuredData(i).std  = 1e4; % Standard deviation in W
    end
    
    % Perform estimation
    WFObj.timestepping(measuredData);
    
    % Save reduced-size solution to an array
    sol              = WFObj.sol;
    sol.site         = WFObj.Wp.site; % Save site info too
    sol_array(sol.k) = sol;
    
    scriptOptions = mergeStruct(scriptOptions,WFObj.scriptOptions);
    [ hFigs,~ ] = WFObs_s_animations( WFObj.Wp,sol_array,WFObj.sys,LESData,scriptOptions,WFObj.strucObs,hFigs );
end