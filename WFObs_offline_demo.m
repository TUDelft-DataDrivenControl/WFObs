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

%% Settings for offline WFObs simulations with a LES database
configName          = 'TORQUE_axi_2turb_alm_turb';

postProcOptions.Animate       = 20;  % Animation frequency
postProcOptions.plotContour   = 1;  % Show flow fields
postProcOptions.plotPower     = 1;  % Plot true and predicted power capture vs. time
postProcOptions.powerForecast = 10; % Plot power forecast (0 = disabled, x = number of steps) (only if plotPower = 1)
postProcOptions.plotError     = 1;  % plot RMS and maximum error vs. time
postProcOptions.savePlots     = 0;  % Save all plots in external files at each time step
postProcOptions.savePath      = ['results/tmp']; % Destination folder of saved files

measPw.enabled      = false;  % Boolean for using power measurements
measPw.turbIds      = [1 2]; % Turbine ids from which measurements are taken
measPw.noiseStd     = 2e4;   % Standard deviation noise added to measurement
measPw.measStd      = 2e4;   % Standard deviation assumed by KF

measFlow.enabled    = true; % Boolean for using flow measurements
measFlow.sensorFile = 'setup_sensors/sensors_2turb_alm.mat';
measFlow.noiseStd   = 1e-1; % Standard deviation noise added to measurement
measFlow.measStd    = 1e-1; % Standard deviation assumed by KF

%% Initialize object
addpath('bin'); % Add the main 'bin' folder
addpath('offline_vis_tools'); % Add the postProcessing folder
WFObj=WFObs_obj(configName); % See './configurations' for options

%% Preload LES measurement data and setup sensors
LESData = load(WFObj.Wp.sim.measurementFile);

%% Execute the WFObs core code
hFigs = [];

% Load sensor file (backwards compatibility, legacy format)
if measFlow.enabled
    sensorInfo = load(measFlow.sensorFile);
end

while WFObj.sol.k < WFObj.Wp.sim.NN
    
    % Load and format measurements from preloaded database
    measuredData = [];
    if measPw.enabled % Setup power measurements
        for i = 1:length(measPw.turbIds)
            measuredData(i).idx   = measPw.turbIds(i); % Turbine number
            measuredData(i).type  = 'P'; % Power measurement
            measuredData(i).value = LESData.turbData.power(WFObj.sol.k+1,...
                                    measPw.turbIds(i)) + measPw.noiseStd*randn(); % With artificial noise
            measuredData(i).std   = measPw.measStd; % Standard deviation in W
        end
    end
    
    if measFlow.enabled % Setup flow measurements
        for jT = [1 2] % jT=1 for 'u', jT=2 for 'v' measurements
            iOffset = length(measuredData);
            for i = 1:sensorInfo.sensors{jT}.N
                measuredData(iOffset+i).idx = sensorInfo.sensors{jT}.loc(i,:); % Sensor location
                measuredData(iOffset+i).std = measFlow.measStd; % Standard deviation in W
                if jT == 1
                    measuredData(iOffset+i).type  = 'u'; % Long. flow measurement
                    measuredData(iOffset+i).value = LESData.u(WFObj.sol.k+1, ...
                        sensorInfo.sensors{1}.grid(i,1),sensorInfo.sensors{1}.grid(i,2)) + ...
                        measFlow.noiseStd*randn(); % Plus artificial noise
                elseif jT == 2
                    measuredData(iOffset+i).type  = 'v'; % Lat. flow measurement
                    measuredData(iOffset+i).value = LESData.v(WFObj.sol.k+1, ...
                        sensorInfo.sensors{2}.grid(i,1),sensorInfo.sensors{2}.grid(i,2)) + ...
                        measFlow.noiseStd*randn(); % Plus artificial noise
                end
            end
        end
        clear iOffset
    end
    
    % Perform estimation
    inputData = WFObj.Wp.turbine.input(WFObj.sol.k+1);
    WFObj.timestepping(inputData,measuredData);
    
    % Save reduced-size solution to an array
    sol = WFObj.sol;
    flowError = [sol.v(:)-vec(LESData.v(sol.k,:,:));sol.u(:)-vec(LESData.u(sol.k,:,:))];
    sol.score.RMSE_flow     = rms(flowError); % Total flow RMSE (m/s)
    sol.score.maxError_flow = max(abs(flowError)); % Maximum error
    sol.site = WFObj.Wp.site; % Save site info too (contains model parameters that may be estimated)
    sol_array(sol.k) = sol;
    
    postProcOptions = mergeStruct(WFObj.scriptOptions,postProcOptions);
    [ hFigs,postProcOptions ] = WFObs_p_animations( WFObj.Wp,sol_array,WFObj.sys,LESData,measuredData,postProcOptions,WFObj.strucObs,hFigs );
end