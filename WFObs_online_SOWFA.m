% Description:
% This example uses zeroMQ to communicate with SOWFA (https://github.com/paulf81/SOWFA). More info on the SuperController in SOWFA
%  can be found in the corresponding pull request https://github.com/NREL/SOWFA/pull/22 .

%
% MATLAB can use zeroMQ, but it is not necessarily so straight-forward. The easiest solution found was
% using "jeroMQ", which can be downloaded from https://github.com/zeromq/jeromq. After installation,
% Update the path below and you should be all set.
%
% For more information, check out:
% https://mathworks.com/matlabcentral/answers/269061-how-do-i-integrate-zeromq-library-with-matlab-i-want-my-matlab-program-to-be-a-subscriber-of-zeromq
% 
% Note: to install jeroMQ, you need to have 'maven' installed. When using Maven to install jeroMQ,
% you may run into an error about the unit testing. If so, disable them and run again using
% 'mvn install -DskipTests'
%
% Recommended Java JDK version: 1.8.0_171 (tested by excluding unit tests)
%
%

% Setup zeroMQ server
zmqServer = zeromqObj('/home/bmdoekemeijer/OpenFOAM/zeroMQ/jeromq-0.4.4-SNAPSHOT.jar',5553,300,true);

% Setup WFSim
addpath('WFSim/layoutDefinitions') % Folder with predefined wind farm layouts
Wp = layoutSet_sowfa_9turb_apc_alm_turbl(); % Choose which scenario to simulate. See 'layoutDefinitions' folder for the full list.
addpath('WFSim/solverDefinitions'); % Folder with model options, solver settings, etc.
modelOptions = solverSet_default(Wp); % Choose model solver options. Default for EnKF/UKF: solverSet_minimal. For ExKF: solverSet_linearmatrices.

% Setup EnKF
addpath('filterDefinitions') % Folder with predefined KF settings
strucObs = filterSet_WES2018(); % Observer/KF settings
addpath('sensorDefinitions') % Folder with sensor settings
measOptions = sensorSet_power_only(Wp); % Measurement options

% Initialize WFObs object
addpath('bin','bin_supplementary/online_zmq')
WFObj = WFObs_obj( Wp,modelOptions,strucObs ); % Initialize WFObj object

% Initial control settings
nTurbs = length(Wp.turbine.Crx);
yawAngleArrayOut   = 270.0*ones(1,nTurbs);
pitchAngleArrayOut = 0.0*ones(1,nTurbs);


disp(['Entering wind farm controller loop...']);
firstRun = true;
measurementVector = [];
while 1
    % Receive information from SOWFA
    dataReceived = zmqServer.receive();
    currentTime  = dataReceived(1,1);
    measurements = dataReceived(1,2:end);
    
    if firstRun
        initialTimeSOWFA = floor(currentTime/100)*100; % Round down by 100s
        dataSend = setupZmqSignal(yawAngleArrayOut,pitchAngleArrayOut); % Initial string
        firstRun = false;
    end
    
    % Synchronize time with WFObj
    if ( (currentTime-initialTimeSOWFA) >= (WFObj.model.sol.time + WFObj.model.Wp.sim.h) )
        disp([datestr(rem(now,1)) '__    Doing estimation cycle.']);
        measuredData = struct();
        for i = 1:length(measOptions.P.turbIds)
            turbId = measOptions.P.turbIds(i);
            measuredData(i).idx   = turbId; % Turbine number
            measuredData(i).type  = 'P'; % Power measurement
            measuredData(i).value = measurements((turbId-1)*3+1);
            measuredData(i).std   = 1e4; % Standard deviation in W
        end
        turbData.CT_prime = 2.0*ones(nTurbs,1); % Greedy
        turbData.phi = 270.0-yawAngleArrayOut';
        
        sol = WFObj.timestepping(turbData,measuredData);
    
        % Optimization
        disp([datestr(rem(now,1)) '__    Doing optimization cycle.']);
        yawAngleArrayOut   = 270.0*ones(1,nTurbs);
        yawAngleArrayOut(1) = 270+15*sin(2*pi*(1/50)*currentTime);
        yawAngleArrayOut(2) = 270+15*sin(2*pi*(1/50)*(currentTime+10));
        yawAngleArrayOut(3) = 270+15*sin(2*pi*(1/50)*(currentTime+20));
        yawAngleArrayOut(4) = 270+15*sin(2*pi*(1/50)*(currentTime+30));
        yawAngleArrayOut(5) = 270+15*sin(2*pi*(1/50)*(currentTime+40));
        
        % Update message string
        disp([datestr(rem(now,1)) '__    Synthesizing message string.']);
        dataSend = setupZmqSignal(yawAngleArrayOut,pitchAngleArrayOut);
    end
    
    % Send a message (control action) back to SOWFA
    zmqServer.send(dataSend);
end
% Close connection
zmqServer.disconnect()

function [dataOut] = setupZmqSignal(yawAngles,pitchAngles)
	dataOut = [];
    for i = 1:length(yawAngles)
        dataOut = [dataOut yawAngles(i) pitchAngles(i)];
    end
end