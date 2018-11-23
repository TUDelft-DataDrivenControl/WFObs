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
addpath('../bin','../bin_supplementary/online_zmq')
zmqServer = zeromqObj('/home/bmdoekemeijer/OpenFOAM/zeroMQ/jeromq-0.4.4-SNAPSHOT.jar',5553,300,true);

% Setup WFSim
Wp = struct('description','9 NREL 5MW turbines SSC demo case');
Wp.sim = struct(...
    'h',1.0,... % timestep (s)
    'startUniform',true ... % Start from a uniform flow field (T) or from a fully developed waked flow field (F).
    );
Wp.turbine = struct(...
    'Crx',[468   1100  1732   468 1100 1732  468    1100   1732],... % X-coordinates of turbines (m)
    'Cry',[320.8 320.8 320.8  700 700  700   1079.2 1079.2 1079.2],... % Y-coordinates of turbines (m)
    'Drotor',126,... % Rotor diameter (m), note that WFSim only supports a uniform Drotor for now
    'powerscale',0.97,... % Turbine power scaling
    'forcescale',2.0 ... % Turbine force scaling
    );
Wp.site = struct(...
    'u_Inf',12.0214,... % Initial long. wind speed in m/s
    'v_Inf',0.0,... % Initial lat. wind speed in m/s
    'p_init',0.0,... % Initial values for pressure terms (Pa)
    'turbul',true,... % Use mixing length turbulence model (true/false)
    'turbModel','WFSim3',...  % Turbulence model of choice
    'lmu',4.0,... %1.20,... % Mixing length in x-direction (m)
    'm',7,... % Turbulence model gridding property
    'n',1,... % Turbulence model gridding property
    'mu',0.0,... % Dynamic flow viscosity
    'Rho',1.20 ... % Air density
    );
Wp.mesh = struct(...
    'gridType','lin',... % Grid type ('lin' the only supported one currently)
    'Lx',2400,... % Domain length in x-direction
    'Ly',1400,... % Domain length in y-direction
    'Nx',80,... % Number of cells in x-direction
    'Ny',42 ... % Number of cells in y-direction
    );
    
addpath('../WFSim/solverDefinitions'); % Folder with model options, solver settings, etc.
modelOptions = solverSet_minimal(Wp); % Choose model solver options. Default for EnKF/UKF: solverSet_minimal. For ExKF: solverSet_linearmatrices.

% Setup EnKF
addpath('../filterDefinitions') % Folder with predefined KF settings
strucObs = filterSet_WES2018(); % Observer/KF settings
addpath('../sensorDefinitions') % Folder with sensor settings
measOptions = sensorSet_power_only(Wp); % Measurement options

% Initialize WFObs object
WFObj = WFObs_obj( Wp,modelOptions,strucObs );

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