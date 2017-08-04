% Add WFSim paths to MATLABs working environment
[WFObsPath, ~, ~] = fileparts(which('WFObs_addpaths.m')); % Get /bin/ path
addpath([WFObsPath '/bin']);                            % Add /bin/ path
addpath([WFObsPath '/configurations'])                  % Add /configurations/ path
addpath([WFObsPath '/setup_sensors/sensor_layouts'])    % Add /sensors/ path
run(    [WFObsPath '/WFSim/WFSim_addpaths.m'])          % Add /WFSim/ paths
clear WFObsPath