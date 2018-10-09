function [measOptions] = sensorSet_flow_downstream(Wp)
measOptions.P.enabled      = false;  % Boolean for using power measurements
% measOptions.P.turbIds      = [1 2]; % Turbine ids from which measurements are taken
% measOptions.P.noiseStd     = 2e4;   % Standard deviation noise added to measurement
% measOptions.P.measStd      = 2e4;   % Standard deviation assumed by KF

delta_x = +2*Wp.turbine.Drotor;
nPoints = 4; % number of point measurements span-wise per turbine

measOptions.U.enabled    = true; % Boolean for using flow measurements
measOptions.U.locations = [];
for i = 1:length(Wp.turbine.Crx)
    for y = Wp.turbine.Cry(i)+Wp.turbine.Drotor*linspace(-0.5,0.5,5)
        measOptions.U.locations = [measOptions.U.locations; Wp.turbine.Crx(i)+delta_x, y];
    end
end
measOptions.U.noiseStd   = 1e-1; % Standard deviation noise added to measurement
measOptions.U.measStd    = 1e-1; % Standard deviation assumed by KF
end