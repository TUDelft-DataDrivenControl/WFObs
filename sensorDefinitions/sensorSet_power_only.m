function [measOptions] = sensorSet_power_only(Wp)
measOptions.P.enabled      = true;  % Boolean for using power measurements
measOptions.P.turbIds      = 1:1:length(Wp.turbine.Crx); % Turbine ids from which measurements are taken
measOptions.P.noiseStd     = 2e4;   % Standard deviation noise added to measurement
measOptions.P.measStd      = 2e4;   % Standard deviation assumed by KF

measOptions.U.enabled    = false; % Boolean for using flow measurements
end