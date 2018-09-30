function [measPwOptions,measFlowOptions] = sensorSet_2turb_alm()
measPwOptions.enabled      = true;  % Boolean for using power measurements
measPwOptions.turbIds      = [1 2]; % Turbine ids from which measurements are taken
measPwOptions.noiseStd     = 2e4;   % Standard deviation noise added to measurement
measPwOptions.measStd      = 2e4;   % Standard deviation assumed by KF

measFlowOptions.enabled    = true; % Boolean for using flow measurements
measFlowOptions.locations  = [ 326.4794  333.3333;
                        326.4794  366.6667;
                        326.4794  400.0000;
                        326.4794  433.3333;
                        326.4794  466.6667;
                        979.4381  333.3333;
                        979.4381  366.6667;
                        979.4381  400.0000;
                        979.4381  433.3333;
                        979.4381  466.6667 ];
measFlowOptions.noiseStd   = 1e-1; % Standard deviation noise added to measurement
measFlowOptions.measStd    = 1e-1; % Standard deviation assumed by KF
end