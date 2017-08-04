function [ measured ] = WFObs_s_loadmeasurements( strucObs, time )
%[ measured ] = WFObs_s_loadmeasurements(sourcepath, datanroffset,time,noise_obs)
%   This script loads SOWFA data files and formats them to the correct
%   format for usage with WFObs. Inputs and outputs are:
%
%    Inputs:
%     *sourcepath      folder where SOWFA files are located
%     *datanroffset    offset for file numbers
%     *time            current time/iteration number
%     *noise_obs       artificial noise added to the velocity measurements
%
%    Outputs:
%     *measured.power  output power of turbines, undisturbed
%     *measured.sol    velocities in format of 'x' in 'Ax=b', disturbed
%     *measured.solq   velocities in format of 'x' in 'Ax=b', undisturbed
%     *measured.u      longitudinal velocities, disturbed
%     *measured.uq     longitudinal velocities, undisturbed
%     *measured.v      lateral velocities, disturbed
%     *measured.vq     lateral velocities, undisturbed

% Import variables
sourcepath   = strucObs.measurementsPath;
datanroffset = strucObs.measurementsOffset;

measured                        = load([sourcepath '/' num2str(datanroffset+time) '.mat']);
measured.uq(isnan(measured.uq)) = 0; % nullify any out-of-range measurements
measured.vq(isnan(measured.vq)) = 0; % nullify any out-of-range measurements

% Add noise to measurements that are to be fed into the observer
measured.u      = measured.uq + strucObs.noise_obs*randn(size(measured.uq)); % Disturbed
measured.v      = measured.vq + strucObs.noise_obs*randn(size(measured.vq)); % Disturbed

% Sort them into a vector of correct size (for Ax=b)
measured.solq  = [vec(measured.uq(3:end-1,2:end-1)'); vec(measured.vq(2:end-1,3:end-1)')];
measured.sol   = [vec(measured.u(3:end-1,2:end-1)') ; vec(measured.v(2:end-1,3:end-1)') ];
end