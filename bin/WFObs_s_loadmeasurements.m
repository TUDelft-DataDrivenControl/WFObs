function [ measuredData ] = WFObs_s_loadmeasurements( LESData, k )
% WFOBS_S_LOADMEASUREMENTS  Loads measurement data for estimation
%
%   SUMMARY
%    This script loads SOWFA data files and formats them to the correct
%    format for usage with WFObs. Inputs and outputs are:
%
%   RELEVANT INPUT/OUTPUT VARIABLES
%     - LESData: this struct contains all the flow fields and turbine data
%                from the LES data.
%
%     - k: current time instant (belonging to index).
%
%     - measured: measurement data struct with the following entries
%        *measuredData.power  output power of turbines, undisturbed
%        *measuredData.sol    velocities in format of 'x' in 'Ax=b', disturbed
%        *measuredData.solq   velocities in format of 'x' in 'Ax=b', undisturbed
%        *measuredData.u      longitudinal velocities, disturbed
%        *measuredData.uq     longitudinal velocities, undisturbed
%        *measuredData.v      lateral velocities, disturbed
%        *measuredData.vq     lateral velocities, undisturbed
%

measuredData.uq    = squeeze(LESData.u(k,:,:));
measuredData.vq    = squeeze(LESData.v(k,:,:));
measuredData.u     = squeeze(LESData.ud(k,:,:));
measuredData.v     = squeeze(LESData.vd(k,:,:));
measuredData.solq  = [vec(measuredData.uq(3:end-1,2:end-1)'); vec(measuredData.vq(2:end-1,3:end-1)')];
measuredData.sol   = [vec(measuredData.u(3:end-1,2:end-1)') ; vec(measuredData.v(2:end-1,3:end-1)') ];
measuredData.power(:,1) = LESData.turbData.power(k,:);
end