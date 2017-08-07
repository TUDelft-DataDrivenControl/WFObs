function [ Wp,sol,strucObs ] = WFObs_o(strucObs,Wp,sys,sol,options)
%[ sol, strucObs ] = WFObs_o(strucObs,Wp,sys,B1,B2,bc,input,measured,sol,k,it,startUniform)
%   This script selects which function to call for state reconstruction:
%   the Extended Kalman Filter (ExKF), Ensemble Kalman Filter (EnKF), or
%   no filter at all (sim).
%
%    Inputs:
%     *strucObs        structure containing observer information for time k-1
%     *Wp              structure containing meshing information
%     *sys             structure containing system matrices
%     *measuredData    structure containing (SOWFA) measurement information
%     *sol             flow fields and system state vector for time k-1
%     *it              current iteration number
%     *options         structure containing model/script option information
%
%    Outputs:
%     *Wp              structure containing meshing information
%     *sol             flow fields and system state vector for time k
%     *strucObs        structure containing observer information for time k

switch lower(strucObs.filtertype)
    case 'exkf'
        [Wp,sol,strucObs] = WFObs_o_exkf(strucObs,Wp,sys,sol,options);
    case 'enkf'
        [Wp,sol,strucObs] = WFObs_o_enkf(strucObs,Wp,sys,sol,options);
    case 'ukf'
        [Wp,sol,strucObs] = WFObs_o_ukf( strucObs,Wp,sys,sol,options);
    case 'sim'
        sol.k    = sol.k - 1; % Necessary since WFSim_timestepping(...) already includes time propagation
        [sol,~]  = WFSim_timestepping(sol,sys,Wp,options);
    otherwise
        error('not a valid filter specified.');
end;

end