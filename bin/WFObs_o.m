function [ sol,Wp,strucObs ] = WFObs_o(strucObs,Wp,sys,B1,B2,bc,input,measured,sol,k,it,options)
%[ sol, strucObs ] = WFObs_o(strucObs,Wp,sys,B1,B2,bc,input,measured,sol,k,it,startUniform)
%   This script selects which function to call for state reconstruction:
%   the Extended Kalman Filter (ExKF), Ensemble Kalman Filter (EnKF), or
%   no filter at all (sim).
%
%    Inputs:
%     *strucObs        structure containing observer information for time k-1
%     *Wp              structure containing meshing information
%     *sys             structure containing system matrices
%     *B1,B2,bc        system matrices related to the boundary conditions
%     *input           structure with turbine control settings (yaw and a)
%     *measured        structure containing (SOWFA) measurement information
%     *sol             flow fields and system state vector for time k-1
%     *k,it            current sample and iteration number, respectively
%     *options         structure containing model/script option information
%
%    Outputs:
%     *sol             flow fields and system state vector for time k
%     *strucObs        structure containing observer information for time k

switch lower(strucObs.filtertype)
    case 'exkf'
        [sol,strucObs] = WFObs_o_exkf(strucObs,Wp,sys,B1,B2,bc,input,measured,sol,k,it,options);
    case 'enkf'
        [sol,Wp,strucObs] = WFObs_o_enkf(strucObs,Wp,sys,B1,B2,bc,input,measured,sol,k,it,options);
    case 'ukf'
        [sol,Wp,strucObs] = WFObs_o_ukf(strucObs,Wp,sys,B1,B2,bc,input,measured,sol,k,it,options);
    case 'sim'
        [sys,power,~,~,~] = Make_Ax_b(Wp,sys,sol,input,B1,B2,bc,k,options); % Create system matrices
        [sol,~]           = Computesol(sys,input,sol,k,it,options);         % Compute solution by standard WFSim function
    otherwise
        error('not a valid filter specified.');
end;
end