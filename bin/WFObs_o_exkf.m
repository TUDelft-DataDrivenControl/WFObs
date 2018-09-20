function [Wp,sol_out,strucObs] = WFObs_o_exkf(strucObs,Wp,sys_in,sol_in,options)
% WFOBS_O_EXKF  Extended KF algorithm for recursive state estimation
%
%   SUMMARY
%    This code performs state estimation using the Extended Kalman filter
%    (ExKF) algorithm. It uses high-fidelity measurements
%    (sol.measuredData) to improve the flow estimation compared to
%    open-loop simulations with WFSim.
%
%   RELEVANT INPUT/OUTPUT VARIABLES
%      see 'WFObs_o.m' for the complete list.
%   

% Import measurement data variable
measuredData = sol_in.measuredData;

if sol_in.k == 1
    if strucObs.measPw
        error(['ExKF currently does not support power measurements. '...
               'Please change to the EnKF/UKF or use flow measurements.']);
    end
    
    if strucObs.pe.enabled
        error(['ExKF currently does not support parameter estimation. '...
               'Please change to the EnKF/UKF or use state estimation only.']);
    end
    
    % Setup covariance and system output matrices
    strucObs.Htt = sparse(eye(Wp.Nu+Wp.Nv+options.exportPressures*Wp.Np));
    strucObs.Htt = strucObs.Htt(strucObs.obs_array,:);    
    strucObs.Pk  = blkdiag(eye(Wp.Nu)*strucObs.se.P0.u,eye(Wp.Nv)*strucObs.se.P0.v);
    strucObs.Qk  = blkdiag(eye(Wp.Nu)*strucObs.se.Qk.u,eye(Wp.Nv)*strucObs.se.Qk.v);
    strucObs.Rk  = blkdiag(eye(Wp.Nu)*strucObs.se.Rk.u,eye(Wp.Nv)*strucObs.se.Rk.v);
    strucObs.Rk  = strucObs.Rk(strucObs.obs_array,strucObs.obs_array);
    if options.exportPressures
        strucObs.Pk = blkdiag(strucObs.Pk,eye(Wp.Np)*strucObs.P0.p);
        strucObs.Qk = blkdiag(strucObs.Qk,eye(Wp.Np)*strucObs.Qk.p);
    end
end;

% ExKF forecast update
soltemp   = sol_in;
soltemp.k = soltemp.k - 1;
[solf,sysf]             = WFSim_timestepping( soltemp, sys_in, Wp, options );       % Forward propagation
Fk(sysf.pRCM,sysf.pRCM) = sysf.A(sysf.pRCM,sysf.pRCM)\sysf.Al(sysf.pRCM,sysf.pRCM); % Linearized A-matrix at time k
Bk(sysf.pRCM,:)         = sysf.A(sysf.pRCM,sysf.pRCM)\sysf.Bl(sysf.pRCM,:);         % Linearized B-matrix at time k

% Neglect pressure terms
if ~options.exportPressures 
    Fk     = Fk(1:strucObs.size_output,1:strucObs.size_output);
    Bk     = Bk(1:strucObs.size_output,:);
    solf.x = solf.x(1:strucObs.size_output);
end;

Pf = Fk*strucObs.Pk*Fk' + strucObs.Qk;  % Covariance matrix P for x(k) knowing y(k-1)

% ExKF analysis update
sol_out     = sol_in; % Copy previous solution before updating x
Kgain       = Pf(:,strucObs.obs_array)*pinv(Pf(strucObs.obs_array,...
                   strucObs.obs_array)+strucObs.Rk); % Kalman gain           
sol_out.x   = solf.x + Kgain*(measuredData.sol(strucObs.obs_array)...
                 -solf.x(strucObs.obs_array)); % Optimally predicted state vector
strucObs.Pk = (eye(size(Pf))-Kgain*strucObs.Htt)*Pf;  % State covariance matrix

% Export new solution from estimation
[sol_out,~]  = MapSolution(Wp,sol_out,Inf,options); % Map solution to flowfields
[~,sol_out]  = Actuator(Wp,sol_out,options);        % Recalculate power after analysis update
end