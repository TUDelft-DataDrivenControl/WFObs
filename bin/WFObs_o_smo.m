function [Wp,sol_out,strucObs] = WFObs_o_smo(strucObs,Wp,sys_in,sol_in,options)
% WFOBS_O_SMO  Sliding mode observer for recursive state estimation
%
%   SUMMARY
%    This code performs state estimation using the Sliding Mode Observer
%    (SMO) algorithm. It uses high-fidelity measurements
%    (sol.measuredData) to improve the flow estimation compared to
%    open-loop simulations with WFSim.
%
%   RELEVANT INPUT/OUTPUT VARIABLES
%      see 'WFObs_o.m' for the complete list.
%   

% Import measurement data variable
measuredData = sol_in.measuredData;

if sol_in.k == 1
    % Setup covariance and system output matrices
    if options.exportPressures
        strucObs.Htt   = sparse(eye(strucObs.size_state));
    else
        strucObs.Htt   = sparse(eye(strucObs.size_output));
    end;
	strucObs.Htt   = strucObs.Htt(strucObs.obs_array,:);
    strucObs.Cinv = sparse(pinv(full(strucObs.Htt)));
end;

% Sliding mode observer forecast
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

% Analysis update
sol_out    = sol_in; % Copy previous solution before updating x
y          = measuredData.sol(strucObs.obs_array);
e          = y - solf.x(strucObs.obs_array);
Kgain      = abs(Fk*strucObs.Cinv*e + strucObs.alpha*strucObs.Cinv*y);
Kgain      = Kgain*Kgain'*strucObs.Cinv;
sol_out.x  = solf.x + Kgain*sign(e);

% Export new solution from estimation
[sol_out,~]  = MapSolution(Wp,sol_out,Inf,options); % Map solution to flowfields
[~,sol_out]  = Actuator(Wp,sol_out,options);        % Recalculate power after analysis update
end