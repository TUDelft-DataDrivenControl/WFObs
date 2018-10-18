function [strucObs,model] = WFObs_o_exkf(strucObs,model)
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

% Setup variables
Wp = model.Wp;
sys_in = model.sys;
sol_in = model.sol;
options = model.modelOptions;
turbInput = sol_in.turbInput;

% Import measurement data variable
measuredData = sol_in.measuredData;

if sol_in.k == 1
    if strucObs.pe.enabled
        error(['ExKF currently does not support parameter estimation. '...
            'Please change to the EnKF/UKF or use state estimation only.']);
    end
    
    % Construct measurement matrix (C-matrix) based on a 'nearest' approach
    obs_array = [];
    for i = 1:length(measuredData)
        addpath('../setup_sensors');
        if strcmp(sol_in.measuredData(i).type,'P')
            error('The ExKF currently does not support the use of power measurements.');
        elseif strcmp(sol_in.measuredData(i).type,'u')
            [~,gridcor(1,1)] = min(abs(measuredData(i).idx(1)-Wp.mesh.ldxx2(:,1)));
            [~,gridcor(1,2)] = min(abs(measuredData(i).idx(2)-Wp.mesh.ldyy(1,:)));
            [ obs_array ] = [obs_array; WFObs_s_sensors_grid2nr( gridcor, Wp , 'u')];
        elseif strcmp(sol_in.measuredData(i).type,'v')
            [~,gridcor(1,1)] = min(abs(measuredData(i).idx(1)-Wp.mesh.ldxx(:,1)));
            [~,gridcor(1,2)] = min(abs(measuredData(i).idx(2)-Wp.mesh.ldyy2(1,:)));
            [ obs_array ] = [obs_array; WFObs_s_sensors_grid2nr( gridcor, Wp , 'v')];
        else
            error('You specified an incompatible measurement. Please use types ''u'', ''v'', or ''P'' (capital-sensitive).');
        end
    end
    strucObs.Htt = sparse(eye(Wp.Nu+Wp.Nv+options.exportPressures*Wp.Np));
    strucObs.Htt = strucObs.Htt(obs_array,:);
    strucObs.obs_array = obs_array;
    
    % Setup covariance and system output matrices
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
[solf,sysf]             = WFSim_timestepping( soltemp, sys_in, Wp, turbInput, options ); % Forward propagation
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
sol_out.x   = solf.x + Kgain*([measuredData.value]'-solf.x(strucObs.obs_array)); % Optimally predicted state vector
strucObs.Pk = (eye(size(Pf))-Kgain*strucObs.Htt)*Pf;  % State covariance matrix

% Export new solution from estimation
[sol_out,~]  = MapSolution(Wp,sol_out,Inf,options); % Map solution to flowfields

% Recalculate power after analysis update
sol_out.turbInput.dCT_prime = zeros(Wp.turbine.N,1);
[~,sol_out]  = Actuator(Wp,sol_out,options); 

% Export variables
model.sol = sol_out;
end