function [ sol ] = WFObs_s_reporting( timerCPU,Wp,sol,strucObs,scriptOptions )
% WFOBS_S_REPORTING  Determine performance scores and print progress
%
%   SUMMARY
%    This code determines several performance scores for the estimation,
%    and prints an important subset of them in the command line window.
%
%   RELEVANT INPUT/OUTPUT VARIABLES
%     - timerCPU: time ('tic') taken at the start of the iteration (s)
%
%     - Wp: this struct contains all the simulation settings related to the
%           wind farm, the turbine inputs, the atmospheric properties, etc.
%
%     - sol: this struct contains the system states at a certain timestep.
%         sol.score.maxErroru: maximum error (m/s) in long. flow direction
%         sol.score.maxErrorv: maximum error (m/s) in lat.  flow direction
%         sol.score.maxError:  maximum error (m/s) in either flow direction
%         sol.score.maxErroruloc: location of max. error in long. flow direction
%         sol.score.maxErrorvloc: location of max. error in lat.  flow direction
%         sol.score.RMSE_flow: RMSE (m/s) of the farm-wide flow velocity
%         sol.score.RMSE_meas: RMSE (m/s) with the measured flow velocities
%         sol.score.RMSE_cline: RMSE (m/s) of the centerline flow velocity
%         sol.score.VAF_cline: Variance accounted for (%) of the centerline flow velocity
%         sol.score.CPUtime: computational cost (s) for estimation update
%
%     - strucObs: this struct contains all the observer settings and
%       (temporary) files used for updates, such as covariance matrices,
%       ensemble/sigma point sets, measurement noise, etc.

%     - scriptOptions: this struct contains all simulation settings, not
%       related to the wind farm itself (solution methodology, outputs, etc.)
%

% Determine the maximum error values and locations
i_corr = find(sol.measuredData.uq ~= 0);
[sol.score.maxErroru,sol.score.maxerroruloc] = max(abs(sol.u(i_corr)-sol.measuredData.uq(i_corr)));
[sol.score.maxErrorv,sol.score.maxerrorvloc] = max(abs(sol.v(i_corr)-sol.measuredData.vq(i_corr)));
[sol.score.maxError ] = max( [sol.score.maxErrorv , sol.score.maxErroru] );

% Determine centerline flow speeds
yend_left  		 = min([Wp.mesh.yline{:}]);
yend_right 		 = max([Wp.mesh.yline{:}]);
centerline_WFSim = mean(sol.u(:,yend_left-1:yend_right),2);
centerline_SOWFA = mean(sol.measuredData.uq(:,yend_left-1:yend_right),2);

% Calculate several scores
sol.score.RMSE_cline = sqrt(mean((centerline_WFSim-centerline_SOWFA).^2));	 % Centerline RMSE (m/s)
sol.score.RMSE_flow  = rms([sol.v(i_corr)-sol.measuredData.vq(i_corr);sol.u(i_corr)-sol.measuredData.uq(i_corr)]);  % Total flow RMSE (m/s)
sol.score.RMSE_meas  = sqrt(mean((sol.measuredData.sol (strucObs.obs_array)-sol.x(strucObs.obs_array)).^2)); % Measured flow RMSE (m/s)
sol.score.VAF_cline  = vaf(centerline_SOWFA,centerline_WFSim); % Variance accounted for (VAF) in (%)
sol.score.CPUtime    = toc(timerCPU); % Computational cost

% Plot progress
if scriptOptions.printProgress
    disp([datestr(rem(now,1)) ' __  t(' num2str(sol.k,['%0' num2str(scriptOptions.tlen) 'd']) ') = ' num2str(sol.time,['%0' num2str(scriptOptions.klen) 'd']) ' s __ Centerline VAF: ' num2str(sol.score.VAF_cline,'%10.2f\n') '%, Centerline RMSE: ' num2str(sol.score.RMSE_cline,'%10.2f\n') ' m/s, Measurement RMSE: ' num2str(sol.score.RMSE_meas,'%10.2f\n') ' m/s.']);
    disp([datestr(rem(now,1)) ' __  t(' num2str(sol.k,['%0' num2str(scriptOptions.tlen) 'd']) ') = ' num2str(sol.time,['%0' num2str(scriptOptions.klen) 'd']) ' s __ Flow RMSE: ' num2str(sol.score.RMSE_flow,'%10.2f\n'), ' m/s, u_Inf: ' num2str(Wp.site.u_Inf,'%10.2f\n') ', v_Inf: ' num2str(Wp.site.v_Inf,'%10.2f\n') ', it. time: ' num2str(sol.score.CPUtime,'%10.2f\n') ' s.']);
    if strcmp(lower(strucObs.filtertype),'enkf') | strcmp(lower(strucObs.filtertype),'ukf')
        if strucObs.tune.est
            for iT = 1:length(strucObs.tune.vars)
                disp([datestr(rem(now,1)) ' __  t(' num2str(sol.k,['%0' num2str(scriptOptions.tlen) 'd']) ') = ' num2str(sol.time,['%0' num2str(scriptOptions.klen) 'd']) ' s __ ' strucObs.tune.vars{iT} ' estimated as ' num2str(Wp.(strucObs.tune.subStruct{iT}).(strucObs.tune.structVar{iT}),'%10.2f\n') '.']);
            end
        end
    end
    disp(' ')
end
end