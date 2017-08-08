%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 'WFObs_s_reporting.m'
%  This script displays the progress and some primary results in the
%  command line. It includes:
%    * The RMS error in flow estimations (m/s)
%    * The maximum error in flow field estimations (m/s)
%    * The iteration time (s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ sol ] = WFObs_s_reporting( timerCPU,Wp,sol,strucObs,scriptOptions )
% WFOBS_S_REPORTING  Display progress and performance of the observer code
    i_corr                   = find(sol.measuredData.uq ~= 0); 
    [sol.score.maxErroru,sol.score.maxerroruloc] = max(abs(sol.u(i_corr)-sol.measuredData.uq(i_corr)));
    [sol.score.maxErrorv,sol.score.maxerrorvloc] = max(abs(sol.v(i_corr)-sol.measuredData.vq(i_corr)));
    [sol.score.maxError ]    = max( [sol.score.maxErrorv , sol.score.maxErroru] );

	% Determine centerline performance
    yend_left  		 = min([Wp.mesh.yline{:}]);
	yend_right 		 = max([Wp.mesh.yline{:}]);
	centerline_WFSim = mean(sol.u(:,yend_left-1:yend_right),2);
	centerline_SOWFA = mean(sol.measuredData.uq(:,yend_left-1:yend_right),2);
			
	% Calculate scores
	sol.score.RMSE_cline = sqrt(mean((centerline_WFSim-centerline_SOWFA).^2));	 % Centerline RMSE (m/s)
    sol.score.RMSE_flow  = rms([sol.v(i_corr)-sol.measuredData.vq(i_corr);sol.u(i_corr)-sol.measuredData.uq(i_corr)]);  % Total flow RMSE (m/s)
	sol.score.RMSE_meas  = sqrt(mean((sol.measuredData.sol (strucObs.obs_array)-sol.x(strucObs.obs_array)).^2)); % Measured flow RMSE (m/s)
	sol.score.VAF_cline  = vaf(centerline_SOWFA,centerline_WFSim); % Variance accounted for (VAF) in (%)
    sol.score.CPUtime    = toc(timerCPU); % Computational cost

    if scriptOptions.printProgress
        disp([datestr(rem(now,1)) ' __  t(' num2str(sol.k,['%0' num2str(scriptOptions.tlen) 'd']) ') = ' num2str(sol.time,['%0' num2str(scriptOptions.klen) 'd']) ' s __ Centerline VAF: ' num2str(sol.score.VAF_cline,'%10.2f\n') '%, Centerline RMSE: ' num2str(sol.score.RMSE_cline,'%10.2f\n') ' m/s, Measurement RMSE: ' num2str(sol.score.RMSE_meas,'%10.2f\n') ' m/s.']);
		disp([datestr(rem(now,1)) ' __  t(' num2str(sol.k,['%0' num2str(scriptOptions.tlen) 'd']) ') = ' num2str(sol.time,['%0' num2str(scriptOptions.klen) 'd']) ' s __ Flow RMSE: ' num2str(sol.score.RMSE_flow,'%10.2f\n'), ' m/s, u_Inf: ' num2str(Wp.site.u_Inf,'%10.2f\n') ', v_Inf: ' num2str(Wp.site.v_Inf,'%10.2f\n') ', it. time: ' num2str(sol.score.CPUtime,'%10.2f\n') ' s.']);
        if strcmp(lower(strucObs.filtertype),'enkf') | strcmp(lower(strucObs.filtertype),'ukf')
            for iT = 1:length(strucObs.tune.vars)
                disp([datestr(rem(now,1)) ' __  t(' num2str(sol.k,['%0' num2str(scriptOptions.tlen) 'd']) ') = ' num2str(sol.time,['%0' num2str(scriptOptions.klen) 'd']) ' s __ ' strucObs.tune.vars{iT} ' estimated as ' num2str(Wp.(strucObs.tune.subStruct{iT}).(strucObs.tune.structVar{iT}),'%10.2f\n') '.']);
            end
        end
        disp(' ')
    end
end