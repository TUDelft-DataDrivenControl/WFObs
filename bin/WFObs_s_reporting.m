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

    sol.score.RMSE    = rms([sol.v(i_corr)-sol.measuredData.vq(i_corr);sol.u(i_corr)-sol.measuredData.uq(i_corr)]);
    sol.score.CPUtime = toc(timerCPU);

    if scriptOptions.printProgress
        disp([datestr(rem(now,1)) ' __  t(' num2str(sol.k,['%0' num2str(scriptOptions.tlen) 'd']) ') = ' num2str(sol.time,['%0' num2str(scriptOptions.klen) 'd']) ' s __ max. error: ' num2str(sol.score.maxError,'%10.2f\n'),' m/s, RMSE: ' num2str(sol.score.RMSE,'%10.2f\n'), ' m/s, it. time: ' num2str(sol.score.CPUtime,'%10.2f\n') ' s.']);
        disp([datestr(rem(now,1)) ' __  t(' num2str(sol.k,['%0' num2str(scriptOptions.tlen) 'd']) ') = ' num2str(sol.time,['%0' num2str(scriptOptions.klen) 'd']) ' s __ u_Inf estimated as ' num2str(Wp.site.u_Inf,'%10.2f\n') '.']);
        disp([datestr(rem(now,1)) ' __  t(' num2str(sol.k,['%0' num2str(scriptOptions.tlen) 'd']) ') = ' num2str(sol.time,['%0' num2str(scriptOptions.klen) 'd']) ' s __ v_Inf estimated as ' num2str(Wp.site.v_Inf,'%10.2f\n') '.']);
        if strcmp(lower(strucObs.filtertype),'enkf') | strcmp(lower(strucObs.filtertype),'ukf')
            for iT = 1:length(strucObs.tune.vars)
                disp([datestr(rem(now,1)) ' __  t(' num2str(sol.k,['%0' num2str(scriptOptions.tlen) 'd']) ') = ' num2str(timeindex,['%0' num2str(scriptOptions.klen) 'd']) ' s __ ' strucObs.tune.vars{iT} ' estimated as ' num2str(Wp.(strucObs.tune.subStruct{iT}).(strucObs.tune.structVar{iT}),'%10.2f\n') '.']);
            end
        end
        disp(' ')
    end
end