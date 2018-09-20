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

try
    % Determine the maximum error values and locations
    [sol.score.maxErroru,sol.score.maxerroruloc] = max(abs(sol.u(:)-sol.measuredData.uq(:)));
    [sol.score.maxErrorv,sol.score.maxerrorvloc] = max(abs(sol.v(:)-sol.measuredData.vq(:)));
    [sol.score.maxError ] = max( [sol.score.maxErrorv , sol.score.maxErroru] );

    % Calculate several scores
    [ sol.cline,sol.measuredData.cline,sol.score.VAF_cline,sol.score.RMSE_cline ] = WFObs_s_cline( Wp,sol ); % Centerline
    sol.score.RMSE_flow  = rms([sol.v(:)-sol.measuredData.vq(:);sol.u(:)-sol.measuredData.uq(:)]);  % Total flow RMSE (m/s)
    sol.score.RMSE_meas  = sqrt(mean((sol.measuredData.sol (strucObs.obs_array)-sol.x(strucObs.obs_array)).^2)); % Measured flow RMSE (m/s)
catch
    sol.score.maxError   = 0;
    sol.score.VAF_cline  = 0;
    sol.score.RMSE_cline = 0;
    sol.score.RMSE_flow  = 0;
    sol.score.RMSE_meas  = 0;
end

sol.score.CPUtime = toc(timerCPU); % Computational cost

% Plot progress
try
    if scriptOptions.printProgress
        if length(sol.score.RMSE_cline) == 1
            disp([datestr(rem(now,1)) ' __  t(' num2str(sol.k,['%0' num2str(scriptOptions.tlen) 'd']) ') = ' num2str(sol.time,['%0' num2str(scriptOptions.klen) 'd']) ' s __ Centerline VAF: ' num2str(sol.score.VAF_cline,'%10.2f\n') '%, Centerline RMSE: ' num2str(sol.score.RMSE_cline,'%10.2f\n') ' m/s, Measurement RMSE: ' num2str(sol.score.RMSE_meas,'%10.2f\n') ' m/s.']);
            disp([datestr(rem(now,1)) ' __  t(' num2str(sol.k,['%0' num2str(scriptOptions.tlen) 'd']) ') = ' num2str(sol.time,['%0' num2str(scriptOptions.klen) 'd']) ' s __ Flow RMSE: ' num2str(sol.score.RMSE_flow,'%10.2f\n'), ' m/s, u_Inf: ' num2str(Wp.site.u_Inf,'%10.2f\n') ', v_Inf: ' num2str(Wp.site.v_Inf,'%10.2f\n') ', it. time: ' num2str(sol.score.CPUtime,'%10.2f\n') ' s.']);
        else
            disp([datestr(rem(now,1)) ' __  t(' num2str(sol.k,['%0' num2str(scriptOptions.tlen) 'd']) ') = ' num2str(sol.time,['%0' num2str(scriptOptions.klen) 'd']) ' s __ Row-avgd cline VAF: ' num2str(mean(sol.score.VAF_cline),'%10.2f\n') '%, Row-avgd cline RMSE: ' num2str(mean(sol.score.RMSE_cline),'%10.2f\n') ' m/s, Measurement RMSE: ' num2str(sol.score.RMSE_meas,'%10.2f\n') ' m/s.']);
            disp([datestr(rem(now,1)) ' __  t(' num2str(sol.k,['%0' num2str(scriptOptions.tlen) 'd']) ') = ' num2str(sol.time,['%0' num2str(scriptOptions.klen) 'd']) ' s __ Flow RMSE: ' num2str(sol.score.RMSE_flow,'%10.2f\n'), ' m/s, u_Inf: ' num2str(Wp.site.u_Inf,'%10.2f\n') ', v_Inf: ' num2str(Wp.site.v_Inf,'%10.2f\n') ', it. time: ' num2str(sol.score.CPUtime,'%10.2f\n') ' s.']);
        end
    end
end
try 
    if scriptOptions.printProgress
        if strcmp(lower(strucObs.filtertype),'enkf') | strcmp(lower(strucObs.filtertype),'ukf')
            if strucObs.pe.enabled
                for iT = 1:length(strucObs.pe.vars)
                    disp([datestr(rem(now,1)) ' __  t(' num2str(sol.k,['%0' num2str(scriptOptions.tlen) 'd']) ') = ' num2str(sol.time,['%0' num2str(scriptOptions.klen) 'd']) ' s __ ' strucObs.pe.vars{iT} ' estimated as ' num2str(Wp.(strucObs.pe.subStruct{iT}).(strucObs.pe.structVar{iT}),'%10.2f\n') '.']);
                end
            end
        end
        disp(' ')
    end
end
end