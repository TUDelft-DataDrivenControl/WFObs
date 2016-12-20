%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 'WFObs_s_reporting.m'
%  This script displays the progress and some primary results in the
%  command line. It includes:
%    * The RMS error in flow estimations (m/s)
%    * The maximum error in flow field estimations (m/s)
%    * The iteration time (s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if k == 1 % filter out all missing SOWFA data
    i_corrected = find(measured.uq ~= 0); 
%     t1xlim = Wp.xline(2)-1;
end;

[maxerroru(k),maxerroruloc(k)] = max(abs(sol.u(i_corrected)-measured.uq(i_corrected)));
[maxerrorv(k),maxerrorvloc(k)] = max(abs(sol.v(i_corrected)-measured.vq(i_corrected)));
[maxerror(k) ]                 = max( [maxerrorv(k) , maxerroru(k)] );

% The next piece of coding shows how "good" the EnKF is.% 'covtracedivbystates' 
% should be comparable to the MSE for good performance.
% Pf_relevant = Pf(1:(Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3),1:(Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3));
% covtracedivbystates(k) = trace(Pf_relevant)/(nrens);
RMSE(k)      = rms([sol.v(i_corrected)-measured.vq(i_corrected);sol.u(i_corrected)-measured.uq(i_corrected)]);
Power(:,k)   = sol.power';
% RMSE_T1(k)   = rms([vec(v(1:t1xlim,:)-measured.vq(1:t1xlim,:)); vec(u(1:t1xlim,:)-measured.uq(1:t1xlim,:))]);
% RMSE_Pw(:,k) = abs(measured.power-Power(:,k)); % Absolute error for turbine power capture

elapsedtime(k) = toc;
disp([datestr(rem(now,1)) ' __  t(' num2str(k,['%0' num2str(tlen) 'd']) ') = ' num2str(timeindex,['%0' num2str(klen) 'd']) ' s __ max. error: ' num2str(maxerror(k),'%10.2f\n'),' m/s, RMSE: ' num2str(RMSE(k),'%10.2f\n'), ' m/s, it. time: ' num2str(elapsedtime(k),'%10.2f\n') ' s.']);

% Reformat power variable and write to file
if strucScript.saveest
    save([strucScript.savepath '\' strucObs.filtertype '_est' num2str(datanroffset+k),'.mat'],'k','sol','Wp'); 
end;