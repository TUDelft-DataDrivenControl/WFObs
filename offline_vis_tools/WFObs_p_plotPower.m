% Check if figure has been closed
if ishandle(hFigs{2}) == 0 % If figure has been closed
    button = questdlg(['You closed Figure 2: plotPower. Reopen?']);
    if strcmp(button,'Yes')
        disp(['You closed Figure 2: plotPower. Reopening for further use.'])
        hFigs{2} = figure;
    else
        postProcOptions.plotPower = false;
        disp(['You closed Figure 2: plotPower. Disabling for further use.'])
    end
end


% Format all power predictions and measurements as a vector
%[pwLES,pwWFSim] = deal(zeros(Wp.turbine.N,length(sol_array)));
solArrayTurb = [sol_array.turbine];
pwWFSim      = [solArrayTurb.power];
timeWFSim    = Wp.sim.time(2:sol.k+1);
timeLES      = Wp.sim.time(1:Wp.sim.NN);
pwLES        = LESData.turbData.power(1:Wp.sim.NN,:)';

if postProcOptions.powerForecast > 0
    k_end = min(Wp.sim.NN,postProcOptions.powerForecast+sol.k);
    disp(['Forecasting power for k = ' num2str(sol.k+1) ' until k = ' num2str(k_end) '.']);
    sol_tmp    = sol;  % Copy current sol
    sol_tmp.uu = sol.u;
    sol_tmp.vv = sol.v;
    sol_tmp.pp = sol.p;
    timeFC = []; pwFC = [];
    while sol_tmp.k < k_end
        [ sol_tmp,~ ] = WFSim_timestepping( sol_tmp, sys, Wp, postProcOptions );
        pwFC(:,sol_tmp.k-sol.k) = sol_tmp.turbine.power; % Extract power
        timeFC(sol_tmp.k-sol.k) = sol_tmp.time;
    end
    pwFC = movmean(pwFC',[10 0])'; % Low-pass filter power forecast
    t_st = min([postProcOptions.powerForecast, 60]); % t-shortterm
    try
        RMSE_FC_shortterm = sqrt(mean((pwLES(:,timeFC(1:t_st))-pwFC(:,1:t_st)).^2,2));
        RMSE_FC_longterm  = sqrt(mean((pwLES(:,timeFC)-pwFC).^2,2));
        %             RMSE_FC_shortterm = mean(abs(pwLES(:,timeFC(1:t_st))-pwFC(:,1:t_st)),2);
        %             RMSE_FC_longterm  = mean(abs(pwLES(:,timeFC)-pwFC),2);
        disp(['Short-term Mean RMSE forecasted vs. true power for k = ' num2str(sol.k+1) ':' num2str(sol.k+t_st) ' is ' num2str(mean(RMSE_FC_shortterm),'%10.2e\n') '.']);
        disp(['Long-term  Mean RMSE forecasted vs. true power for k = ' num2str(sol.k+1) ':' num2str(k_end) ' is ' num2str(mean(RMSE_FC_longterm),'%10.2e\n') '.']);
    end
end

% Plot results
set(0,'CurrentFigure',hFigs{2}); clf;
subplotDim = numSubplots(Wp.turbine.N); % Determine optimal layout for subplots
for j = 1:Wp.turbine.N
    subplot(subplotDim(1),subplotDim(2),j);
    plot(timeLES,pwLES(j,:),'k-','lineWidth',0.25,'DisplayName', ['LES']); hold on
    plot(timeWFSim,pwWFSim(j,:),'-','lineWidth',1.0,'DisplayName', ['WFObs: ' strucObs.filtertype]); hold on;
    if postProcOptions.powerForecast > 0 && length(timeFC) > 0
        plot(timeFC,pwFC(j,:),'-','lineWidth',0.75,'DisplayName',['WFObs: ' strucObs.filtertype ' (FC)']); hold on;
    end
    legend('-DynamicLegend');
    xlabel('Time (s)');
    xlim([0 Wp.sim.NN]);
    ylim([0 1e6*ceil(max([pwWFSim(:); pwLES(:)])/1e6)]);
    grid on;
    ylabel('Power (W)');
    title(['Turbine ' num2str(j) '']);
end;