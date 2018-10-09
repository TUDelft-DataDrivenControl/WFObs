% Check if figure has been closed
if ishandle(hFigs{3}) == 0 % If figure has been closed
    button = questdlg(['You closed Figure 3: plotError. Reopen?']);
    if strcmp(button,'Yes')
        disp(['You closed Figure 3: plotError. Reopening for further use.'])
        hFigs{3} = figure;
        
    else
        postProcOptions.plotError = false;
        disp(['You closed Figure 3: plotError. Disabling for further use.'])
    end
end

if postProcOptions.plotError
    % Format all error scores as a vector
    [RMSE,maxError] = deal(zeros(1,length(sol_array)));
    for jt = 1:length(sol_array)
        time(jt)     = sol_array(jt).time;
        RMSE(jt)     = sol_array(jt).score.RMSE_flow;
        maxError(jt) = sol_array(jt).score.maxError_flow;
    end;
    % Plot results
    set(0,'CurrentFigure',hFigs{3}); clf;
    plot(time,maxError,'DisplayName','Max. error (m/s)'); hold on;
    plot(time,RMSE,'--','DisplayName','RMS error (m/s)'); title('Error in velocity estimates');
    xlabel('Time (s)');
    ylabel('Error (m/s)');
    grid on;
    ylim([0 5]);
    % xlim([0 Wp.sim.NN]);
    legend('-DynamicLegend');
    
    if postProcOptions.savePlots
        drawnow;
        saveas(hFigs{3},[postProcOptions.savePath '/' strucObs.filtertype '_errorplot.png']);
    end;
end