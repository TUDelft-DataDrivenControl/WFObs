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
    time     = [sol_array.time];
    RMSE     = [sol_array.RMSE_flow];
    maxError = [sol_array.maxError_flow];
    
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