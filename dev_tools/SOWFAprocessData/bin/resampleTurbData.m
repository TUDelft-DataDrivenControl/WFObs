function [ turbDataOut ] = resampleTurbData( turbData, tau, time_target )

fieldNamesListTurb = fieldnames(turbData);

% Check sampling time from raw data
if var(diff(turbData.time)) < 1e-5
     dtRawTurb = mean(diff(turbData.time));
else
    disp('Detected time-varying timestep sizes in turbData. Resampling at smallest dt.');
    dtRawTurb = min(diff(turbData.time));
    time_new  = turbData.time(1):dtRawTurb:turbData.time(end);
    for jField = 1:length(fieldNamesListTurb)
        jField = fieldNamesListTurb{jField};
        for jTurb = 1:size(turbData.(jField),2)
            turbDataOut.(jField)(:,jTurb) = interp1(turbData.time,turbData.(jField)(:,jTurb),time_new);
        end
    end
    turbData = turbDataOut; 
    clear turbDataOut
end

% Create low pass filter function
t_lpf = [1:length(turbData.time)]*dtRawTurb;
LPF   = c2d(1/(tau*tf('s')+1),dtRawTurb,'tustin');

orderSubPlots      = numSubplots(length(fieldNamesListTurb)-1);

% Low-pass filter all variables for each turbine
figure % create empty figure
jPlot = 0;
for j = 1:length(fieldNamesListTurb)
    jField = fieldNamesListTurb{j};
    if strcmp(jField,'time') % add any non-filtered fields here
        turbDataOut.(jField) = time_target;
    else
        jPlot = jPlot+1;
        subplot(orderSubPlots(1),orderSubPlots(2),jPlot);
        for jTurb = 1:size(turbData.(jField),2)
            y_in  = turbData.(jField)(:,jTurb);
            y_fil = lsim(LPF,y_in,t_lpf);
            
            % Resample to desired time
            y_out = interp1(turbData.time,y_fil,time_target,'linear');
            turbDataOut.(fieldNamesListTurb{j})(:,jTurb) = y_out;
            
            % Plot results
            plot(t_lpf,y_in,'displayName',['raw T' num2str(jTurb)]); hold on;
            plot(time_target,y_out,'--','displayName',['filtered T' num2str(jTurb)]);
            legend('-DynamicLegend'); grid on;
            xlabel('Time (s)');
            ylabel(fieldNamesListTurb{j});
        end
    end
end
drawnow;
end