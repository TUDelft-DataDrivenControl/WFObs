function [ turbDataOut ] = resampleTurbData( turbDataIn, filterSettings, time_target )

fieldNamesListTurb = fieldnames(turbDataIn);
NF                 = length(fieldNamesListTurb); % Number of datapoints
NT                 = size(turbDataIn.power,2);   % Number of turbines
turbDataIn.time    = turbDataIn.time-turbDataIn.time(1); % Start at 0

% Resample data to smallest dt, if necessary
if var(diff(turbDataIn.time)) < 1e-6
    dtRaw       = mean(diff(turbDataIn.time));
    turbDataRaw = turbDataIn;
else
    disp('Detected time-varying timestep sizes in turbData. Resampling at smallest dt.');
    dtRaw = min(diff(turbDataIn.time));
    time_raw  = turbDataIn.time(1):dtRaw:turbDataIn.time(end);
    for jField = 1:NF
        jField = fieldNamesListTurb{jField};
        if strcmp(jField,'time')
            turbDataRaw.time = time_raw;
        else
            for jTurb = 1:NT
                turbDataRaw.(jField)(:,jTurb) = interp1(turbDataIn.time,turbDataIn.(jField)(:,jTurb),time_raw);
            end
        end
    end
end

% create empty figure and calculate optimal dimensions
figure; jPlot = 0; 
np = numSubplots(NF-1);

% average and resample data
kl = ceil(filterSettings.tL/dtRaw); % width of sliding window (left/backward)
kr = ceil(filterSettings.tR/dtRaw); % width of sliding window (right/forward)
turbDataOut.time = time_target;
for j = 1:NF
    jField = fieldNamesListTurb{j};
    if strcmp(jField,'time') == false
        jPlot = jPlot+1;
        subplot(np(1),np(2),jPlot);
        for jTurb = 1:NT
            y_in  = turbDataRaw.(jField)(:,jTurb);
            if filterSettings.MM 
                y_in  = movmean(y_in,[kl kr]); % moving mean average
            end
            y_out = interp1(turbDataRaw.time,y_in,time_target,'linear'); % Resample to desired time
            turbDataOut.(jField)(:,jTurb) = y_out;
            
            % Plot results
            plot(turbDataRaw.time,y_in, '-', 'displayName',['raw T' num2str(jTurb)]); hold on;
            plot(turbDataOut.time,y_out,'--','displayName',['filtered T' num2str(jTurb)]);
            legend('-DynamicLegend'); grid on;
            xlabel('Time (s)');
            ylabel(jField);
        end
    end
end
drawnow

end