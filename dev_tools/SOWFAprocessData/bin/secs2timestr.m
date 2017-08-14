function [ timestring ] = secs2timestr( t )
% Convert seconds to *h *m *s format
    t     = round(t); % Remove any ms
    hours = floor(t / 3600);
    t     = t - hours * 3600;
    mins  = floor(t / 60);
    secs  = t - mins * 60;
    
    % Construct string
    timestring = [num2str(secs,'%02d') 's'];
    if mins > 0
        timestring = [num2str(mins,'%02d') 'm ' timestring];
    end    
    if hours > 0
        timestring = [num2str(hours) 'h ' timestring];
    end
end