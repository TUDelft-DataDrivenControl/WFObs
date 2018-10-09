function [LESData] = loadLESdata(fileLocation)
    LESData = load(fileLocation);
    
    % Set-up interpolants
    [T,XU,YU] = ndgrid(LESData.flow.time',LESData.flow.xu(:,1),LESData.flow.yu(1,:)');
    [T,XV,YV] = ndgrid(LESData.flow.time',LESData.flow.xv(:,1),LESData.flow.yv(1,:)');
    LESData.uInterpolant = griddedInterpolant(T,XU,YU,LESData.flow.u,'linear','none');
    LESData.vInterpolant = griddedInterpolant(T,XV,YV,LESData.flow.v,'linear','none');
    
    LESData.PwrInterpolant = @(currentTime) LinearPwrInterpolation(currentTime, LESData);
    
    function PowerOut = LinearPwrInterpolation(currentTime,LESData)
        timeIndex = interp1(LESData.turb.rawData.time,1:length(LESData.turb.rawData.time),currentTime,'nearest','none');
        PowerOut = LESData.turb.rawData.power(timeIndex,:);
    end
end