function [inputData,measuredData] = getDataLES_LegacyFormat(LESData,measOptions,currentTime)

    % Load and format measurements from preloaded LES database (legacy format)
    measuredData = [];
    if measOptions.P.enabled % Setup power measurements
        allTurbinePowers = LESData.PwrInterpolant(currentTime);
        for i = 1:length(measOptions.P.turbIds)
            measuredData(i).idx   = measOptions.P.turbIds(i); % Turbine number
            measuredData(i).type  = 'P'; % Power measurement
            measuredData(i).value = allTurbinePowers(measuredData(i).idx) + ...
                                    measOptions.P.noiseStd*randn(); % With artificial noise
            measuredData(i).std   = measOptions.P.measStd; % Standard deviation in W
        end
    end
    
    if measOptions.U.enabled % Setup flow measurements
        for jT = [1 2] % jT=1 for 'u', jT=2 for 'v' measurements
            iOffset = length(measuredData);
            for i = 1:size(measOptions.U.locations,1)
                measuredData(iOffset+i).idx = measOptions.U.locations(i,:); % Sensor location
                measuredData(iOffset+i).std = measOptions.U.measStd; % Standard deviation in W
                
                if jT == 1
                    measuredData(iOffset+i).type  = 'u'; % Long. flow measurement
                    measuredData(iOffset+i).value = LESData.uInterpolant(currentTime,...
                        measOptions.U.locations(i,1), measOptions.U.locations(i,2)) ...
                        + measOptions.U.noiseStd*randn(); % Plus artificial noise
                elseif jT == 2
                    measuredData(iOffset+i).type  = 'v'; % Lat. flow measurement
                    measuredData(iOffset+i).value = LESData.vInterpolant(currentTime,...
                        measOptions.U.locations(i,1), measOptions.U.locations(i,2)) ...
                        + measOptions.U.noiseStd*randn(); % Plus artificial noise
                end
            end
        end
    end
    
    % Gather turbine inputs from legacy LES file database
    timeIndex = interp1([LESData.turb.inputs.t],1:length(LESData.turb.inputs),currentTime,'nearest');
    inputData = LESData.turb.inputs(timeIndex);
end

