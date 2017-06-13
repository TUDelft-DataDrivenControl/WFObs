clear all; tic;
TargetFolder = 'D:\Documents\MATLAB\WFObs\WFSim\data_SOWFA\WithPrecursor\2turb_50x25_lin';
SuperCONpath = 'D:\Documents\MATLAB\WFObs\WFSim\data_SOWFA\WithPrecursor\superCONOUT.csv';
dataOffset   = 20000;
dataRange    = 1:1998;
dataPrefix   = '';

%% Core code
disp('Importing SCO data...');
SCO=importSuperCONOUT(SuperCONpath);

% Format sensorlist to compatible variable names
sensorList = SCO.sensorList;
sensorList = strrep(sensorList,' ','');  % remove blank spaces
sensorList = strrep(sensorList,'-','_'); % replace dashed lines

% Extract sensor names and units
disp('Extracting sensor names and units from SCO data...');
for j = 1:length(sensorList)
    unitsIdx       = strfind(sensorList{j},'(');
    if length(unitsIdx) > 0
        sensorListName{j} = sensorList{j}(1:unitsIdx-1);
        sensorListUnit{j} = sensorList{j}(unitsIdx+1:end-1);
    else
        sensorListName{j} = sensorList{j};
        sensorListUnit{j} = '_';
    end;
end;

disp('Appending files with turb data...'); tic;
for k = dataRange
    clear turb
    turb.units = sensorListUnit;
    
    % Iterate over all field values from SCO
    for j = 1:length(sensorListName)
        
        % Write all turbine measurements at time k to one struct
        for jj = 1:length(SCO.data)
            turb.data.(sensorListName{j})(jj) = SCO.data{jj}(k*50,j);
        end;
    end;
    
    % Append to existing mat files
    save([TargetFolder '\' dataPrefix num2str(k+dataOffset) '.mat'],'turb','-append');
    
end;
toc; disp(['Appended turb data to ' num2str(length(dataRange)) ' files.']);
   