function [ flowDataRaw,turbDataOut ] = loadSOWFAdata( filesInFolder,scriptOptions )
% Load first .vtk slice to get the mesh settings and atmospheric settings
[dataType,cellCenters,cellData] = importVTK(filesInFolder{1});
flowDataRaw.u = cellData(:,1); % in SOWFA, left to right is 'u'
flowDataRaw.v = cellData(:,2); % in SOWFA, bottom to top 'v'

% this way, x is vertical and y is horizontal
[flowDataRaw.xu,flowDataRaw.xv] = deal(cellCenters(:,1)); % x is horizontal in SOWFA
[flowDataRaw.yu,flowDataRaw.yv] = deal(cellCenters(:,2)); % y is vertical in SOWFA

z = cellCenters(:,3);
if abs(min(z) - max(z)) > 1 % This is more than just a horizontal slice
    error('Horizontal slice expected.');
end

% Initialize flowDataRaw struct()
NN = nnz(cell2mat(strfind(filesInFolder,'.vtk')));
tmpString =regexp(filesInFolder{1},'/','split');
vtkName = tmpString{end};
vtkName = strrep(vtkName,'.vtk',''); % Remove the VTK part
vtkId = str2num(vtkName);
dtFlow = rem(vtkId,1000);
disp('Assuming the flow is sampled in steps of 5 seconds (based on VTK names.');
flowDataRaw.time  = dtFlow*(1:NN'); % Usually sampled at 1 Hz
flowDataRaw.zu_3d = unique(z);

% load superCONOUT file
try
    tmp_sco = strfind(lower(filesInFolder),'uperconout');
    jFile = filesInFolder{find(~cellfun(@isempty,tmp_sco))};
    turbDataRaw = importSuperCONOUT(jFile);

    disp('Imported everything. Process turbDataRaw')
    turbDataOut.time  = turbDataRaw.time-turbDataRaw.time(1); % Start at 0
    for j = 1:length(turbDataRaw.data)
        turbDataOut.power(:,j)    = turbDataRaw.data{j}(:,2);
        turbDataOut.Mz(:,j)       = turbDataRaw.data{j}(:,37);
        turbDataOut.phi(:,j)      = turbDataRaw.data{j}(:,21);
%         turbDataOut.yawerror(:,j) = turbDataRaw.data{j}(:,22);
%         turbDataOut.gentorq(:,j)  = turbDataRaw.data{j}(:,4);
%         turbDataOut.genspeed(:,j) = turbDataRaw.data{j}(:,3);
%         turbDataOut.pitch(:,j)    = mean(turbDataRaw.data{j}(:,6:8),2); % Collective/avg
    end
catch
    disp('No superCONOUT file found. Checking for ADMR files.');
    
    if isfield(scriptOptions,'turbOutADM')
        disp('Found turbOutADM! If this is wrong, comment out ''scriptOptions.turbOutADM = ...'' in your configuration file.');
        nacYaw = importADMR([scriptOptions.turbOutADM '/nacYaw']);
        power =  importADMR([scriptOptions.turbOutADM '/powerGenerator']);
        turbDataOut.time = unique(nacYaw.Times);
        turbDataOut.time(isnan(turbDataOut.time)) = [];
        turbDataOut.time = turbDataOut.time-turbDataOut.time(1);
        for iTurb = 1:(max(nacYaw.Turbine)+1)
            turbDataOut.phi(:,iTurb) = nacYaw.data(nacYaw.Turbine==iTurb-1,:) - 270.0;
            turbDataOut.power(:,iTurb) = power.data(nacYaw.Turbine==iTurb-1,:);
        end
        disp('WARNING: ASSUMING ALL CT_PRIME VALUES TO BE GREEDY: 2.0');
        turbDataOut.CT_prime = 2.0*ones(size(turbDataOut.phi));
    else
        disp('No superCONOUT file found. Exporting flow only.');
        turbDataOut = [];
    %     error(['No superCONOUT file found. Please place it' ...
    %            ' in your source folder next to the .vtk files and call it ''superCONOUT.csv''.']);
    end
end
end