function [ flowDataRaw,turbDataOut ] = loadSOWFAdata( filesInFolder,scriptOptions )
% Load first .vtk slice to get the mesh settings and atmospheric settings
[dataType,cellCenters,cellData] = importVTK(filesInFolder{1});
flowDataRaw.u = cellData(:,2); % in SOWFA, bottom to top 'v'
flowDataRaw.v = cellData(:,1); % in SOWFA, left to right is 'u'

% this way, x is vertical and y is horizontal
[flowDataRaw.xu,flowDataRaw.xv] = deal(cellCenters(:,2)); % x is horizontal in SOWFA
[flowDataRaw.yu,flowDataRaw.yv] = deal(cellCenters(:,1)); % y is vertical in SOWFA

z = cellCenters(:,3);
if abs(min(z) - max(z)) > 1 % This is more than just a horizontal slice
    error('Horizontal slice expected.');
end
    
% Initialize flowDataRaw struct()
NN                = nnz(cell2mat(strfind(filesInFolder,'.vtk')));
flowDataRaw.time  = (0:NN-1)'; % Usually sampled at 1 Hz
flowDataRaw.zu_3d = unique(z);

% load superCONOUT file
try 
    tmp_sco = strfind(lower(filesInFolder),'uperconout');
    jFile = filesInFolder{find(~cellfun(@isempty,tmp_sco))};
    turbDataRaw = importSuperCONOUT(jFile);
catch
    error(['No superCONOUT file found. Please place it' ...
           ' in your source folder next to the .vtk files and call it ''superCONOUT.csv''.']);
end

disp('Imported everything. Process turbDataRaw')
turbDataOut.time  = turbDataRaw.time-turbDataRaw.time(1); % Start at 0
for j = 1:length(turbDataRaw.data)
    turbDataOut.power(:,j)    = turbDataRaw.data{j}(:,2);    
    turbDataOut.Mz(:,j)       = turbDataRaw.data{j}(:,37);
    turbDataOut.phi(:,j)      = turbDataRaw.data{j}(:,21);
    turbDataOut.yawerror(:,j) = turbDataRaw.data{j}(:,22);
    turbDataOut.gentorq(:,j)  = turbDataRaw.data{j}(:,4);
    turbDataOut.genspeed(:,j) = turbDataRaw.data{j}(:,3);
    turbDataOut.pitch(:,j)    = mean(turbDataRaw.data{j}(:,6:8),2); % Collective/avg
end
end