function [ flowDataRaw,turbDataOut ] = loadSOWFAdata( filesInFolder )
addpath('natsortfiles'); % Natural sorter library
filesInFolder = natsortfiles(filesInFolder); % Sort numerically

% Import and verify grid (does not change over time)
[dataType,cellCenters,~] = importVTK(filesInFolder{1});
[flowDataRaw.xu,flowDataRaw.xv] = deal(cellCenters(:,2)); % x is horizontal in SOWFA
[flowDataRaw.yu,flowDataRaw.yv] = deal(cellCenters(:,1)); % y is vertical in SOWFA
% this way, x is vertical and y is horizontal
z = cellCenters(:,3);
if abs(min(z) - max(z)) > 1 % This is more than just a horizontal slice
    error('Horizontal slice expected.');
end
    
% Initialize flowDataRaw struct()
NN = nnz(cell2mat(strfind(filesInFolder,'.vtk')));
Nu = length(flowDataRaw.xu);
[flowDataRaw.u,flowDataRaw.v] = deal(zeros(NN,Nu));
flowDataRaw.time              = (1:1:NN)'; % Usually sampled at 1 Hz
flowDataRaw.zu_3d             = unique(z);

% Load all .vtk and SuperCONOUT files and reorganize
for j = 1:length(filesInFolder)
    jFile = filesInFolder{j};
    disp(['Loading SOWFA file ' jFile ' ...'])
    if ~isempty(strfind(jFile,'.vtk')) % if contains .vtk, load as flow file
        [~,~,cellData]     = importVTK(jFile);
        flowDataRaw.u(j,:) = cellData(:,2); % in SOWFA, bottom to top 'v'
        flowDataRaw.v(j,:) = cellData(:,1); % in SOWFA, left to right is 'u'
        % Formatted this way: u is bottom -> top, v is left -> right
    elseif ~isempty(strfind(lower(jFile),'uperconout')) % superCONOUT file
        if exist('SCO')
            error('Multiple SuperCONOUT files detected.')
        else
            turbDataRaw = importSuperCONOUT(jFile);
        end
    else
        disp(['Skipping file with name ' jFile '.'])
    end
end

if exist('turbDataRaw') == 0
    error(['No superCONOUT file found. Please place it' ...
           ' in your source folder, next to the .vtk files.']);
end
disp('Imported everything. Process turbDataRaw')
turbDataOut.time  = turbDataRaw.time;
for j = 1:length(turbDataRaw.data)
%     turbDataOut.Ur(:,j)    = turbDataRaw.data{j}(:,1);
%     turbDataOut.CT         = M(:,4);
%     turbDataOut.a          = M(:,5);
    turbDataOut.phi(:,j)      = turbDataRaw.data{j}(:,21);
    turbDataOut.yawerror(:,j) = turbDataRaw.data{j}(:,22);
    turbDataOut.gentorq(:,j)  = turbDataRaw.data{j}(:,4);
    turbDataOut.genspeed(:,j) = turbDataRaw.data{j}(:,3);
    turbDataOut.power(:,j)    = turbDataRaw.data{j}(:,2);
    turbDataOut.pitch(:,j)    = mean(turbDataRaw.data{j}(:,6:8),2); % Collective/avg
end
end