function [ flowDataRaw,turbDataOut ] = loadSOWFAdata( filesInFolder,saveMemory )
addpath('bin/natsortfiles'); % Natural sorter library
filesInFolder = natsortfiles(filesInFolder); % Sort numerically

% Import and verify grid (does not change over time)
[dataType,cellCenters,~] = importVTK(filesInFolder{1});

if saveMemory
    % this way, x is vertical and y is horizontal
    [flowDataRaw.xu,flowDataRaw.xv] = deal(cellCenters(1:2:end,2)); % x is horizontal in SOWFA
    [flowDataRaw.yu,flowDataRaw.yv] = deal(cellCenters(1:2:end,1)); % y is vertical in SOWFA
    z = cellCenters(1:2:end,3);
else
    % this way, x is vertical and y is horizontal
    [flowDataRaw.xu,flowDataRaw.xv] = deal(cellCenters(:,2)); % x is horizontal in SOWFA
    [flowDataRaw.yu,flowDataRaw.yv] = deal(cellCenters(:,1)); % y is vertical in SOWFA
    z = cellCenters(:,3);
end

if abs(min(z) - max(z)) > 1 % This is more than just a horizontal slice
    error('Horizontal slice expected.');
end
    
% Initialize flowDataRaw struct()
NN = nnz(cell2mat(strfind(filesInFolder,'.vtk')));
Nu = length(flowDataRaw.xu);
[flowDataRaw.u,flowDataRaw.v] = deal(zeros(NN,Nu));
flowDataRaw.time              = (1:1:NN)'; % Usually sampled at 1 Hz
flowDataRaw.zu_3d             = unique(z);

% Check if SuperCONOUT file exists
if nnz(cell2mat(strfind(lower(filesInFolder),'uperconout'))) <= 0
    error(['No superCONOUT file found. Please place it' ...
           ' in your source folder, next to the .vtk files.']);
end

% Load all .vtk and SuperCONOUT files and reorganize
ticIt = tic;
NF = length(filesInFolder);
for j = 1:NF
    jFile = filesInFolder{j};
    elapsedTime = toc(ticIt); % Timer and ETA calculator
    if j == 1
        ETA = 'N/A';
    else
        ETA = secs2timestr((NF-j+1)*(elapsedTime/(j-1)));
    end
    disp(['Loading SOWFA file ' num2str(j) '/' num2str(NF) '(' num2str(round(100*(j-1)/NF)) '%) with filename "' jFile '". ETA: ' ETA '.'])
    if ~isempty(strfind(lower(jFile),'.vtk')) % if contains .vtk, load as flow file
        [~,~,cellData]     = importVTK(jFile);
        
        if saveMemory
            flowDataRaw.u(j,:) = cellData(1:2:end,2); % in SOWFA, bottom to top 'v'
            flowDataRaw.v(j,:) = cellData(1:2:end,1); % in SOWFA, left to right is 'u'
            % Formatted this way: u is bottom -> top, v is left -> right            
        else
            flowDataRaw.u(j,:) = cellData(:,2); % in SOWFA, bottom to top 'v'
            flowDataRaw.v(j,:) = cellData(:,1); % in SOWFA, left to right is 'u'
            % Formatted this way: u is bottom -> top, v is left -> right
        end
        
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
    turbDataOut.Mz(:,j)       = turbDataRaw.data{j}(:,37);
    turbDataOut.phi(:,j)      = turbDataRaw.data{j}(:,21);
    turbDataOut.yawerror(:,j) = turbDataRaw.data{j}(:,22);
    turbDataOut.gentorq(:,j)  = turbDataRaw.data{j}(:,4);
    turbDataOut.genspeed(:,j) = turbDataRaw.data{j}(:,3);
    turbDataOut.power(:,j)    = turbDataRaw.data{j}(:,2);
    turbDataOut.pitch(:,j)    = mean(turbDataRaw.data{j}(:,6:8),2); % Collective/avg
end
end