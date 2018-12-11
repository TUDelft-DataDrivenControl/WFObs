function [filesInFolder,flowData,turbData,dataSource] = loadLESData(scriptOptions,rawTurbData)
% Sort and import files from data folder
disp('Sorting and importing files from source folder...')
rawFileList = dir(scriptOptions.sourcePath); % Gather all files from PALM/SOWFA folder
if length(rawFileList) < 3
    error('The specified directory does not exist/is empty.');
end

% If applicable: cd into subfolder called sliceDataInstantaneous
if any(strcmp({rawFileList.name},'sliceDataInstantaneous'))
    rawFileList = dir([scriptOptions.sourcePath filesep 'sliceDataInstantaneous'])
end

% Sort files/folder names and determine how many folders are in directory
nrFolders = 0;
for j = 1:length(rawFileList)-2
    filesInFolder{j} = [rawFileList(j+2).folder filesep rawFileList(j+2).name];   % Remove '.' and '..'
    if exist(filesInFolder{j}) == 7
        nrFolders = nrFolders + 1;
    end
end

% More than 90% are folders: must be subdirectories...?
bufferNrFldrs = 0.90;
if nrFolders >= bufferNrFldrs * length(filesInFolder)
    disp([num2str(bufferNrFldrs*100) '% of the objects in specified folder are folders. Searching subdirectories for VTKs.'])
    rawFileListSubdir = dir([filesInFolder{2} filesep '*.vtk']); % Taking 2nd directory
    if length(rawFileListSubdir) == 1
        horVtkIndx = 1;
    else
        horVtkIndx = find(~cellfun(@isempty,strfind({rawFileListSubdir.name},'hor')));
        if length(horVtkIndx) ~= 1
            error('  Could not find a unique VTK file with ''hor'' in its name.');
        end        
    end
    vtkFileName = rawFileListSubdir(horVtkIndx).name;
    disp(['  Assuming the resp. .vtk filename is ''' vtkFileName ''' .']);
    disp(['  Sorting VTK files by indexing the subdirectories.']);
    filesInFolderNew = {};
    for j = 1:length(filesInFolder)
        fnameTmp = [filesInFolder{j} filesep vtkFileName];
        if exist(fnameTmp) == 2
            filesInFolderNew{end+1} = fnameTmp;
        end
    end
    filesInFolder = filesInFolderNew; % overwrite
end

filesInFolder = natsortfiles(filesInFolder); % Sort numerically
% Decide whether this is SOWFA or PALM data
if nnz(cell2mat(strfind(filesInFolder,'.nc'))) > 0
    disp('Found a .nc file. Assuming this is PALM data.')
    % Load all flow data and turbine data at once
    [ flowData,turbData ] = loadPALMdata(filesInFolder,rawTurbData.hubHeight);
    dataSource            = 'palm';
elseif ( nnz(cell2mat(strfind(filesInFolder,'.vtk'))) > 0 )
    disp('Found a .vtk file. Assuming this is SOWFA data.')
    % Load only flow data mesh and complete turbine data
    [ flowData,turbData ] = loadSOWFAdata(filesInFolder,scriptOptions);
    dataSource            = 'sowfa';
else
    error('Did not find any SOWFA/PALM data in the folder.');
end
end