function [ flowData,turbData,CT_given ] = loadLESData( scriptOptions,HH )
% Sort and import files from PALM data folder
disp('Sorting and importing files from source folder...')
rawFileList = dir(scriptOptions.sourcePath); % Gather all files from PALM/SOWFA folder
if length(rawFileList) < 3
    error('The specified directory does not exist/has no LES files.');
end
for j = 1:length(rawFileList)-2
    filesInFolder{j} = [rawFileList(j+2).folder '/' rawFileList(j+2).name];   % Remove '.' and '..'
end

% Decide whether this is SOWFA or PALM data
if nnz(cell2mat(strfind(filesInFolder,'.nc'))) > 0
    disp('Found a .nc file. Assuming this is PALM data.')
    [ flowData,turbData ] = loadPALMdata( filesInFolder, HH );
    CT_given              = true;
elseif nnz(cell2mat(strfind(filesInFolder,'.vtk'))) > 0
    % load SOWFA data
    disp('Found a .vtk file. Assuming this is SOWFA data.')
    [ flowData,turbData ] = loadSOWFAdata(filesInFolder,scriptOptions.saveMemory);
    CT_given              = false;
else
    error('Did not find any SOWFA/PALM data in the folder.');
end
end