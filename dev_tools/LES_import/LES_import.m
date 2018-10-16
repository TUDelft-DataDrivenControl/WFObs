clear all; clc; close all;

% SETUP: Configuration file location
configurationName = 'palm_4turb_adm_turb'; % use '-all' to do a batch run over all configurations


% CORE: Loop over all specified configuration files
addpath('bin');
mkdir('exportedData');
if strcmp(lower(configurationName),'-all')
    % Organize configurations into a list of filenames
    configurations = dir('configurations');
    configurations = {configurations(3:end).name};
else
    configurations = {configurationName};
end
for i = 1:length(configurations)
    disp(['Simulating configuration ' num2str(i) '/' num2str(length(configurations)) ': ' configurations{i} '.'])
    run(['configurations/' configurations{i}])
    
    % Execute export for configuration [i]
    clear u v turbData meshingOut
    [ time,u,v,turbData,meshingOut ] = LES_import_core( scriptOptions,rawTurbData,filterSettings,meshSetup );

    % Save output data for configuration [i]
    disp('Saving flow and meshing information...');
    save(['exportedData/' scriptOptions.outputFilename '_data.mat'],'time','u','v','turbData')
    save(['exportedData/' scriptOptions.outputFilename '_meshing.mat'],'-struct','meshingOut')
end