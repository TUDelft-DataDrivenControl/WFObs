clear all; clc; close all;

% SETUP: Configuration file location
configurationName = 'sowfa_9turb_ACC2019'; % use '-all' to do a batch run over all configurations

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
    [ flow,turb ] = LES_import_core( scriptOptions,rawTurbData,filterSettings,meshSetup );

    % Save output data for configuration [i]
    disp('Saving exported data...');
    save(['exportedData/' scriptOptions.outputFilename],'flow','turb')
end