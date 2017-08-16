clear all; clc; close all;

configurationName = '-all'; % use '-all' to do a batch run over all configurations

% Organize file names
if strcmp(lower(configurationName),'-all')
    configurations = dir('configurations');
    configurations = {configurations(3:end).name};
else
    configurations = {configurationName};
end
 
addpath('bin');
mkdir('exportedData');
for i = 1:length(configurations)
    disp(['Simulating configuration ' num2str(i) '/' num2str(length(configurations)) ': ' configurations{i} '.'])
    run(['configurations/' configurations{i}])
    
    % Core code
    [ time,u,v,turbData,meshingOut ] = LES_import_core( scriptOptions,rawTurbData,meshSetup );

    % Save output data
    disp('Saving flow and meshing information...');
    save(['exportedData/' scriptOptions.outputFilename '_data.mat'],'time','u','v','turbData')
    save(['exportedData/' scriptOptions.outputFilename '_meshing.mat'],'-struct','meshingOut')
end