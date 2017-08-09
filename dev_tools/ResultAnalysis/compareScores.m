clear all;

data = {};

%% OPEN-LOOP SIMS
% data{end+1} = struct(...
%     'name','Sim (no param tuning)',...
%     'path','../../results/sim_noParamTuning/workspace.mat');

%% Periodic parameter tuning
% data{end+1} = struct(...
%     'name','Sim (with param tuning)',...
%     'path','../../results/sim_paramTuning/workspace.mat');
data{end+1} = struct(...
    'name','EnKF (with param tuning)',...
    'path','../../results/YawCase3_EnKF_n50_periodicParamEstimation/workspace.mat');
data{end+1} = struct(...
    'name','EnKF (nominal)',...
    'path','../../results/YawCase3_EnKF_n50_nominal/workspace.mat');


%% UKF: FEW MEASUREMENTS
% data{end+1} =  struct(...
%     'name','UKF + no param. estimation + 20% measurements',...
%     'path','../../results/YawCase3_ukf_noParamEst/workspace.mat');
% data{end+1} =  struct(...
%     'name','UKF + param. estimation + 20% measurements',...
%     'path','../../results/YawCase3_ukf_paramEst/workspace.mat');
% data{end+1} =  struct(...
%     'name','UKF + no state estimation + 20% measurements',...
%     'path','../../results/YawCase3_ukf_noStateEst/workspace.mat');

%% UKF: MANY MEASUREMENTS
% data{end+1} =  struct(...
%     'name','UKF + no param. estimation + 20% measurements',...
%     'path','../../results/results_manyMeasurements/YawCase3_ukf_noParamEst/workspace.mat');
% data{end+1} =  struct(...
%     'name','UKF + param. estimation + 20% measurements',...
%     'path','../../results/results_manyMeasurements/YawCase3_ukf_paramEst/workspace.mat');
% data{end+1} =  struct(...
%     'name','UKF + no state estimation + 20% measurements',...
%     'path','../../results/results_manyMeasurements/YawCase3_ukf_noStateEst/workspace.mat');

%% ENKF
%% N = 50 and 2% (downstream) measurements
% data{end+1} =  struct(...
%     'name','EnKF (n=50) + no param. estimation + 2% measurements',...
%     'path','../../results/YawCase3_enkf_n50_noParamEst/workspace.mat');
% data{end+1} =  struct(...
%     'name','EnKF (n=50) + param. estimation + 2% measurements',...
%     'path','../../results/YawCase3_enkf_n50_paramEst/workspace.mat');
% data{end+1} =  struct(...
%     'name','EnKF (n=50) + no state estimation + 2% measurements',...
%     'path','../../results/YawCase3_enkf_n50_noStateEst/workspace.mat');
%% N = 200 and 2% (downstream) measurements
% data{end+1} =  struct(...
%     'name','EnKF (n=200) + no param. estimation + 2% measurements',...
%     'path','../../results/YawCase3_enkf_n200_noParamEst/workspace.mat');
% data{end+1} =  struct(...
%     'name','EnKF (n=200) + param. estimation + 2% measurements',...
%     'path','../../results/YawCase3_enkf_n200_paramEst/workspace.mat');
% % data{end+1} =  struct(...
%     'name','EnKF (n=200) + no state estimation + 2% measurements',...
%     'path','../../results/YawCase3_enkf_n200_noStateEst/workspace.mat');
%% N = 50 and 20% measurements
% data{end+1} =  struct(...
%     'name','EnKF (n=50) + no param. estimation + 20% measurements',...
%     'path','../../results/results_manyMeasurements/YawCase3_enkf_n50_noParamEst/workspace.mat');
% data{end+1} =  struct(...
%     'name','EnKF (n=50) + param. estimation + 20% measurements',...
%     'path','../../results/results_manyMeasurements/YawCase3_enkf_n50_paramEst/workspace.mat');
% data{end+1} =  struct(...
%     'name','EnKF (n=50) + no state estimation + 20% measurements',...
%     'path','../../results/results_manyMeasurements/YawCase3_enkf_n50_noStateEst/workspace.mat');
%% N = 200 and 20% measurements
% data{end+1} =  struct(...
%     'name','EnKF (n=200) + no param. estimation + 20% measurements',...
%     'path','../../results/results_manyMeasurements/YawCase3_enkf_n200_noParamEst/workspace.mat');
% data{end+1} =  struct(...
%     'name','EnKF (n=200) + param. estimation + 20% measurements',...
%     'path','../../results/results_manyMeasurements/YawCase3_enkf_n200_paramEst/workspace.mat');
% data{end+1} =  struct(...
%     'name','EnKF (n=200) + no state estimation + 20% measurements',...
%     'path','../../results/results_manyMeasurements/YawCase3_enkf_n200_noStateEst/workspace.mat');

for di = 1:length(data)
    WS_tmp = load(data{di}.path);
    
    % collect all fieldnames and save to data{*} struct
    scoreFields = fieldnames(WS_tmp.sol_array{1}.score);
    for si = 1:length(scoreFields)
        data{di}.(scoreFields{si}) = zeros(1,length(WS_tmp.sol_array));
        for t = 1:length(WS_tmp.sol_array) % collect all data
            data{di}.(scoreFields{si})(t) = WS_tmp.sol_array{t}.score.(scoreFields{si});
        end
    end
    
    clear WS_tmp
end


% Produce figure and print mean values
close all; clc;
for si = 1:length(scoreFields)
    hFigs{si} = figure;
    for di = 1:length(data)
        plot(data{di}.(scoreFields{si}),'DisplayName',data{di}.name);
        hold on;
    end
    title(['Time series for ' scoreFields{si}])
    legend('-dynamicLegend');
    xlabel('Time (k)');
    grid on;
    
    disp(['Mean(' scoreFields{si} '):'])
    for di = 1:length(data)
        disp(['  ' data{di}.name ': ' num2str(mean(data{di}.(scoreFields{si}))) '.']);
    end
    disp(' ')
end
