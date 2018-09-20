clear all; close all;

data = {};

%% OPEN-LOOP SIMS
% data{end+1} = struct(...
%     'name','sim',...
%     'path','../../results/WE2017/axi_2turb_alm_turb_sim/workspace.mat');
% data{end+1} = struct(...
%     'name','ExKF',...
%     'path','../../results/WE2017/axi_2turb_alm_turb_measFlow_exkf_stateEst/workspace.mat');
% data{end+1} = struct(...
%     'name','EnKF',...
%     'path','../../results/WE2017/axi_2turb_alm_turb_measFlow_enkf_stateEst/workspace.mat');
% data{end+1} = struct(...
%     'name','UKF',...
%     'path','../../results/WE2017/axi_2turb_alm_turb_measFlow_ukf_stateEst/workspace.mat');

mainDir = 'D:\bmdoekemeijer\My Documents\SurfDrive\PhD\Dissemination\2017 Wind Energy\MATLAB_out';
data{end+1} = struct(...
    'name','sim',...
    'path',[mainDir '/apc_9turb_alm_turb_simPoorLmu/workspace.mat');
data{end+1} = struct(...
    'name','ExKF',...
    'path','../../results/WE2017/axi_2turb_alm_turb_measFlow_exkf_stateEst/workspace.mat');
data{end+1} = struct(...
    'name','EnKF',...
    'path','../../results/WE2017/axi_2turb_alm_turb_measFlow_enkf_stateEst/workspace.mat');
data{end+1} = struct(...
    'name','UKF',...
    'path','../../results/WE2017/axi_2turb_alm_turb_measFlow_ukf_stateEst/workspace.mat');

for di = 1:length(data)
    WS_tmp = load(data{di}.path);
    
    % collect all fieldnames and save to data{*} struct
    scoreFields = fieldnames(WS_tmp.sol_array(1).score);
    for si = 1:length(scoreFields)
        data{di}.(scoreFields{si}) = zeros(1,length(WS_tmp.sol_array));
        for t = 1:length(WS_tmp.sol_array) % collect all data
            data{di}.(scoreFields{si})(t) = WS_tmp.sol_array(t).score.(scoreFields{si});
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
