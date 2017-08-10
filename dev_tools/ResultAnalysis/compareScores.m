clear all;

data = {};

%% OPEN-LOOP SIMS
% data{end+1} = struct(...
%     'name','Sim (nominal)',...
%     'path','../../results/YawCase3_sim_nominal/workspace.mat');
% data{end+1} = struct(...
%     'name','Sim (tuned)',...
%     'path','../../results/YawCase3_sim_tuned/workspace.mat');

% Periodic parameter tuning
data{end+1} = struct(...
    'name','UKF (nominal)',...
    'path','../../results/YawCase3_UKF_nominal/workspace.mat');
data{end+1} = struct(...
    'name','UKF (tuned)',...
    'path','../../results/YawCase3_UKF_tuned/workspace.mat');
% data{end+1} = struct(...
%     'name','EnKF (nominal)',...
%     'path','../../results/YawCase3_EnKF_nominal/workspace.mat');
% data{end+1} = struct(...
%     'name','EnKF (tuned)',...
%     'path','../../results/YawCase3_EnKF_tuned/workspace.mat');



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
