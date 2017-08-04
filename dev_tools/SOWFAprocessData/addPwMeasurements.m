clear all; tic;
TargetFolder        = 'D:\bmdoekemeijer\My Documents\MATLAB\WFObs\WFSim\Data_SOWFA\WithPrecursor';
TargetFolder_subdir = '2turb_50x25_lin';
dataOffset          = 20000;
dataRange           = 1:1998;
dataPrefix          = '';

SCO=importSuperCONOUT([TargetFolder '/SuperCONOUT.csv']);
for j = 1:length(SCO.data)
    powerArray(j,:) = SCO.data{j}(50:50:end,2);
end;
for j = dataRange
    power = powerArray(:,j);
    save([TargetFolder '\' TargetFolder_subdir '\' dataPrefix num2str(j+dataOffset) '.mat'],'power','-append');
end;
toc; disp(['Appended power to ' num2str(length(dataRange)) ' files']);