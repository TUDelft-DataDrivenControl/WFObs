clear all; close all; clc;

F_array  = 0.5:0.1:2.0;
mu_array = 0.0:0.1:2.0;
    
% Do GS simulations
mRMSE = zeros(length(F_array),length(mu_array));
parfor i = 1:length(F_array) 
    F_in = F_array(i);
    disp(['Evaluating F_in = ' num2str(F_in)]);
    mRMSE_tmp = zeros(1,length(mu_array));
    for j = 1:length(mu_array)
        mu_in = mu_array(j);
        disp(['Evaluating mu_in = ' num2str(mu_in)]);
        mRMSE_tmp(j) = WFObs_gs(mu_in,F_in);
    end;
    mRMSE(i,:) = mRMSE_tmp;
end;
save('gsOut_workspace.mat');