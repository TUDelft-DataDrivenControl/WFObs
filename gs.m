clear all; close all; clc;

F_array  = 0.6:0.05:0.8;
mu_array = [0.1 1 2 5 10 25];
    
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

[minValue,minIndex]=min(mRMSE);
[minValue2,minIndex2]=min(minValue);
min(min(mRMSE));
disp(['Optimal control settings: min(mRMSE) = ' num2str(min(min(mRMSE))) ', with F=' num2str(F_array(minIndex(minIndex2))) ', Mu=' num2str(mu_array(minIndex2))]);