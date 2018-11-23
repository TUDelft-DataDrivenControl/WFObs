clear all; close all; clc;
addpath('../../WFSim/bin/core');
addpath('../../bin_supplementary/offline_tools');

sourceFolder = 'GS_out_palm6turb';
outputFile = 'GS_postProcessed_palm6turb.mat';
LESDataFile = '../../data_LES/LESData_palm_6turb_adm_turbl.mat'; % Specify LES data file

% List all files
filesSource = dir(sourceFolder);
if length(filesSource) < 3
    error('Please specify a valid directory.');
end
for i = 3:length(filesSource)
    fileList{i-2} = [filesSource(i).folder '/' filesSource(i).name];
end
clear i 

% Post-processing
disp('Setting up LES data...')
load(fileList{1}); % Load first file to get WpPar/meshing
Wp = meshing(WpPar,0,0);
LESData = loadLESdata(LESDataFile); % Load pregenerated LES data
for i = 1:length(sol_array)
    solTrue(i) = struct(...
        'u',reshape(...
                LESData.uInterpolant(...
                    sol_array(i).time*ones(size(Wp.mesh.ldxx2(:))),...
                    Wp.mesh.ldxx2(:),Wp.mesh.ldyy(:)),...
                Wp.mesh.Nx,Wp.mesh.Ny), ...
        'v',reshape(...
                LESData.vInterpolant(sol_array(i).time*ones(size(Wp.mesh.ldxx(:))),...
                	Wp.mesh.ldxx(:),Wp.mesh.ldyy2(:)),...
                Wp.mesh.Nx,Wp.mesh.Ny), ...
        'P',LESData.PwrInterpolant(sol_array(i).time)' );
end
clear WpPar sol_array i

% Determine centerline locations
disp('Setting up centerline calculations...')
yTurbs = sort(Wp.turbine.Cry);
threshold = 0.5*Wp.turbine.Drotor;
thresholdRel = threshold/max(Wp.turbine.Cry);
yTurbsUnique = uniquetol(yTurbs,thresholdRel);
for i = 1:length(yTurbsUnique)
    xCL = Wp.mesh.ldxx2(1):5:Wp.mesh.ldxx(end);
    yCL = mean(yTurbs(abs((yTurbs-yTurbsUnique(i)))<threshold));
    yCL = yCL + [-Wp.turbine.Drotor:5:Wp.turbine.Drotor];
    [X,Y] = ndgrid(xCL,yCL);
    centerline(i) = struct('X',X,'Y',Y,'U_LES',@(t) mean(LESData.uInterpolant(t*ones(size(X)),X,Y),2));
end
% Create flow interpolant of WFSim
flowInterpolant = griddedInterpolant(Wp.mesh.ldxx2,Wp.mesh.ldyy,zeros(size(Wp.mesh.ldxx2)));
clear i Wp X Y yTurbs threshold thresholdRel yTurbsUnique xCL yCL LESData

% Process each GS output file
disp('Starting file processing.')
parpool(40)
parfor j = 1:length(fileList)
    loadedData = load(fileList{j});
    WpPar     = loadedData.WpPar;
    sol_array = loadedData.sol_array;
    
    % Determine flow error
    flowError = abs([vec([sol_array.uEst])-vec([solTrue.u]);...
                     vec([sol_array.vEst])-vec([solTrue.v])]);
    flowError(isnan(flowError))=[]; % remove NaN entries
    
    % Find optimal power scaling
    powerscaleOpt = mean(mean([solTrue.P]./[sol_array.PEst]));
    for i = 1:length(sol_array)
        sol_array(i).PEst = sol_array(i).PEst * powerscaleOpt;
    end
    
    % Determine power error
    powerLES = [solTrue.P];
    powerWFSim = [sol_array.PEst];    
    powerError = abs(powerWFSim-powerLES);
    powerVAF = zeros(1,length(WpPar.turbine.Crx));
    for i = 1:length(WpPar.turbine.Crx) % = number of turbines
        powerVAF(i) = vaf(powerLES(i,:),powerWFSim(i,:)); % var. accounted for (%)
    end
    
    % Determine centerline error
    [RMSE_clines,VAF_clines] = deal(zeros(length(sol_array),length(centerline)));
    for ii = 1:length(sol_array)
        tmp = ii;
        flowInterpolantTmp = flowInterpolant;
        flowInterpolantTmp.Values = sol_array(ii).uEst;
        for ic = 1:length(centerline)
            cline_WFSim = mean(flowInterpolantTmp(centerline(ic).X,centerline(ic).Y),2);
            cline_LES   = centerline(ic).U_LES(sol_array(ii).time);
            RMSE_clines(ii,ic) = sqrt(mean(cline_WFSim-cline_LES).^2);
            VAF_clines(ii,ic) = vaf(cline_LES,cline_WFSim);
        end
    end    
    
    % Write to score output struct
    scoreOut(j).mRMSE_flow  = rms(flowError);
    scoreOut(j).mRMSE_cline = mean(RMSE_clines,1); % time-averaged mRMSE
    scoreOut(j).mVAF_cline  = mean(VAF_clines,1); % time-averaged VAF
    scoreOut(j).forcescale  = WpPar.turbine.forcescale;
    scoreOut(j).powerscale  = powerscaleOpt;
    scoreOut(j).mRMSE_power = mean(powerError,2)';
    scoreOut(j).mmRMSE_power = mean((mean(powerError,2)));
    scoreOut(j).VAF_power   = powerVAF;
    scoreOut(j).mVAF_power  = mean(powerVAF);
    scoreOut(j).Wp          = WpPar;
end
clear j 

save(outputFile,'scoreOut');
disp('Finished.')