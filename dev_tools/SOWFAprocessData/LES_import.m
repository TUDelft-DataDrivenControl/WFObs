clear all; clc; close all;
%  %  %  %  %  %  %  %  %  %  %  %  %  %

%
%
%
%
%
%% Set-up
% Source files
scriptOptions.plotFrequency = 200; % Plot mapping every * instances (will always plot k == 1, set to high value for no plots after k == 1)
% scriptOptions.sourcePath  = ['D:/bmdoekemeijer/My Documents/MATLAB/WFObs/WFSim/data_PALM/2turb_adm_long'];
scriptOptions.sourcePath  = ['D:/Yawexcitationcase3/sliceDataInstant'];%firstTenSlices'];
scriptOptions.saveMemory  = false; % turn on if you are having memory issues (SOWFA data only)

% Turbine properties directly from PALM or SOWFA. The reference frame is 
%   x (vertical, upwards pos.) - y (horizontal, rightwards pos.).
% rawTurbData = struct('Crx',[5700.0, 6456.0],'Cry',[1175.0, 1175.0]); % PALM 2-turb case
% rawTurbData = struct('Crx',[1342.0, 1658.0],'Cry',[1226.3, 1773.7]); % SOWFA MS thesis case
rawTurbData = struct('Crx',[1279.5, 1720.5],'Cry',[1118.1, 1881.9]); % SOWFA YawCase 1-3

rawTurbData.Drotor    = [126.4, 126.4]; % Rotor diameter in (m)
rawTurbData.hubHeight = 90.0;           % Hub height in (m)
rawTurbData.tau       = 3;              % Time constant tau in low pass filter 1/(tau*s+1)

% Desired output settings
meshSetup.dt          = 1.0 ; % Timestep in seconds
meshSetup.rho         = 1.20; % Air density (kg m^-3)
meshSetup.distance_S  = 180 ; % distance (m) upwind   first  turbine to export
meshSetup.distance_N  = 800;  % distance (m) downwind  last  turbine to export
meshSetup.distance_W  = 270 ; % distance (m) west most left  turbine (from hub) to export
meshSetup.distance_E  = 270 ; % distance (m) east most right turbine (from hub) to export
meshSetup.Nx          = 50;   % Number of grid points in x-direction (-)
meshSetup.Ny          = 25;   % Number of grid points in y-direction (-)



%% Core code
addpath('bin');
tic;
% Sort and import files from PALM data folder
[flowData,turbData,CT_given] = loadLESData( scriptOptions,rawTurbData.hubHeight );

% Filter and resample turbine data
disp('Filtering and resampling turbine and flow data...')
time            = meshSetup.dt:meshSetup.dt:min([flowData.time(end),turbData.time(end)]);
turbData        = resampleTurbData(turbData,rawTurbData.tau,time);
turbData.Crx    = rawTurbData.Crx;
turbData.Cry    = rawTurbData.Cry;
turbData.Drotor = rawTurbData.Drotor;
turbData.HH     = rawTurbData.hubHeight;

% Resample flow data
k_raw = interp1([0; flowData.time],[1 1:length(flowData.time)],time,'nearest');
flowData.time = time;
for j = {'u','v'}
    flowData.(j{1}) = flowData.(j{1})(k_raw,:);
end
clear j k_raw

% Rotate and translate mesh to align wind with vertical axis
[ flowData,turbData,u_Inf,v_Inf ] = rotateTranslate( flowData,turbData,meshSetup );

% Determine target meshing (simplified meshing.m code)
Wp = simpleMeshing(turbData,meshSetup);

% Perform remesh for every timestep
NN = length(time);
for k = 1:NN
    if k == 1 
        ticIt      = tic;      % Timer for iterations (s)
        ETA        = 'N/A';    % Estimated time left
        h2         = figure;   % Create empty figure
        tri        = delaunay(flowData.yu,flowData.xu); % Necessary for plotting raw data
        intpMethod = 'linear'; % Interpolation method for flow
        k_len      = floor(log10(NN))+1;
        [u,v]      = deal(zeros(NN,meshSetup.Nx,meshSetup.Ny));
        
        % Settings for exporting input settings
        if CT_given == false
            rotorPts.N = 50; % number of interpolation points over rotor
            for j = 1:length(turbData.Crx)
                rotorPts.x{j} = ones(1,rotorPts.N)*turbData.Crx(j);  
                rotorPts.y{j} = turbData.Cry(j)+turbData.Drotor(j)*linspace(-0.5,0.5,rotorPts.N);
            end
            CT_prime   = zeros(NN,length(turbData.Crx));  
        end
    else
        elapsedTime = toc(ticIt);  % Timer and ETA calculator
        ETA = secs2timestr((NN-k+1)*(elapsedTime/(k-1))); 
    end;
    disp(['Performing remesh for k = ' num2str(k,['%0' num2str(k_len) 'd']) ...
          '. Progress: ' num2str(floor(100*(k-1)/NN),'%02d') '%. ETA: ' ETA '.']);

    % Grab velocity component at time k
    uk_raw      = flowData.u(k,:);   
    vk_raw      = flowData.v(k,:);  
    uk_remeshed = griddata(flowData.yu,flowData.xu,uk_raw,Wp.ldyy(:),Wp.ldxx2(:), intpMethod);
    vk_remeshed = griddata(flowData.yv,flowData.xv,vk_raw,Wp.ldyy2(:),Wp.ldxx(:), intpMethod);   
    
    if k == 1 % Check if data points are even valid
        if nnz(isnan(uk_remeshed)) + nnz(isnan(vk_remeshed)) > 0
            w = warndlg('The specified meshing falls outside of the source data range. Will use the nearest-neighbour approach for extrapolation.');
            drawnow; waitfor(w); intpMethod = 'nearest'; % Change interpolation method to 'nearest' to avoid NaNs
            uk_remeshed = griddata(flowData.yu,flowData.xu,uk_raw,Wp.ldyy(:),Wp.ldxx2(:), intpMethod);
            vk_remeshed = griddata(flowData.yv,flowData.xv,vk_raw,Wp.ldyy2(:),Wp.ldxx(:), intpMethod);
        end  
    end
    
    % Check if at least 'a' or/and 'CT' exist. If not: calculate from flow data
    if CT_given == false
        for j = 1:length(turbData.Crx)
            rotorVelocity     = griddata(flowData.yu,flowData.xu,uk_raw,rotorPts.y{j},rotorPts.x{j}, 'linear');
            rotorVelocityMean = mean(rotorVelocity);
            CT_prime(k,j)     = 2*turbData.Mz(k,j)/(meshSetup.rho*rotorVelocityMean^2 ...
                                 *0.25*pi*turbData.Drotor(j)^2*turbData.HH);
        end
        clear rotorVelocity rotorVelocityMean        
    end
    
    % save to tensors
    u(k,:,:) = reshape(uk_remeshed,Wp.Nx,Wp.Ny);
    v(k,:,:) = reshape(vk_remeshed,Wp.Nx,Wp.Ny);

    % Produce figures
    if ~rem(k,scriptOptions.plotFrequency) | k == 1
        set(0,'CurrentFigure',h2); clf
        subplot(1,2,1);
        trisurf(tri, flowData.yu, flowData.xu, uk_raw); % vk_raw
        lighting none; shading flat; colorbar; axis equal;
        light('Position',[-50 -15 29]); view(0,90);
        caxis([min(min(uk_raw)) max(max(uk_raw))+.01]);  hold on;
        plot3([0, Wp.Ly Wp.Ly 0 0],[0 0 Wp.Lx Wp.Lx 0], [1e3*ones(5,1)],'k--' )  % Draw submesh area
        for jTurb = 1:length(turbData.Crx)  % Draw turbines
            plot3(turbData.Cry(jTurb)+[-0.5,0.5]*turbData.Drotor(jTurb),turbData.Crx(jTurb)*[1,1],[1e3 1e3],'k-');
        end
        ylim([min([-10;flowData.xu]),max([Wp.Lx+10;flowData.xu])]);
        xlim([min([-10;flowData.yu]),max([Wp.Ly+10;flowData.yu])]);        
        ylabel('x-direction (m)'); xlabel('y-direction (m)');
        title('RAW $u$ [m/s]','interpreter','latex');
        
        subplot(1,2,2);
        contourf(Wp.ldyy, Wp.ldxx2,reshape(uk_remeshed,Wp.Nx,Wp.Ny),'Linecolor','none'); 
        axis equal tight; colorbar; hold on;
        caxis([min(min(uk_raw)) max(max(uk_raw))+.01]); 
        for jTurb = 1:length(turbData.Crx)
            plot(turbData.Cry(jTurb)+[-0.5,0.5]*turbData.Drotor(jTurb),turbData.Crx(jTurb)*[1,1],'k-');
        end
        ylabel('x-direction (m)'); xlabel('y-direction (m)');
        title('REMESHED $u$ [m/s]','interpreter','latex');
        drawnow
    end
end

if CT_given
    % Determine Ct' (Ct' = Ct/(1-a)^2, by Goit et al.)
    CT_prime = turbData.CT./((1-turbData.a).^2);    
end

% % Format inputs
turbDataFields = fieldnames(turbData);
for j = 1:NN
    for jField = 1:length(turbDataFields)
        jField = turbDataFields{jField};
        inputData(j).t(:,1)        = turbData.time(j);
        inputData(j).phi(:,1)      = turbData.phi(j,:);
        inputData(j).CT_prime(:,1) = CT_prime(j,:);
    end
end

% Write properties for meshing
meshingOut.type   = 'lin';
meshingOut.Lx     = Wp.Lx;
meshingOut.Ly     = Wp.Ly;
meshingOut.Nx     = Wp.Nx;
meshingOut.Ny     = Wp.Ny;
meshingOut.Crx    = turbData.Crx;
meshingOut.Cry    = turbData.Cry;
meshingOut.h      = meshSetup.dt;
meshingOut.L      = time(end);
meshingOut.Drotor = turbData.Drotor;
meshingOut.input  = inputData;
meshingOut.u_Inf  = u_Inf;
meshingOut.v_Inf  = v_Inf;
meshingOut.Rho    = meshSetup.rho;

clear turbDataFields Nx_raw Ny_raw k uk_raw uk_remeshed ...
      vk_raw vk_remeshed jTurb l h tri triremeshed j jField

  
%% Save output data
[~,namePH,~]=fileparts(scriptOptions.sourcePath); % Placeholder filename
namePH = [namePH '_' num2str(meshSetup.Nx) 'x' num2str(meshSetup.Ny) '_data.mat'];
disp('Saving flow and turbine data...');
[fname1,pth] = uiputfile(namePH,'Save flow data');
save([pth fname1],'time','u','v','turbData')
if length(strfind(fname1,'data.mat')) > 0
    namePH_2 = strrep(fname1,'data.mat','meshing.mat'); 
else
    namePH_2 = [namePH '_meshing.mat']; 
end
disp('Saving meshing information...');
[fname2,pth] = uiputfile(namePH_2,'Save meshing data');
save([pth fname2],'-struct','meshingOut')
toc