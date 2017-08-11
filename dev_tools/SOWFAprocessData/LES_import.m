clear all; clc; close all;

addpath bin

%% Set-up
% Source files
scriptOptions.plotFrequency = 5; % Plot mapping every * instances (will always plot k == 1, set to high value for no plots after k == 1)
% scriptOptions.sourcePath  = ['D:/bmdoekemeijer/My Documents/MATLAB/WFObs/WFSim/data_PALM/2turb_adm_long'];
scriptOptions.sourcePath  = ['D:/Yawexcitationcase3/sliceDataInstant'];%firstTenSlices'];
scriptOptions.saveMemory  = true; % turn on if you are having memory issues

% Turbine properties directly from PALM
% (x,y)   = [1226.3, 1342.0; 1773.7, 1658.0];  % MS Thesis case original (x,y)
% (x,y)   = [1118.1, 1279.5; 1881.9, 1720.5]; % Yaw case 1-3 original (x,y)
% rawTurbData.Crx       = [5700, 6456];   % Raw (inertial frame) in (m) PALM
% rawTurbData.Cry       = [1175, 1175];   % Raw (inertial frame) in (m) PALM
rawTurbData.Drotor    = [126.4, 126.4];
rawTurbData.Crx       = [1279.5, 1720.5];   % Raw turbine locations in (m) SOWFA
rawTurbData.Cry       = [1118.1, 1881.9];   % Raw turbine locations in (m) SOWFA
rawTurbData.hubHeight = 90.0;               % Hub height in (m)
rawTurbData.tau       = 3; % Time constant tau in low pass filter 1/(tau*s+1)
% Frame is
%
%    ^
%  x |
%    |_ _ >  y
%


% Desired output settings
meshSetup.name        = 'lin50x25';
meshSetup.dt          = 1.0 ; % Timestep in seconds
meshSetup.distance_S  = 300 ; % distance (m) upwind   first  turbine to export
meshSetup.distance_N  = 740;  % distance (m) downwind  last  turbine to export
meshSetup.distance_W  = 280 ; % distance (m) west most left  turbine (from hub) to export
meshSetup.distance_E  = 280 ; % distance (m) east most right turbine (from hub) to export
meshSetup.Nx          = 50;   % Number of grid points in x-direction (-)
meshSetup.Ny          = 25;   % Number of grid points in y-direction (-)



%% Core code
tic

% Sort and import files from PALM data folder
disp('Sorting and importing files from source folder...')
addpath(scriptOptions.sourcePath); % Add source folder to directory
filesInFolder = dir(scriptOptions.sourcePath); % Gather all files from PALM/SOWFA folder
filesInFolder = {filesInFolder(3:end).name};   % Remove '.' and '..'

% Decide whether this is SOWFA or PALM data
if nnz(cell2mat(strfind(filesInFolder,'.nc'))) > 0
    disp('Found a .nc file. Assuming this is PALM data.')
    addpath(genpath('mexcdf'));
    [ flowData,turbData ] = loadPALMdata( filesInFolder, rawTurbData.hubHeight );
    CT_given              = true;
    
elseif nnz(cell2mat(strfind(filesInFolder,'.vtk'))) > 0
    % load SOWFA data
    disp('Found a .vtk file. Assuming this is SOWFA data.')
    [ flowData,turbData ] = loadSOWFAdata(filesInFolder,scriptOptions.saveMemory);
    CT_given              = false;
else
    error('Did not find any SOWFA/PALM data in the folder.');
end

% Low-pass filter and resample turbine data
disp('Filtering and resampling turbine and flow data...')
time         = meshSetup.dt:meshSetup.dt:min([flowData.time(end),turbData.time(end)]);
turbData     = resampleTurbData(turbData,rawTurbData.tau,time);
turbData.Crx    = rawTurbData.Crx;
turbData.Cry    = rawTurbData.Cry;
turbData.Drotor = rawTurbData.Drotor;
turbData.HH     = rawTurbData.hubHeight;

k_raw = interp1([0; flowData.time],[1 1:length(flowData.time)],time,'nearest');
flowData.time = time;
for j = {'u','v'}
    flowData.(j{1}) = flowData.(j{1})(k_raw,:);
end
clear j k_raw rawTurbData

% Determine freestream conditions flow field by checking all corners
disp('Rotating and translating grid according to u_Inf and v_Inf...')
u_Inf = median(flowData.u(:));
v_Inf = median(flowData.v(:));
WD = atan(v_Inf/u_Inf); % Wind direction [rad]
if abs(WD) > deg2rad(2.5) % Only rotate if mismatch > 2.5 degrees 
    [flowData,turbData] = rotateTranslate(flowData,turbData,WD);
end

% Translation
[~,UpstrIndx] = min(turbData.Crx); % Find most upstream turbine
flowData.xv = flowData.xv   - turbData.Crx(UpstrIndx) + meshSetup.distance_S;
flowData.xu = flowData.xu   - turbData.Crx(UpstrIndx) + meshSetup.distance_S;
flowData.yu = flowData.yu   - turbData.Cry(UpstrIndx) + meshSetup.distance_W;
flowData.yv = flowData.yv   - turbData.Cry(UpstrIndx) + meshSetup.distance_W;
turbData.Crx = turbData.Crx - turbData.Crx(UpstrIndx) + meshSetup.distance_S;
turbData.Cry = turbData.Cry - turbData.Cry(UpstrIndx) + meshSetup.distance_W;
clear UpstrIndx

% Determine target meshing (simplified meshing.m code)
Wp.Lx    = max(turbData.Crx) + meshSetup.distance_N;
Wp.Ly    = max(turbData.Cry) + meshSetup.distance_E;
Wp.ldx   = linspace(0,Wp.Lx,meshSetup.Nx);
Wp.ldy   = linspace(0,Wp.Ly,meshSetup.Ny);
Wp.ldxx  = repmat(Wp.ldx',1,meshSetup.Ny);
Wp.ldyy  = repmat(Wp.ldy,meshSetup.Nx,1);
Wp.ldx2  = 0.5*(Wp.ldx(1:end-1)+Wp.ldx(2:end));
Wp.ldx2  = [Wp.ldx2 2*Wp.ldx2(end)-Wp.ldx2(end-1)]; 
Wp.ldy2  = 0.5*(Wp.ldy(1:end-1)+Wp.ldy(2:end));
Wp.ldy2  = [Wp.ldy2 2*Wp.ldy2(end)-Wp.ldy2(end-1)];
Wp.ldxx2 = repmat(Wp.ldx2',1,meshSetup.Ny);
Wp.ldyy2 = repmat(Wp.ldy2,meshSetup.Nx,1);
Wp.Nx    = meshSetup.Nx;
Wp.Ny    = meshSetup.Ny;

% Perform remesh for every timestep
NN         = length(time);
[u,v]      = deal(zeros(NN,meshSetup.Nx,meshSetup.Ny));
intpMethod = 'linear';
figure;
for k = 1:NN
    disp(['Performing remesh for k = ' num2str(k) '...']);
    
    % Grab velocity component at time k
    uk_raw      = flowData.u(k,:);   
    vk_raw      = flowData.v(k,:);  
    uk_remeshed = griddata(flowData.yu,flowData.xu,uk_raw,Wp.ldyy(:),Wp.ldxx2(:), intpMethod);
    vk_remeshed = griddata(flowData.yv,flowData.xv,vk_raw,Wp.ldyy2(:),Wp.ldxx(:), intpMethod);   
    
    if k == 1 % Check if data is available 
        if nnz(isnan(uk_remeshed)) + nnz(isnan(vk_remeshed)) > 0
            w = warndlg('The specified meshing falls outside of the source data range. Will use the nearest-neighbour approach for extrapolation.');
            drawnow; waitfor(w); intpMethod = 'nearest';
            uk_remeshed = griddata(flowData.yu,flowData.xu,uk_raw,Wp.ldyy(:),Wp.ldxx2(:), intpMethod);
            vk_remeshed = griddata(flowData.yv,flowData.xv,vk_raw,Wp.ldyy2(:),Wp.ldxx(:), intpMethod);   
        end  
        % Delaunay calculations for plotting
        tri         = delaunay(flowData.yu,flowData.xu);
        triremeshed = delaunay(Wp.ldyy2(:),Wp.ldxx(:));
        
        % Settings for exporting input settings
        if CT_given == false
            rho        = 1.20;
            rotorPts.N = 50; % number of interpolation points over rotor
            rotorPts.y = turbData.Cry(j)+turbData.Drotor(j)*linspace(-0.5,0.5,rotorPts.N);
            rotorPts.x = ones(1,rotorPts.N)*turbData.Crx(j);                      
            CT_prime   = zeros(NN,length(turbData.Crx));  
        end
    end
    
    % Check if at least 'a' or/and 'CT' exist. If not: calculate from flow data
    if CT_given == false
        for j = 1:length(turbData.Crx)
            rotorVelocity = griddata(flowData.yu,flowData.xu,uk_raw,rotorPts.y,rotorPts.x, 'linear');
            rotorVelocityMean = mean(rotorVelocity);
            CT_prime(k,j) = 2*turbData.Mz(k,j)/(rho*rotorVelocityMean^2 ...
                             *0.25*pi*turbData.Drotor(j)^2*turbData.HH);
        end
        clear rotorVelocity rotorVelocityMean        
    end
    
    % save to tensors
    u(k,:,:) = reshape(uk_remeshed,Wp.Nx,Wp.Ny);
    v(k,:,:) = reshape(vk_remeshed,Wp.Nx,Wp.Ny);

    % Produce figures
    if ~rem(k,scriptOptions.plotFrequency) | k == 1
        clf; subplot(1,2,1);
        h = trisurf(tri, flowData.yu, flowData.xu, uk_raw); % vk_raw
        lighting none; shading flat; colorbar; axis equal;
        l = light('Position',[-50 -15 29]); view(0,90);
        caxis([min(min(uk_raw)) max(max(uk_raw))+.01]);  hold on;
        plot3([0, Wp.Ly Wp.Ly 0 0],[0 0 Wp.Lx Wp.Lx 0], [1e3*ones(5,1)],'k--' )  % Draw submesh area
        for jTurb = 1:length(turbData.Crx)  % Draw turbines
            plot3(turbData.Cry(jTurb)+[-60,60],turbData.Crx(jTurb)*[1,1],[1e3 1e3],'k-');
        end
        ylim([min([-10;flowData.xu]),max([Wp.Lx+10;flowData.xu])]);
        xlim([min([-10;flowData.yu]),max([Wp.Ly+10;flowData.yu])]);        
        ylabel('x-direction (m)'); xlabel('y-direction (m)');
        title('RAW $u$ [m/s]','interpreter','latex');
        
        subplot(1,2,2);
        h = trisurf(triremeshed, Wp.ldyy2(:), Wp.ldxx(:), uk_remeshed);
        lighting none; shading flat; axis equal tight; colorbar;
        l = light('Position',[-50 -15 29]); view(0,90);
        caxis([min(min(uk_raw)) max(max(uk_raw))+.01]);  hold on;
        for jTurb = 1:length(turbData.Crx)
            plot3(turbData.Cry(jTurb)+[-60,60],turbData.Crx(jTurb)*[1,1],[10 10],'k-');
        end
        ylabel('x-direction (m)'); xlabel('y-direction (m)');
        title('REMESHED $u$ [m/s]','interpreter','latex');
        drawnow
    end
end
clear Nx_raw Ny_raw k uk_raw uk_remeshed vk_raw vk_remeshed jTurb l h

if CT_given
    % Determine Ct' (Ct' = Ct/(1-a)^2, by Goit et al.)
    turbData.CT_prime = turbData.CT./((1-turbData.a).^2);    
else
    % copy from what was calculated in iterations
    turbData.CT_prime = CT_prime;
    clear CT_prime
end

% % Save output data
disp('Saving output data...');
[fname,pth] = uiputfile('.mat');
save([pth fname],'time','u','v','xu','yu','xv','yv','turbData','Wp')
toc