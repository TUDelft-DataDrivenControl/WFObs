clear all; clc; close all;

%% Set-up
% Source files
scriptOptions.plotMapping = true; % Plot mapping at time k == 1  (validation)
% scriptOptions.sourcePath  = ['D:/bmdoekemeijer/My Documents/MATLAB/WFObs/WFSim/data_PALM/2turb_adm_long'];
scriptOptions.sourcePath  = ['D:/Yawexcitationcase3/firstTenSlices'];%sliceDataInstant'];

% Turbine properties directly from PALM
% (x,y)   = [1226.3, 1342.0; 1773.7, 1658.0];  % MS Thesis case original (x,y)
% (x,y)   = [1118.1, 1279.5; 1881.9, 1720.5]; % Yaw case 1-3 original (x,y)
% rawTurbData.Crx       = [5700, 6456];   % Raw (inertial frame) in (m) PALM
% rawTurbData.Cry       = [1175, 1175];   % Raw (inertial frame) in (m) PALM
rawTurbData.Cry       = [1118.1, 1881.9];   % Raw turbine locations in (m) SOWFA
rawTurbData.Crx       = [1279.5, 1720.5];   % Raw turbine locations in (m) SOWFA
% Frame is
%
%    ^
%  x |
%    |_ _ >  y
%
rawTurbData.hubHeight = 90.0;               % Hub height in (m)

% Desired output settings
meshSetup.name        = 'lin50x25';
meshSetup.dt          = 1.0 ; % Timestep in seconds
meshSetup.distance_S  = 300 ; % distance (m) upwind   first  turbine to export
meshSetup.distance_N  = 740;  % distance (m) downwind  last  turbine to export
meshSetup.distance_W  = 280 ; % distance (m) west most left  turbine to export
meshSetup.distance_E  = 280 ; % distance (m) east most right turbine to export
meshSetup.Nx          = 50;
meshSetup.Ny          = 25;


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
    [ flowDataRaw,turbDataRaw ] = loadPALMdata( filesInFolder, rawTurbData.hubHeight );
    
elseif nnz(cell2mat(strfind(filesInFolder,'.vtk'))) > 0
    % load SOWFA data
    disp('Found a .vtk file. Assuming this is SOWFA data.')
    [ flowDataRaw,turbDataRaw ] = loadSOWFAdata(filesInFolder);
else
    error('Did not find any SOWFA/PALM data in the folder.');
end

% Resample data
disp('Resampling data in time...')
time  = meshSetup.dt:meshSetup.dt:min([flowDataRaw.time(end),turbDataRaw.time(end)]);
k_raw = interp1([0; flowDataRaw.time],[1 1:length(flowDataRaw.time)],time,'nearest');
flowDataResampled      = flowDataRaw;
flowDataResampled.time = time;
for j = {'u','v'} % Resample flow fields
    flowDataResampled.(j{1}) = flowDataRaw.(j{1})(k_raw,:,:,:);
end
k_raw = interp1(turbDataRaw.time(:,1),1:length(turbDataRaw.time(:,1)),time,'nearest');
fieldNamesListTurb = fieldnames(turbDataRaw);
for j = 1:length(fieldNamesListTurb)
    turbDataResampled.(fieldNamesListTurb{j}) = turbDataRaw.(fieldNamesListTurb{j})(k_raw,:);
end
turbDataResampled.t = time;
clear j k_raw fieldNamesListTurb flowDataRaw turbDataRaw


% Normalize variables
flowData = flowDataResampled;
turbData = turbDataResampled;
clear flowDataResampled turbDataResampled

turbData.Crx    = rawTurbData.Crx;% - min(flowData.xu);
turbData.Cry    = rawTurbData.Cry;% - min(flowData.yu);
clear rawTurbData

% Determine freestream conditions flow field by checking all corners
disp('Rotating and translating grid according to u_Inf and v_Inf...')
u_Inf = median(flowData.u(:));
v_Inf = median(flowData.v(:));

% Rotate wind field
U_inf = sqrt(u_Inf^2+v_Inf^2);
WD = atan(v_Inf/u_Inf); % in radians
if abs(WD) > deg2rad(2.5) 
    [flowData,turbData] = rotateTranslate(flowData,turbData,WD);
end

% Preprocessing for translation: check if even possible
if meshSetup.distance_S > min(turbData.Crx)
    disp(' ')
    disp('WARNING: Specified output settings for x exceed available data.')
    disp('Shrinking to a smaller distance_S that matches the available data.')
    meshSetup.distance_S = min(turbData.Crx)
end
yMax  = max(flowData.yu)-min(flowData.yu);
xMax  = max(flowData.xv)-min(flowData.xv);
xTurbSeperation = max(turbData.Crx)-min(turbData.Crx);
Wp.Lx = meshSetup.distance_S + meshSetup.distance_N + xTurbSeperation;
Wp.Ly = meshSetup.distance_W + meshSetup.distance_E;
if Wp.Lx > xMax
    disp(' ')
    disp('WARNING: Specified output settings for x exceed available data.')
    disp('Resizing to a symmetric mesh that matches the available data.')
    [meshSetup.distance_S,meshSetup.distance_N] = deal((xMax - xTurbSeperation) / 2)
end
if Wp.Ly > yMax
    disp('WARNING: Specified output settings for y exceed available data.');
    disp('Resizing to a symmetric mesh that matches the available data.');
    [meshSetup.distance_W,meshSetup.distance_E] = deal(yMax / 2)
end
clear xMax yMax xTurbSeperation

% Translation
[~,UpstrIndx] = min(turbData.Crx); % Find most upstream turbine
flowData.xv = flowData.xv - turbData.Crx(UpstrIndx) + meshSetup.distance_S;
flowData.xu = flowData.xu - turbData.Crx(UpstrIndx) + meshSetup.distance_S;
flowData.yu = flowData.yu - turbData.Cry(UpstrIndx) + meshSetup.distance_W;
flowData.yv = flowData.yv - turbData.Cry(UpstrIndx) + meshSetup.distance_W;
turbData.Crx = turbData.Crx - turbData.Crx(UpstrIndx) + meshSetup.distance_S;
turbData.Cry = turbData.Cry - turbData.Cry(UpstrIndx) + meshSetup.distance_W;
clear UpstrIndx

% Determine target meshing (simplified meshing.m code)
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


%% Perform remesh for every timestep
NN      = length(time);
[u,v]   = deal(zeros(NN,meshSetup.Nx,meshSetup.Ny));
for k = 1:NN
    disp(['Performing remesh for k = ' num2str(k) '...']);
    
    % u-velocity
    uk_raw      = flowData.u(k,:);   
    uk_remeshed = griddata(flowData.yu,flowData.xu,uk_raw,Wp.ldyy(:),Wp.ldxx2(:), 'nearest');
    uk_remeshed(isnan(uk_remeshed)) = 0; % remove NaNs
    
    % v-velocity
    vk_raw      = flowData.v(k,:);  
    vk_remeshed = griddata(flowData.yv,flowData.xv,vk_raw,Wp.ldyy2(:),Wp.ldxx(:), 'nearest');
    vk_remeshed(isnan(uk_remeshed)) = 0; % remove NaNs    
    
    % save to tensors
    u(k,:,:) = reshape(uk_remeshed,Wp.Nx,Wp.Ny);
    v(k,:,:) = reshape(vk_remeshed,Wp.Nx,Wp.Ny);
        
    if scriptOptions.plotMapping
        clf;
        subplot(1,2,1);
        tri = delaunay(flowData.yu,flowData.xu);
        h = trisurf(tri, flowData.yu, flowData.xu, uk_raw);
%         h = trisurf(tri, flowData.yu, flowData.xu, vk_raw);
        lighting none; shading flat
        l = light('Position',[-50 -15 29]); view(0,90);
        caxis([min(min(uk_raw)) max(max(uk_raw))+.01]);  hold on; 
%         caxis([min(vk_raw) max(vk_raw)+.01]);  hold on; 
        axis equal; axis tight; colorbar; 
        % Draw submesh area
        plot3([0, Wp.Ly Wp.Ly 0 0],[0 0 Wp.Lx Wp.Lx 0], [1e3*ones(5,1)],'k--' ) 
        % Draw turbines
        for jTurb = 1:length(turbData.Crx)
            plot3(turbData.Cry(jTurb)+[-60,60],turbData.Crx(jTurb)*[1,1],[1e3 1e3],'k-');
        end
        ylabel('x-direction (m)'); xlabel('y-direction (m)');
        %title('RAW $u$ [m/s]','interpreter','latex');
        title('RAW $v$ [m/s]','interpreter','latex');

        subplot(1,2,2);
        tri = delaunay(Wp.ldyy2(:),Wp.ldxx(:));
        h = trisurf(tri, Wp.ldyy2(:), Wp.ldxx(:), uk_remeshed);
        lighting none; shading flat
        l = light('Position',[-50 -15 29]); view(0,90);
        caxis([min(min(uk_raw)) max(max(uk_raw))+.01]);  hold on; 
        axis equal tight; colorbar; 
        for jTurb = 1:length(turbData.Crx)
            plot3(turbData.Cry(jTurb)+[-60,60],turbData.Crx(jTurb)*[1,1],[10 10],'k-');
        end
        ylabel('x-direction (m)'); xlabel('y-direction (m)');
        title('REMESHED $u$ [m/s]','interpreter','latex');
        
        clear jTurb tri l 
        drawnow
    end
end
clear Nx_raw Ny_raw k uk_raw uk_remeshed vk_raw vk_remeshed

% % Save output data
% disp('Saving output data...');
% [fname,pth] = uiputfile('.mat');
% save([pth fname],'time','u','v','xu','yu','xv','yv','turbData','Wp')
toc