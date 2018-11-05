function [ time,u,v,turbData,meshingOut ] = LES_import_core( scriptOptions,rawTurbData,filterSettings,meshSetup )
%% Core code
addpath('bin');
addpath('bin/natsortfiles'); % Natural sorter library
tic;

% Load LES data (PALM: all,  SOWFA: first slice and turbine data.)
[filesInFolder,flowData,turbDataRaw,dataSource] = loadLESData(scriptOptions,rawTurbData);

% Filter and resample turbine data
if isempty(turbDataRaw)
    disp('Filtering and resampling flow data...')
    time = 0:meshSetup.dt:flowData.time(end);
else
    disp('Filtering and resampling turbine and flow data...')
    time            = 0:meshSetup.dt:min([flowData.time(end),turbDataRaw.time(end)]);
    turbData        = resampleTurbData(turbDataRaw,filterSettings.turbData,time);
end
turbData.Crx    = rawTurbData.Crx;
turbData.Cry    = rawTurbData.Cry;
turbData.Drotor = rawTurbData.Drotor;
turbData.HH     = rawTurbData.hubHeight;

% Create resample vector for flow data
k_flow = interp1(flowData.time,1:length(flowData.time),time,'nearest',1);

% Rotate and translate mesh to align wind with vertical axis
[ flowData,turbData,u_Inf,v_Inf,WD,cornerPoints ] = rotateTranslate( flowData,turbData,meshSetup );

% Determine target meshing (simplified meshing.m code)
Wp = simpleMeshing(turbData,meshSetup);

% Perform remesh for every timestep
NN = length(time);
Nt = length(turbData.Crx);
for k = 1:NN
    if k == 1
        ticIt      = tic;      % Timer for iterations (s)
        ETA        = 'N/A';    % Estimated time left
        h2         = figure;   % Create empty figure
        tri        = delaunay(flowData.yu,flowData.xu); % Necessary for plotting raw data
        k_len      = floor(log10(NN))+1;
        [u,v]      = deal(zeros(NN,meshSetup.Nx,meshSetup.Ny));
        
        % Setup interpolation functions with initial values = 0
        uInterpolant = scatteredInterpolant(flowData.yu,flowData.xu,zeros(size(flowData.yu)),'linear');
        vInterpolant = scatteredInterpolant(flowData.yv,flowData.xv,zeros(size(flowData.yu)),'linear');
        
        % Settings for exporting input settings
        CT_given = nnz(strcmp(fieldnames(turbData),'CT')) > 0 ;
        if CT_given == false
            for j = 1:Nt
                rotorPts.x{j} = ones(1,filterSettings.Ur.nPts)*turbData.Crx(j);
                rotorPts.y{j} = turbData.Cry(j)+turbData.Drotor(j)*linspace(-0.5,0.5,filterSettings.Ur.nPts);
            end
            UrMean_unf = zeros(NN,Nt);
        end
        
        % Determine .vtk file list (SOWFA only)
        if strcmp(dataSource,'sowfa')
            vtkFiles = filesInFolder(find(~cellfun(@isempty,strfind(filesInFolder,'.vtk'))));
            vtkFiles = natsortfiles(vtkFiles); % Sort numerically
        end
    else        
        elapsedTime = toc(ticIt);  % Timer and ETA calculator
        ETA = secs2timestr((NN-k+1)*(elapsedTime/(k-1)));
    end;

    % Load file (SOWFA) / flow slice (PALM)
    if strcmp(dataSource,'sowfa')
        disp(['Loading VTK file from: ''' vtkFiles{k_flow(k)} '''.'])
        [~,~,cellData] = importVTK(vtkFiles{k_flow(k)});   
        flowDataInt.uk = cellData(:,2);
        flowDataInt.vk = cellData(:,1);
    else
        flowDataInt.uk = flowData.u(k_flow(k),:)';
        flowDataInt.vk = flowData.v(k_flow(k),:)';
    end
    disp(['Performing remesh for k = ' num2str(k,['%0' num2str(k_len) 'd']) ...
        '. Progress: ' num2str(floor(100*(k-1)/NN),'%02d') '%. ETA: ' ETA '.']);
    
    % Rotated flow fields
    flowData.uk  = cos(WD)*flowDataInt.uk+sin(WD)*flowDataInt.vk;
    flowData.vk  = cos(WD)*flowDataInt.vk-sin(WD)*flowDataInt.uk;
    
    % Update values of interpolant
    uInterpolant.Values = flowData.uk;
    vInterpolant.Values = flowData.vk;
    
    % Interpolate and find velocity
    uk_remeshed = uInterpolant(Wp.ldyy(:),Wp.ldxx2(:));
    vk_remeshed = vInterpolant(Wp.ldyy2(:),Wp.ldxx(:)); 
    
    % Calculate CT_prime from flow data
    if CT_given == false && ~isempty(turbDataRaw)
        for j = 1:Nt
            UrMean_unf(k,j) = mean(uInterpolant(rotorPts.y{j},rotorPts.x{j}));
        end
    end
       
    % Produce figures
    if ~rem(k,scriptOptions.plotFrequency) | k == 1
        set(0,'CurrentFigure',h2); clf
        subplot(1,2,1);
        trisurf(tri, flowData.yu, flowData.xu, flowData.uk); % vk_raw
        lighting none; shading flat; colorbar; axis equal;
        light('Position',[-50 -15 29]); view(0,90); hold on;
        plot3([0, Wp.Ly Wp.Ly 0 0],[0 0 Wp.Lx Wp.Lx 0], [1e3*ones(5,1)],'k--' )  % Draw submesh area
        for jTurb = 1:Nt  % Draw turbines
            plot3(turbData.Cry(jTurb)+[-0.5,0.5]*turbData.Drotor(jTurb),turbData.Crx(jTurb)*[1,1],[1e3 1e3],'k-');
            text(turbData.Cry(jTurb)+0.5*turbData.Drotor(jTurb),turbData.Crx(jTurb),1e3,['T' num2str(jTurb)]);
        end
        ylim([min([-10;flowData.xu]),max([Wp.Lx+10;flowData.xu])]);
        xlim([min([-10;flowData.yu]),max([Wp.Ly+10;flowData.yu])]);
        ylabel('x-direction (m)'); xlabel('y-direction (m)');
        title('RAW $u$ [m/s]','interpreter','latex');
        
        subplot(1,2,2);
        contourf(Wp.ldyy, Wp.ldxx2,reshape(uk_remeshed,Wp.Nx,Wp.Ny),'Linecolor','none');
        axis equal tight; colorbar; hold on;
        for jTurb = 1:Nt
            plot(turbData.Cry(jTurb)+[-0.5,0.5]*turbData.Drotor(jTurb),turbData.Crx(jTurb)*[1,1],'k-');
        end
        ylabel('x-direction (m)'); xlabel('y-direction (m)');
        title('REMESHED $u$ [m/s]','interpreter','latex');
        drawnow
    end
    
    % Check if data points are even valid
    if nnz(isnan(uk_remeshed)) > 0 || nnz(isnan(vk_remeshed)) > 0
        disp('ERROR ENCOUNTERED'); disp(' ');
        disp('Your meshing is out of reach.');
        disp(['Outer (source) meshing is given by (y,x): ']);
        disp(cornerPoints)
        disp(['Inner (desired) meshing is given by (y,x): ']);
        cornerPointsInner = round([0, 0; Wp.Ly, 0; Wp.Ly, max(Wp.ldxx2(:)); 0, max(Wp.ldxx2(:))]);
        disp(cornerPointsInner);
        disp(' ')
        disp('PLOTTING COMMAND:');
        disp(['figure; plot([' num2str(cornerPoints([1:end 1],1)') '],[' num2str(cornerPoints([1:end 1],2)') ']);' ]);
        disp(['hold on; plot([' num2str(cornerPointsInner([1:end 1],1)') '],[' num2str(cornerPointsInner([1:end 1],2)') ']);' ]);
        disp(['xlabel(''y-distance (m)''); ylabel(''x-distance (m)''); axis equal;']);
        disp(['legend(''Source data'',''Target meshing''); grid on;']);
        disp(' ');
        error('Please reduce your desired grid dimensions in meshSetup.');
    end  
    
    % save to tensors
    u(k,:,:) = reshape(uk_remeshed,Wp.Nx,Wp.Ny);
    v(k,:,:) = reshape(vk_remeshed,Wp.Nx,Wp.Ny);
end

%% Post-processing
if ~isempty(turbDataRaw)
    CT_prime = zeros(NN,Nt);
    if CT_given % Ct' (Ct' = Ct/(1-a)^2, by Goit et al.)
        CT_prime_unf = turbData.CT./((1-turbData.a).^2);
        CT_prime     = CT_prime_unf;
    else
        UrMean = UrMean_unf;
        if filterSettings.Ur.MM % Filter average flow velocity on rotor
            mmLimits = ceil([filterSettings.Ur.tL filterSettings.Ur.tR]/meshSetup.dt);
            UrMean   = movmean(UrMean,mmLimits,1);
        end
        for j = 1:Nt % Determine CT_prime from turbine oop bending moment
            % Calculate raw/direct signal without any filtering
            CT_prime_unf(:,j) = 2*interp1(turbDataRaw.time,turbDataRaw.Mz(:,j),time)'./ ...
                (meshSetup.rho*UrMean_unf(:,j).^2*0.25*pi*turbData.Drotor(j)^2*turbData.HH);
            
            % Signal using filtered Umean
            CT_prime(:,j) = 2*turbData.Mz(:,j)./(meshSetup.rho*UrMean(:,j).^2 ...
                *0.25*pi*turbData.Drotor(j)^2*turbData.HH);
        end
    end
    if filterSettings.CTp.MM % Filter CT_prime directly
        mmLimits = ceil([filterSettings.CTp.tL filterSettings.CTp.tR]/meshSetup.dt);
        CT_prime = movmean(CT_prime,mmLimits,1);
    end;
    
    % Compare filtered with nonfiltered signal
    figure;
    figLayout = numSubplots(Nt);
    for j = 1:Nt
        subplot(figLayout(1),figLayout(2),j);
        plot(time,CT_prime_unf(:,j),'-','displayName','Unfiltered');
        hold on;
        plot(time,CT_prime(:,j),'--','displayName','Filtered');
        grid on;
        xlabel('Time (s)')
        ylabel(['CT'' (-)'])
        title(['Turbine ' num2str(j)])
    end
    legend('Unfiltered','Filtered')
    
    % Format inputs
    turbDataFields = fieldnames(turbData);
    for j = 1:NN
        for jField = 1:length(turbDataFields)
            jField = turbDataFields{jField};
            inputData(j).t(:,1)        = turbData.time(j);
            inputData(j).phi(:,1)      = turbData.phi(j,:);
            inputData(j).CT_prime(:,1) = CT_prime(j,:);
            inputData(j).CT_prime_unfiltered(:,1) = CT_prime_unf(j,:);
        end
    end
end

% Write properties for meshing
flow = struct(...
    'time',time,...
    'xu',
    'u',u,...
    'v',v,...
    
meshingOut.gridType  = 'lin';
meshingOut.Lx        = Wp.Lx;
meshingOut.Ly        = Wp.Ly;
meshingOut.Nx        = Wp.Nx;
meshingOut.Ny        = Wp.Ny;
meshingOut.Crx       = turbData.Crx;
meshingOut.Cry       = turbData.Cry;
meshingOut.h         = meshSetup.dt;
meshingOut.L         = time(end);
meshingOut.Drotor    = turbData.Drotor;
meshingOut.turbInput = inputData;
meshingOut.u_Inf     = u_Inf;
meshingOut.v_Inf     = v_Inf;
meshingOut.Rho       = meshSetup.rho;

toc
end