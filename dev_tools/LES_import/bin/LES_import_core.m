function [ time,u,v,turbData,meshingOut ] = LES_import_core( scriptOptions,rawTurbData,meshSetup )
%% Core code
addpath('bin');
tic;
% Sort and import files from PALM data folder
[flowData,turbData,CT_given] = loadLESData( scriptOptions,rawTurbData.hubHeight );

% Filter and resample turbine data
disp('Filtering and resampling turbine and flow data...')
time            = meshSetup.dt:meshSetup.dt:min([flowData.time(end),turbData.time(end)]);
turbData        = resampleTurbData(turbData,time);
turbData.Crx    = rawTurbData.Crx;
turbData.Cry    = rawTurbData.Cry;
turbData.Drotor = rawTurbData.Drotor;
turbData.HH     = rawTurbData.hubHeight;

% Resample flow data
k_raw = interp1(flowData.time,1:length(flowData.time),time,'nearest');
flowData.time = time;
for j = {'u','v'}
    flowData.(j{1}) = flowData.(j{1})(k_raw,:);
end
clear j k_raw

% Rotate and translate mesh to align wind with vertical axis
[ flowData,turbData,u_Inf,v_Inf,cornerPoints ] = rotateTranslate( flowData,turbData,meshSetup );

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
        cmin       = min(flowData.u(:));
        cmax       = max(flowData.u(:));
        
        % Setup interpolation functions
        uInterpolant = scatteredInterpolant(flowData.yu,flowData.xu,flowData.u(k,:)',intpMethod);
        vInterpolant = scatteredInterpolant(flowData.yv,flowData.xv,flowData.v(k,:)',intpMethod);
        
        % Settings for exporting input settings
        if CT_given == false
            rotorPts.N = 50; % number of interpolation points over rotor
            for j = 1:length(turbData.Crx)
                rotorPts.x{j} = ones(1,rotorPts.N)*turbData.Crx(j);
                rotorPts.y{j} = turbData.Cry(j)+turbData.Drotor(j)*linspace(-0.5,0.5,rotorPts.N);
            end
            CT_prime = zeros(NN,length(turbData.Crx));
        end
    else
        % Update values of interpolant
        uInterpolant.Values = flowData.u(k,:)';
        vInterpolant.Values = flowData.v(k,:)';
        
        elapsedTime = toc(ticIt);  % Timer and ETA calculator
        ETA = secs2timestr((NN-k+1)*(elapsedTime/(k-1)));
    end;
    disp(['Performing remesh for k = ' num2str(k,['%0' num2str(k_len) 'd']) ...
        '. Progress: ' num2str(floor(100*(k-1)/NN),'%02d') '%. ETA: ' ETA '.']);
    
    % Interpolate and find velocity
    uk_remeshed = uInterpolant(Wp.ldyy(:),Wp.ldxx2(:));
    vk_remeshed = vInterpolant(Wp.ldyy2(:),Wp.ldxx(:)); 
    
    % Calculate CT_prime from flow data
    if CT_given == false
        for j = 1:length(turbData.Crx)
            rotorVelocity     = uInterpolant(rotorPts.y{j},rotorPts.x{j});
            rotorVelocityMean = mean(rotorVelocity);
            CT_prime(k,j)     = 2*turbData.Mz(k,j)/(meshSetup.rho*rotorVelocityMean^2 ...
                *0.25*pi*turbData.Drotor(j)^2*turbData.HH);
        end
        clear rotorVelocity rotorVelocityMean
    end
       
    % Produce figures
    if ~rem(k,scriptOptions.plotFrequency) | k == 1
        set(0,'CurrentFigure',h2); clf
        subplot(1,2,1);
        trisurf(tri, flowData.yu, flowData.xu, flowData.u(k,:)); % vk_raw
        lighting none; shading flat; colorbar; axis equal;
        light('Position',[-50 -15 29]); view(0,90);
        caxis([cmin cmax]);  hold on;
        plot3([0, Wp.Ly Wp.Ly 0 0],[0 0 Wp.Lx Wp.Lx 0], [1e3*ones(5,1)],'k--' )  % Draw submesh area
        for jTurb = 1:length(turbData.Crx)  % Draw turbines
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
        caxis([cmin cmax]);  hold on;
        for jTurb = 1:length(turbData.Crx)
            plot(turbData.Cry(jTurb)+[-0.5,0.5]*turbData.Drotor(jTurb),turbData.Crx(jTurb)*[1,1],'k-');
        end
        ylabel('x-direction (m)'); xlabel('y-direction (m)');
        title('REMESHED $u$ [m/s]','interpreter','latex');
        drawnow
    end
    
    % Check if data points are even valid
    if nnz(isnan(uk_remeshed)) + nnz(isnan(vk_remeshed)) > 0
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