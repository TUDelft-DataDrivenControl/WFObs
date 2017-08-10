clear all; clc; close all;
addpath '..\export_fig'

%% Script settings
giveupdates = 1;        % Display progress in command line
export_ids  = 1:1999;   % Define file number(s) to export
dataoffset  = 0;        % Offset in dataset names
showmapping = 0;        % Plot mappings from SOWFA to WFSim

sourcepath   = ['K:\dcsc\DataDriven\Data\SOWFA\sampleForPODexcitationPRBSpositiveNoPrecursor\data_'];
SuperCONpath = ['K:\dcsc\DataDriven\Data\SOWFA\sampleForPODexcitationPRBSpositiveNoPrecursor\superCONOUT.csv'];
savepath     = ['..\SOWFA_export\NoPrecursor\']; % be sure to end with backslash


%% Meshing settings
distance_S  = 400  ;     % distance (m) upwind   first  turbine to export
distance_N  = 1200 ;     % distance (m) downwind  last  turbine to export
distance_W  = 700  ;     % distance (m) west most left  turbine to export
distance_E  = 700  ;     % distance (m) east most right turbine to export

meshsetup{1} = struct('name',       'lin50x50',...          % Subfoldername for saved files
                      'type',       'lin'   , ...           % Meshing type ('lin' or 'exp')
                      'Nx',         50 ,...                 % Number of points length-wise
                      'Ny',         25 );                   % Number of points width-wise

                  
%% Turbine properties (location in original axis system)
turbine.MaxCp   = 0.4866;
turbine.Drotor  = 126.3992;    
turbine.locs   = [1226.3, 1342.0; 1773.7, 1658.0];  % MS Thesis case original (x,y)
%turbine.locs    = [1118.1, 1279.5; 1881.9, 1720.5]; % Yaw case 1-3 original (x,y)
% turbine.locs    = [1086.7, 1554.9 ; ... % turbine 1 % APC 9 turbine case
%                    1634.0, 1870.9 ; ... % turbine 2
%                    1465.9, 898.1  ; ... % turbine 3
%                    1276.3, 1226.5 ; ... % turbine 4
%                    2013.2, 1214.1 ; ... % turbine 5
%                    1823.6, 1542.5 ; ... % turbine 6
%                    2560.5, 1530.1 ; ... % turbine 7
%                    2370.9, 1858.5 ; ... % turbine 8
%                    2181.3, 2186.9 ];    % turbine 9

% Please specify such that turbine 1 is upstream, with turbine 2 downstream
% of turbine 1 (aligned) in the rotated frame. The remaining turbines
% should by definition then be north and east of these 2 turbines (in rotated frame).
Nt = size(turbine.locs,1); % Number of turbines

% Create a subdirectory for each meshing resolution
if giveupdates; disp([datestr(now,'HH.MM') ': Creating subdirectory for each meshing resolution.']); end;
for j = 1:size(meshsetup,1)
    foldername{j} = [savepath meshsetup{j}.name];
    if (exist(foldername{j})~=7)
        mkdir(foldername{j});
    end;
end;

% Create a file with all defined settings
save([savepath 'export_settings.mat']);

% Include turbine data for measurements
SCO = importSuperCONOUT(SuperCONpath);

%% Loop
for i = export_ids
    tic; [dataType,cellCenters,cellData] = importVTK([sourcepath num2str(i) '.vtk']);
    
    for j = 1:size(SCO.data{1},2)
        fieldCompatibleString = strrep(SCO.sensorList(j),' ',''); % Remove spaces
        fieldCompatibleString = strrep(fieldCompatibleString,'-','_'); % Replace spaces with underscore
        kt = strfind(fieldCompatibleString{1},'(');
        if length(kt) > 0
            fieldCompatibleString = fieldCompatibleString{1}(1:kt-1); % Remove units
        else
            fieldCompatibleString = fieldCompatibleString{1};
        end;
        turb.(fieldCompatibleString) = zeros(1,Nt);
        itimeindex = find(SCO.time==i-dataoffset);
        for jt = 1:Nt
            turb.(fieldCompatibleString)(jt) = SCO.data{jt}(itimeindex,j);
        end;
    end;
    if i == export_ids(1)
        % Check if all data is on same altitude
        if nnz(cellCenters(:,3)-cellCenters(1,3))>0
            error('Data contains multiple vertical slices. Incompatible as of now.');
        end;
        
        % Calculate transformations
        if giveupdates; disp([datestr(now,'HH.MM') ': Calculating transformations for grid mapping.']); end;
        if(showmapping)
            hf = figure; subplot(1,3,1);
            plot(cellCenters(:,1),cellCenters(:,2),'k.','DisplayName','Mesh');
            title('Original SOWFA mesh'); hold on; axis equal tight;
            for j = 1:size(turbine.locs,1)
                plot(turbine.locs(j,1),turbine.locs(j,2),'*','DisplayName',['Turbine ' num2str(j)]);
            end;
            drawnow;
        end;
        
        % Calculate new grid coordinates after rotation and translation
        rotangle= pi/2 - atan((turbine.locs(2,2)-turbine.locs(1,2))/(turbine.locs(2,1)-turbine.locs(1,1)));
        x = (cellCenters(:,1)-turbine.locs(1,1))*cos(rotangle)-sin(rotangle)*(cellCenters(:,2)-turbine.locs(1,2))+distance_W;
        y = (cellCenters(:,1)-turbine.locs(1,1))*sin(rotangle)+cos(rotangle)*(cellCenters(:,2)-turbine.locs(1,2))+distance_S;
        
        % Calculate new turbine coordinates after rotation and translation
        for j = 1:Nt
            turbines_rot(j,1) = (turbine.locs(j,1)-turbine.locs(1,1))*cos(rotangle)-sin(rotangle)*(turbine.locs(j,2)-turbine.locs(1,2))+distance_W;
            turbines_rot(j,2) = (turbine.locs(j,1)-turbine.locs(1,1))*sin(rotangle)+cos(rotangle)*(turbine.locs(j,2)-turbine.locs(1,2))+distance_S;
        end;
        
        % Calculate u_inf for plot
        u = cellData(:,1)*sin(rotangle)+cos(rotangle)*cellData(:,2); % SOWFA velocity S->N
        clabelmax = ceil(max(u)); clear u;
        
        if(showmapping)
            figure(hf); subplot(1,3,2);
            plot(x,y,'k.','DisplayName','Mesh');
            title('Rotated and translated SOWFA mesh'); hold on; axis equal tight;
            for j = 1:Nt
                plot(turbines_rot(j,1),turbines_rot(j,2),'*','DisplayName',['Turbine ' num2str(j)]);
            end;
            drawnow;
        end;
        
        % Calculate new WFSim meshes
        if giveupdates; disp([datestr(now,'HH.MM') ': Calculating new meshings.']); end;
        Lx = max(turbines_rot(:,2))+distance_N;
        Ly = max(turbines_rot(:,1))+distance_E;
        
        for j = 1:length(meshsetup)
            if strcmp(lower(meshsetup{j}.type), 'lin')
                % linear gridding
                Wp{j}.ldx  = linspace(0,Lx,meshsetup{j}.Nx);
                Wp{j}.ldy  = linspace(0,Ly,meshsetup{j}.Ny);
                Wp{j}.ldxx = repmat(Wp{j}.ldx',1,meshsetup{j}.Ny);
                Wp{j}.ldyy = repmat(Wp{j}.ldy,meshsetup{j}.Nx,1);
                
                % Create secondary grid from primary grid
                Wp{j}.ldx2  = 0.5*(Wp{j}.ldx(1:end-1)+Wp{j}.ldx(2:end));
                Wp{j}.ldx2  = [Wp{j}.ldx2 2*Wp{j}.ldx2(end)-Wp{j}.ldx2(end-1)]; % add extra cells
                Wp{j}.ldy2  = 0.5*(Wp{j}.ldy(1:end-1)+Wp{j}.ldy(2:end));
                Wp{j}.ldy2  = [Wp{j}.ldy2 2*Wp{j}.ldy2(end)-Wp{j}.ldy2(end-1)]; % add extra cells
                Wp{j}.ldxx2 = repmat(Wp{j}.ldx2',1,meshsetup{j}.Ny);
                Wp{j}.ldyy2 = repmat(Wp{j}.ldy2,meshsetup{j}.Nx,1);
            else
                % Very compact form of meshing2.m function
                ldx  = expvec( Lx, meshsetup{j}.cellSize(1), uniquetol(turbines_rot(:,2),1e-2), meshsetup{j}.R_con(1), meshsetup{j}.N_con(1), meshsetup{j}.dx_min(1) );
                ldy  = expvec( Ly, meshsetup{j}.cellSize(2), uniquetol(turbines_rot(:,1),1e-2), meshsetup{j}.R_con(2), meshsetup{j}.N_con(2), meshsetup{j}.dx_min(2) );

                Wp{j}.ldxx  = repmat(ldx',1,length(ldy));
                Wp{j}.ldyy  = repmat(ldy,length(ldx),1);
                ldx2        = 0.5*(ldx(1:end-1)+ldx(2:end));
                ldx2        = [ldx2 2*ldx2(end)-ldx2(end-1)]; % add extra cells
                ldy2        = 0.5*(ldy(1:end-1)+ldy(2:end));
                ldy2        = [ldy2 2*ldy2(end)-ldy2(end-1)]; % add extra cells
                Wp{j}.ldxx2 = repmat(ldx2',1,length(ldy));
                Wp{j}.ldyy2 = repmat(ldy2,length(ldx),1);
                clear ldx ldy ldx2 ldy2
            end;
        end;
        
        if(showmapping)
            figure(hf); subplot(1,3,3);
            plot(x,y,'k.','DisplayName','SOWFA mesh'); hold on;
            plot(Wp{1}.ldyy(:),Wp{1}.ldxx(:),'g.','DisplayName','WFSim mesh');
            title('Extraction area'); hold on; axis equal tight;
            for j = 1:Nt
                plot(turbines_rot(j,1),turbines_rot(j,2),'*','DisplayName',['Turbine ' num2str(j)]);
            end;
            legend('-DynamicLegend');
            drawnow;
            
            hg = figure;
            tri = delaunay(x,y);
            
            for j = 1:size(meshsetup,1)
                triq1{j} = delaunay(Wp{j}.ldxx2,Wp{j}.ldyy );  % Needed for plotting u
                triq2{j} = delaunay(Wp{j}.ldxx, Wp{j}.ldyy2);  % Needed for plotting v
            end;
        end;
        
        % Describe how to import this new mesh to WFSims meshing.m file
        if giveupdates; disp([datestr(now,'HH.MM') ': Exporting relevant parameters for your meshing.m file.']); end;
        for j = 1:length(meshsetup)
            fid = fopen([foldername{j} '\meshing.txt'],'wt');
            fprintf(fid,['case {''' strrep(lower(foldername{j}), '\', '_') '''}\n' ...
                '   Lx     = ' num2str(Lx,8) ';\n' ...
                '   Ly     = ' num2str(Ly,8) ';\n' ...
                '   Nx     = ' num2str(meshsetup{j}.Nx) ';\n' ...
                '   Ny     = ' num2str(meshsetup{j}.Ny) ';\n' ...
                '   Crx    = ' mat2str(turbines_rot(:,2),7) ';\n' ...
                '   Cry    = ' mat2str(turbines_rot(:,1),7) ';\n' ...
                '   turbine.Drotor = ' num2str(turbine.Drotor,7) ';\n' ...
                '   turbine.MaxCp  = ' num2str(turbine.MaxCp,7) ';\n' ...
                '   Crxs   = ' mat2str(ones(size(turbine.locs,1),1)) ';   %% Obsolete option\n' ...
                '   Qx     = 1;        %% Irrelevant for linear grids\n' ...
                '   Qy     = 1;        %% Irrelevant for linear grids\n' ...
                '   sigmax = 1;        %% Irrelevant for linear grids\n' ...
                '   sigmay = 1;        %% Irrelevant for linear grids']);
            fclose(fid);
        end;
        if giveupdates; disp([datestr(now,'HH.MM') ': Note that these settings are exclusively valid for a linear meshing.']); end;
    end;
    
    % Transform velocity measurements to right coordinate system
    u = cellData(:,1)*sin(rotangle)+cos(rotangle)*cellData(:,2); % SOWFA velocity S->N
    v = cellData(:,1)*cos(rotangle)-sin(rotangle)*cellData(:,2); % SOWFA velocity W->E
    
    % Project data onto downscaled grids from WFSim
    for j = 1:size(meshsetup,1)
        uw{j} = griddata(x,y,u,Wp{j}.ldyy, Wp{j}.ldxx2, 'linear');
        vw{j} = griddata(x,y,v,Wp{j}.ldyy2, Wp{j}.ldxx, 'linear');
        
        uq = uw{j}; ux = Wp{j}.ldxx2; uy = Wp{j}.ldyy;
        vq = vw{j}; vx = Wp{j}.ldxx;  vy = Wp{j}.ldyy2;
        %save([foldername{j} '\' num2str(i)],'ux','uy','uq','vx','vy','vq','power');
        save([foldername{j} '\' num2str(i)],'ux','uy','uq','vx','vy','vq','turb');
    end;
    
    
    % Plot results in comparison
    if(showmapping)
        figure(hg); clf;
        subplot(1,size(meshsetup,1)+1,1);
        h = trisurf(tri, x, y, u); hold on;
        plot3(turbines_rot(:,1),turbines_rot(:,2),100*ones(size(turbines_rot,1),1),'r*');
        ylabel('x-direction'); xlabel('y-direction'); 
        title(['Original (i = ' num2str(i) ')']);
        lighting none; shading flat
        l = light('Position',[-50 -15 29]);
        caxis([0 clabelmax]); view(0,90); hold on;
        hold off; axis equal
        xlim([0 Ly]); ylim([0 Lx]);
        
        for j = 1:size(meshsetup,1)
            subplot(1,size(meshsetup,1)+1,j+1);
            h = trisurf(triq1{j}, Wp{j}.ldyy, Wp{j}.ldxx2, uw{j}); hold on;
            plot3(turbines_rot(:,1),turbines_rot(:,2),100*ones(size(turbines_rot,1),1),'r*');
            xlabel('y-direction'); %zlabel(['2D Velocity (i = ' num2str(i) ')']);
            title(meshsetup{j}.name);
            lighting none; shading flat
            l = light('Position',[-50 -15 29]);
            caxis([0 clabelmax]); view(0,90); hold on;
            %set(gca,'yticklabel',{[]})
            hold off; axis equal
            xlim([0 Ly]); ylim([0 Lx]);
        end;
        colorbar; clabel([],'Velocity (m/s)'); drawnow;
        
    end;
    disp([datestr(now,'HH.MM') ': Data #' num2str(i) '. Iteration time: ' num2str(toc,3) ' s']);
end;